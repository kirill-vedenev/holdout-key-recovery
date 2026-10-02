// Exact GF(2) elimination. No floating-point arithmetic, FFT, normal equations,
// or assumption about rank is used to accept a kernel element.
#include <cuda_runtime.h>
#include "engine.hpp"

namespace tii {
#define CUDA_TRY(call) do { auto status_=(call); if(status_!=cudaSuccess) \
    throw std::runtime_error(std::string(#call)+": "+cudaGetErrorString(status_)); } while(false)
__constant__ u64 device_binom[(saarinen::MAX_K+1)*(saarinen::MAX_D+1)];
__device__ inline u64 colex_device(u64 mask) {
    u64 rank=0;int i=1;
    while(mask) { int bit=__ffsll(static_cast<long long>(mask))-1; mask&=mask-1;
        rank+=device_binom[bit*(saarinen::MAX_D+1)+i++]; }
    return rank;
}
inline unsigned blocks_for(u64 n) { return unsigned(std::max<u64>(1,std::min<u64>(65535,(n+255)/256))); }
__global__ void assemble_kernel(u64* a,const Row* rows,u64 count,u64 words,int degree) {
    for(u64 row=u64(blockIdx.x)*blockDim.x+threadIdx.x;row<count;row+=u64(gridDim.x)*blockDim.x) {
        const u64 derivative=rows[row].derivative;
        u64 allowed=rows[row].point&~derivative;
        const int take=degree-__popcll(derivative);
        int n=__popcll(allowed); if(n<take) continue;
        unsigned char pos[63];int j=0;
        while(allowed) { pos[j++]=__ffsll(static_cast<long long>(allowed))-1;allowed&=allowed-1; }
        u64 subset=(u64(1)<<take)-1,limit=u64(1)<<n;
        while(subset<limit) {
            u64 v=subset,mon=derivative;
            while(v) { mon|=u64(1)<<pos[__ffsll(static_cast<long long>(v))-1];v&=v-1; }
            u64 c=colex_device(mon);a[row*words+c/64]|=u64(1)<<(c%64);
            u64 low=subset&(~subset+1),next=subset+low;
            if(next>=limit) break;
            subset=next|(((next^subset)>>2)/low);
        }
    }
}
__global__ void gather_kernel(const u64* a,u64 rows,u64 words,u64 column,u64* out) {
    for(u64 i=u64(blockIdx.x)*blockDim.x+threadIdx.x;i<rows;i+=u64(gridDim.x)*blockDim.x)
        out[i]=a[i*words+column];
}
struct PanelBits { int count; unsigned char bit[8]; };
__global__ void panel_index_kernel(const u64* a,u64 rows,u64 words,u64 start,PanelBits p,unsigned char* indices) {
    for(u64 row=u64(blockIdx.x)*blockDim.x+threadIdx.x;row<rows;row+=u64(gridDim.x)*blockDim.x) {
        u64 v=a[row*words+start];unsigned index=0;
        for(int i=0;i<p.count;++i) index|=unsigned((v>>p.bit[i])&1)<<i;
        indices[row]=index;
    }
}
__global__ void table_kernel(const u64* panel,u64 words,u64 start,u64 tail,u64 entries,u64* table) {
    for(u64 i=u64(blockIdx.x)*blockDim.x+threadIdx.x;i<entries*tail;i+=u64(gridDim.x)*blockDim.x) {
        u64 mask=i/tail,w=i%tail,value=0;
        while(mask) {int bit=__ffsll(static_cast<long long>(mask))-1;mask&=mask-1;value^=panel[u64(bit)*words+start+w];}
        table[i]=value;
    }
}
__global__ void eliminate_kernel(u64* a,u64 rows,u64 words,u64 start,u64 tail,
                                 const unsigned char* indices,const u64* table) {
    for(u64 i=u64(blockIdx.x)*blockDim.x+threadIdx.x;i<rows*tail;i+=u64(gridDim.x)*blockDim.x) {
        u64 row=i/tail,w=i%tail;unsigned index=indices[row];
        if(index) a[row*words+start+w]^=table[u64(index)*tail+w];
    }
}
__global__ void extract_kernel(const u64* a,u64 rows,u64 words,const u64* free_mask,
                               const u64* free_values,const u32* row_pivot,u64 sw,u64* result) {
    const unsigned lane=threadIdx.x%32,warp=threadIdx.x/32,warps=blockDim.x/32;
    for(u64 row=blockIdx.x;row<rows;row+=gridDim.x) {
        u32 pivot=row_pivot[row];if(pivot==NONE) continue;
        for(u64 j=warp;j<sw;j+=warps) {
            u64 value=0;
            for(u64 w=lane;w<words;w+=32) {
                u64 bits=a[row*words+w]&free_mask[w];
                while(bits) {u64 c=w*64+__ffsll(static_cast<long long>(bits))-1;bits&=bits-1;value^=free_values[c*sw+j];}
            }
            for(int offset=16;offset;offset>>=1)
                value^=static_cast<u64>(__shfl_xor_sync(0xffffffffU,static_cast<unsigned long long>(value),offset));
            if(lane==0) result[u64(pivot)*sw+j]=value;
        }
    }
}
struct Shard {
    int id;u64 first,count;
    u64 *matrix=nullptr,*column=nullptr,*panel=nullptr,*table=nullptr;
    unsigned char* indices=nullptr;
    Row* rows=nullptr;
    explicit Shard(int device,u64 begin,u64 n):id(device),first(begin),count(n) {}
    ~Shard() {
        cudaSetDevice(id);cudaFree(matrix);cudaFree(column);cudaFree(panel);cudaFree(table);cudaFree(indices);cudaFree(rows);
    }
};
struct CudaBackend final:Backend {
    std::vector<std::unique_ptr<Shard>> shards;
    explicit CudaBackend(const Request& request,const Options& opt):Backend(request) {
        int available=0;CUDA_TRY(cudaGetDeviceCount(&available));
        for(std::size_t j=0;j<opt.devices.size();++j) {
            int id=opt.devices[j];require(id<available,"requested GPU is not visible");
            CUDA_TRY(cudaSetDevice(id));cudaDeviceProp prop{};CUDA_TRY(cudaGetDeviceProperties(&prop,id));
            require(prop.major>=6,"GPU compute capability is unsupported");
            u64 first=r.rows.size()*j/opt.devices.size(),end=r.rows.size()*(j+1)/opt.devices.size();
            if(first==end) continue;
            auto s=std::make_unique<Shard>(id,first,end-first);
            std::size_t free=0,total=0;CUDA_TRY(cudaMemGetInfo(&free,&total));
            u64 needed=s->count*r.words*8+s->count*(sizeof(Row)+9)+(256+8)*r.words*8;
            require(needed+(u64(1)<<30)<=free,"insufficient free GPU memory for matrix, tables and reserve");
            CUDA_TRY(cudaMemcpyToSymbol(device_binom,saarinen::HOST_CHOOSE.values.data(),sizeof(saarinen::HOST_CHOOSE.values)));
            CUDA_TRY(cudaMalloc(reinterpret_cast<void**>(&s->matrix),s->count*r.words*8));
            CUDA_TRY(cudaMalloc(reinterpret_cast<void**>(&s->column),s->count*8));
            CUDA_TRY(cudaMalloc(reinterpret_cast<void**>(&s->indices),s->count));
            CUDA_TRY(cudaMalloc(reinterpret_cast<void**>(&s->panel),8*r.words*8));
            CUDA_TRY(cudaMalloc(reinterpret_cast<void**>(&s->table),256*r.words*8));
            CUDA_TRY(cudaMalloc(reinterpret_cast<void**>(&s->rows),s->count*sizeof(Row)));
            CUDA_TRY(cudaMemcpy(s->rows,r.rows.data()+first,s->count*sizeof(Row),cudaMemcpyHostToDevice));
            std::cout<<"{\"event\":\"gpu_allocation\",\"device\":"<<id<<",\"name\":\""<<escape(prop.name)
                     <<"\",\"matrix_bytes\":"<<s->count*r.words*8<<",\"planned_base_bytes\":"<<needed
                     <<",\"free_bytes_before\":"<<free<<"}"<<std::endl;
            shards.push_back(std::move(s));
        }
    }
    std::string name() const override {return "cuda_exact_panel_elimination_"+std::to_string(shards.size())+"_devices";}
    void sync() {
        for(auto& s:shards) {CUDA_TRY(cudaSetDevice(s->id));CUDA_TRY(cudaDeviceSynchronize());}
    }
    void initialize() override {
        for(auto& s:shards) {
            CUDA_TRY(cudaSetDevice(s->id));CUDA_TRY(cudaMemsetAsync(s->matrix,0,s->count*r.words*8));
            assemble_kernel<<<blocks_for(s->count),256>>>(s->matrix,s->rows,s->count,r.words,r.d);
            CUDA_TRY(cudaGetLastError());
        }
        sync();
    }
    void gather(u64 word,std::vector<u64>& columns) override {
        columns.resize(r.rows.size());
        for(auto& s:shards) {
            CUDA_TRY(cudaSetDevice(s->id));
            gather_kernel<<<blocks_for(s->count),256>>>(s->matrix,s->count,r.words,word,s->column);
            CUDA_TRY(cudaGetLastError());
        }
        for(auto& s:shards) {
            CUDA_TRY(cudaSetDevice(s->id));
            CUDA_TRY(cudaMemcpy(columns.data()+s->first,s->column,s->count*8,cudaMemcpyDeviceToHost));
        }
    }
    void read_rows(u64 first,u64 count,u64* dst) override {
        for(auto& s:shards) {
            u64 lo=std::max(first,s->first),hi=std::min(first+count,s->first+s->count);if(lo>=hi) continue;
            CUDA_TRY(cudaSetDevice(s->id));
            CUDA_TRY(cudaMemcpy(dst+(lo-first)*r.words,s->matrix+(lo-s->first)*r.words,(hi-lo)*r.words*8,cudaMemcpyDeviceToHost));
        }
    }
    void write_rows(u64 first,u64 count,const u64* src) override {
        for(auto& s:shards) {
            u64 lo=std::max(first,s->first),hi=std::min(first+count,s->first+s->count);if(lo>=hi) continue;
            CUDA_TRY(cudaSetDevice(s->id));
            CUDA_TRY(cudaMemcpy(s->matrix+(lo-s->first)*r.words,src+(lo-first)*r.words,(hi-lo)*r.words*8,cudaMemcpyHostToDevice));
        }
    }
    void eliminate(u64 begin,const std::vector<u32>& pivots,const std::vector<u64>& panel) override {
        u64 start=begin/64,tail=r.words-start,entries=u64(1)<<pivots.size();
        PanelBits bits{};bits.count=pivots.size();for(int j=0;j<bits.count;++j) bits.bit[j]=pivots[j]%64;
        for(auto& s:shards) {
            CUDA_TRY(cudaSetDevice(s->id));
            CUDA_TRY(cudaMemcpy(s->panel,panel.data(),panel.size()*8,cudaMemcpyHostToDevice));
            panel_index_kernel<<<blocks_for(s->count),256>>>(s->matrix,s->count,r.words,start,bits,s->indices);
            CUDA_TRY(cudaGetLastError());
            table_kernel<<<blocks_for(entries*tail),256>>>(s->panel,r.words,start,tail,entries,s->table);
            CUDA_TRY(cudaGetLastError());
            eliminate_kernel<<<blocks_for(s->count*tail),256>>>(s->matrix,s->count,r.words,start,tail,s->indices,s->table);
            CUDA_TRY(cudaGetLastError());
        }
        sync();
    }
    std::vector<u64> extract(const std::vector<u64>& values,const std::vector<u64>& free_mask,
                             const std::vector<u32>& row_pivot,u32 count) override {
        const u64 sw=(u64(count)+63)/64,bytes=r.N*sw*8;
        struct Temporary {
            int id;u64 *values=nullptr,*mask=nullptr,*out=nullptr;u32* pivots=nullptr;
            ~Temporary(){cudaSetDevice(id);cudaFree(values);cudaFree(mask);cudaFree(out);cudaFree(pivots);}
        };
        std::vector<std::unique_ptr<Temporary>> tmp;
        for(auto& s:shards) {
            CUDA_TRY(cudaSetDevice(s->id));auto t=std::make_unique<Temporary>();t->id=s->id;
            CUDA_TRY(cudaMalloc(reinterpret_cast<void**>(&t->values),bytes));
            CUDA_TRY(cudaMalloc(reinterpret_cast<void**>(&t->out),bytes));
            CUDA_TRY(cudaMalloc(reinterpret_cast<void**>(&t->mask),r.words*8));
            CUDA_TRY(cudaMalloc(reinterpret_cast<void**>(&t->pivots),s->count*4));
            CUDA_TRY(cudaMemcpy(t->values,values.data(),bytes,cudaMemcpyHostToDevice));
            CUDA_TRY(cudaMemcpy(t->mask,free_mask.data(),r.words*8,cudaMemcpyHostToDevice));
            CUDA_TRY(cudaMemcpy(t->pivots,row_pivot.data()+s->first,s->count*4,cudaMemcpyHostToDevice));
            CUDA_TRY(cudaMemsetAsync(t->out,0,bytes));
            extract_kernel<<<unsigned(std::min<u64>(65535,s->count)),128>>>(s->matrix,s->count,r.words,t->mask,t->values,t->pivots,sw,t->out);
            CUDA_TRY(cudaGetLastError());tmp.push_back(std::move(t));
        }
        sync();auto result=values;std::vector<u64> buffer(values.size());
        for(auto& t:tmp) {
            CUDA_TRY(cudaSetDevice(t->id));CUDA_TRY(cudaMemcpy(buffer.data(),t->out,bytes,cudaMemcpyDeviceToHost));
            for(std::size_t j=0;j<result.size();++j) result[j]^=buffer[j];
        }
        return result;
    }
};
} // namespace tii
int main(int argc,char** argv) {
    return tii::cli(argc,argv,[](const tii::Request& r,const tii::Options& o) {
        return std::make_unique<tii::CudaBackend>(r,o);
    });
}
