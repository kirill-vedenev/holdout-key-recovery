// Exact toy-size polynomial holdout solver. Only public points enter this program.
#include <m4ri/m4ri.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <vector>
using Mask = uint64_t;
static double seconds() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}
static std::vector<Mask> subsets(int k, int d, Mask allowed=~Mask(0)) {
    std::vector<Mask> out;
    std::function<void(int,int,Mask)> visit = [&](int a,int left,Mask m) {
        if (!left) { out.push_back(m); return; }
        for (int b=a;b<=k-left;++b)
            if ((allowed>>b)&1) visit(b+1,left-1,m|(Mask(1)<<b));
    };
    visit(0,d,0); return out;
}
static bool odd_choose(int n,int r) { return (r & ~n)==0; }
static std::vector<int> orders(int d,int s,bool reduced) {
    std::vector<int> keep;
    for (int u=s-1;u>=0;--u) {
        bool implied=false;
        for (int v:keep) if (odd_choose(d-u,v-u)) implied=true;
        if (!reduced || !implied) keep.push_back(u);
    }
    std::reverse(keep.begin(),keep.end()); return keep;
}

extern "C" int local_rows(int k,int d,uint64_t point,uint8_t* output) {
    if(k<1 || k>63 || d<1 || d>k) return -1;
    auto mons=subsets(k,d); size_t stride=(mons.size()+7)/8;
    std::unordered_map<Mask,int> index; int nr=0;
    for(int u=0;u<=2;++u) for(Mask s:subsets(k,u)) index[s]=nr++;
    std::fill(output,output+nr*stride,0);
    for(size_t c=0;c<mons.size();++c) {
        Mask a=mons[c];
        for(Mask s=a;;s=(s-1)&a) {
            if(__builtin_popcountll(s)<=2 && !((a^s)&~point))
                output[index.at(s)*stride+c/8] |= 1<<(c%8);
            if(!s)break;
        }
    }
    return 0;
}

// Truncated products in lexicographic square-free monomial order. This helper
// has no kernel/secret knowledge; the caller supplies arbitrary coordinate arcs.
extern "C" int series_products(int k,int d,int m,int modulus,int depth,
                               const uint8_t* arcs,uint8_t* output) {
    if(k<1 || k>63 || d<1 || d>k || m<1 || m>8 || depth<0) return -1;
    const int q=1<<m, stride=depth+1;
    std::vector<uint8_t> mul(q*q);
    for(int a=0;a<q;++a) for(int b=0;b<q;++b) {
        int x=a,y=b,c=0;
        while(y) { if(y&1)c^=x; y>>=1; x<<=1; if(x&q)x^=modulus; }
        mul[a*q+b]=c;
    }
    std::vector<std::vector<uint8_t>> stack(d+1,std::vector<uint8_t>(stride));
    stack[0][0]=1; size_t col=0;
    std::function<void(int,int)> visit=[&](int start,int level) {
        if(level==d) {
            std::copy(stack[level].begin(),stack[level].end(),output+(col++)*stride);
            return;
        }
        for(int a=start;a<=k-(d-level);++a) {
            auto& dest=stack[level+1]; std::fill(dest.begin(),dest.end(),0);
            const uint8_t* arc=arcs+a*stride;
            for(int u=0;u<=depth;++u) if(stack[level][u])
                for(int v=0;v<=depth-u;++v) if(arc[v])
                    dest[u+v]^=mul[stack[level][u]*q+arc[v]];
            visit(a+1,level+1);
        }
    };
    visit(0,0); return 0;
}

// Incremental formal products: cache all prefixes, and recompute only the
// coefficients changed by a newly chosen jet. This avoids cubic work in depth.
struct SeriesCache {
    struct Level { std::vector<Mask> masks; std::vector<int> parent,last; std::vector<uint8_t> values; };
    int k,d,q,stride; std::vector<uint8_t> arcs,mul; std::vector<Level> levels;
    SeriesCache(int kk,int dd,int m,int modulus,int depth):k(kk),d(dd),q(1<<m),stride(depth+1),arcs(k*stride),mul(q*q),levels(d+1) {
        for(int a=0;a<q;++a) for(int b=0;b<q;++b) {
            int x=a,y=b,c=0;
            while(y) { if(y&1)c^=x; y>>=1; x<<=1; if(x&q)x^=modulus; }
            mul[a*q+b]=c;
        }
        levels[0].masks={0}; levels[0].values.resize(stride); levels[0].values[0]=1;
        for(int j=1;j<=d;++j) {
            auto& level=levels[j]; level.masks=subsets(k,j);
            std::unordered_map<Mask,int> index;
            for(size_t a=0;a<levels[j-1].masks.size();++a) index[levels[j-1].masks[a]]=a;
            level.values.resize(level.masks.size()*stride);
            for(Mask mon:level.masks) {
                int last=63-__builtin_clzll(mon);
                level.last.push_back(last); level.parent.push_back(index.at(mon^(Mask(1)<<last)));
            }
        }
    }
    void evaluate(int r,uint8_t* output) {
        for(int j=1;j<=d;++j) {
            auto& level=levels[j]; const auto& previous=levels[j-1];
            for(size_t a=0;a<level.masks.size();++a) {
                const uint8_t* p=previous.values.data()+size_t(level.parent[a])*stride;
                const uint8_t* f=arcs.data()+level.last[a]*stride; uint8_t c=0;
                for(int u=0;u<=r;++u) if(p[u] && f[r-u]) c^=mul[p[u]*q+f[r-u]];
                level.values[a*stride+r]=c;
            }
        }
        for(size_t a=0;a<levels[d].masks.size();++a) output[a]=levels[d].values[a*stride+r];
    }
};
extern "C" void* series_cache_create(int k,int d,int m,int modulus,int depth) {
    if(k<1 || k>63 || d<1 || d>k || m<1 || m>8 || depth<0) return nullptr;
    try { return new SeriesCache(k,d,m,modulus,depth); } catch(...) { return nullptr; }
}
extern "C" void series_cache_free(void* ptr) { delete static_cast<SeriesCache*>(ptr); }
extern "C" int series_cache_set(void* ptr,int r,const uint8_t* values) {
    auto& s=*static_cast<SeriesCache*>(ptr); if(r<0 || r>=s.stride)return -1;
    for(int a=0;a<s.k;++a) { if(values[a]>=s.q)return -1; s.arcs[a*s.stride+r]=values[a]; }
    return 0;
}
extern "C" int series_cache_eval(void* ptr,int r,uint8_t* output) {
    auto& s=*static_cast<SeriesCache*>(ptr); if(r<0 || r>=s.stride)return -1;
    s.evaluate(r,output); return 0;
}

#ifndef SHARED_ONLY
int main(int argc,char** argv) {
  try {
    if(argc!=4) throw std::runtime_error("usage: holdout public.txt kernel.bin max_GiB");
    const double start=seconds();
    std::ifstream in(argv[1]); int k,d,points;
    if(!(in>>k>>d>>points) || k<1 || k>63 || d<2 || d>k)
        throw std::runtime_error("invalid public header (require 2<=d<=k<=63)");
    std::vector<Mask> p(points); std::vector<int> mult(points);
    for(int j=0;j<points;++j)
        if(!(in>>p[j]>>mult[j]) || (p[j]>>k) || mult[j]<0 || mult[j]>=d)
            throw std::runtime_error("invalid public point/multiplicity");
    auto mons=subsets(k,d); const int N=mons.size();
    std::unordered_map<Mask,int> mon_index;
    for(int c=0;c<N;++c) mon_index[mons[c]]=c;
    std::vector<std::vector<Mask>> indices(d);
    for(int u=0;u<d;++u) indices[u]=subsets(k,u);
    std::vector<std::unordered_map<Mask,int>> slots(points);
    int rows=0; uint64_t full_rows=0;
    for(int j=0;j<points;++j) {
        for(int u:orders(d,mult[j],true)) for(Mask s:indices[u]) slots[j][s]=rows++;
        for(int u=0;u<mult[j];++u) full_rows+=indices[u].size();
    }
    double bytes=double(rows)*((N+63)/64)*8;
    if(bytes>std::stod(argv[3])*1073741824.) throw std::runtime_error("matrix exceeds memory cap");
    std::cerr<<"ASSEMBLE rows="<<rows<<" columns="<<N<<" packed_GiB="<<bytes/1073741824.<<std::endl;
    mzd_t* A=mzd_init(rows,N); uint64_t nnz=0;
    for(int j=0;j<points;++j) for(int c=0;c<N;++c) {
        Mask a=mons[c];
        for(Mask s=a;;s=(s-1)&a) {
            if(!((a^s)&~p[j])) {
                auto it=slots[j].find(s);
                if(it!=slots[j].end()) { mzd_write_bit(A,it->second,c,1); ++nnz; }
            }
            if(!s) break;
        }
    }
    const double built=seconds(); std::cerr<<"ELIMINATE nnz="<<nnz<<std::endl;
    int rank=mzd_echelonize_m4ri(A,1,0);
    const double reduced=seconds(); int nullity=N-rank;
    std::vector<int> pivot(rank),free; std::vector<bool> is_pivot(N,false);
    int cursor=0;
    for(int r=0;r<rank;++r) {
        while(cursor<N && !mzd_read_bit(A,r,cursor)) ++cursor;
        if(cursor==N) throw std::runtime_error("invalid echelon output");
        pivot[r]=cursor; is_pivot[cursor++]=true;
    }
    for(int c=0;c<N;++c) if(!is_pivot[c]) free.push_back(c);
    mzd_t* KT=mzd_init(N,nullity);
    for(int f=0;f<nullity;++f) {
        mzd_write_bit(KT,free[f],f,1);
        for(int r=0;r<rank;++r) if(mzd_read_bit(A,r,free[f])) mzd_write_bit(KT,pivot[r],f,1);
    }
    mzd_free(A);
    const double extracted=seconds();
    std::cerr<<"VERIFY full Hasse constraints; kernel_dimension="<<nullity<<std::endl;
    // Independently assemble by derivative index and extensions, not monomial
    // submasks. Verify every returned vector against every original order.
    for(int j=0;j<points;++j) {
        int nr=0; for(int u=0;u<mult[j];++u) nr+=indices[u].size();
        if(!nr || !nullity) continue;
        mzd_t* B=mzd_init(nr,N); int row=0;
        for(int u=0;u<mult[j];++u) for(Mask s:indices[u]) {
            for(Mask extra:subsets(k,d-u,p[j]&~s))
                mzd_write_bit(B,row,mon_index.at(s|extra),1);
            ++row;
        }
        mzd_t* check=mzd_mul(nullptr,B,KT,0);
        if(!mzd_is_zero(check)) throw std::runtime_error("full Hasse verification failed");
        mzd_free(check); mzd_free(B);
    }
    const double verified=seconds();
    std::ofstream out(argv[2],std::ios::binary);
    out.write("WHK1",4);
    uint32_t header[4]={uint32_t(k),uint32_t(d),uint32_t(N),uint32_t(nullity)};
    out.write(reinterpret_cast<char*>(header),sizeof(header));
    std::vector<uint8_t> packed((N+7)/8);
    for(int f=0;f<nullity;++f) {
        std::fill(packed.begin(),packed.end(),0);
        for(int c=0;c<N;++c) if(mzd_read_bit(KT,c,f)) packed[c/8]|=1<<(c%8);
        out.write(reinterpret_cast<char*>(packed.data()),packed.size());
    }
    if(!out) throw std::runtime_error("kernel output failed");
    mzd_free(KT);
    std::cout<<"{\"rows\":"<<rows<<",\"full_rows\":"<<full_rows
        <<",\"columns\":"<<N<<",\"rank\":"<<rank<<",\"kernel_dimension\":"<<nullity
        <<",\"nnz\":"<<nnz<<",\"packed_matrix_bytes\":"<<uint64_t(bytes)
        <<",\"full_hasse_verified\":true,\"build_seconds\":"<<built-start
        <<",\"elimination_seconds\":"<<reduced-built
        <<",\"extraction_seconds\":"<<extracted-reduced
        <<",\"verification_seconds\":"<<verified-extracted
        <<",\"total_seconds\":"<<seconds()-start<<"}"<<std::endl;
  } catch(const std::exception& e) { std::cerr<<e.what()<<std::endl; return 1; }
}
#endif
