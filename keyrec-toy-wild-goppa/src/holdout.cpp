#include "public.hpp"
#include <m4ri/m4ri.h>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <unordered_map>

using namespace toy;using namespace wild;
using Clock=std::chrono::steady_clock;
static double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
struct Free {void operator()(mzd_t* a) const {if(a)mzd_free(a);}};
using Packed=std::unique_ptr<mzd_t,Free>;
static void write_block(mzd_t* a,int r,int c,Byte value) {
    // Multiplication by value in the basis (1,w), w^2=w+1.
    Byte first=value,second=mul4[value][2];
    if(first&1)mzd_write_bit(a,2*r,2*c,1);
    if(first&2)mzd_write_bit(a,2*r+1,2*c,1);
    if(second&1)mzd_write_bit(a,2*r,2*c+1,1);
    if(second&2)mzd_write_bit(a,2*r+1,2*c+1,1);
}
int main(int argc,char** argv) {
    try {
        require(argc>=3 && argc<=5,"usage: holdout PUBLIC KERNEL [PANEL_SIZE=256] [MATRIX_GIB=16]");
        auto start=Clock::now();Public p=read_public(argv[1]);int cap=argc>=4?std::stoi(argv[3]):256;
        double max_gib=argc>=5?std::stod(argv[4]):16;require(cap>=1 && cap<=4096,"panel size outside 1..4096");
        auto mons=masks(p.k,p.d);auto points=p.retained_nonunits();auto orders=derivative_orders(p.d,p.s);
        std::vector<std::uint64_t> derivatives;for(int r:orders){auto v=masks(p.k,r);derivatives.insert(derivatives.end(),v.begin(),v.end());}
        std::unordered_map<std::uint64_t,int> slots;for(int r=0;r<int(derivatives.size());++r)slots[derivatives[r]]=r;
        int N=int(mons.size()),rows=int(points.size()*derivatives.size());
        std::uint64_t bytes=std::uint64_t(2*rows)*((2*N+63)/64)*8;
        require(bytes<=max_gib*1073741824.,"packed matrix exceeds allocation cap");
#if !__M4RI_HAVE_OPENMP
        require(bytes<1073741824ull || std::getenv("WILD_ALLOW_SERIAL")!=nullptr,
                "large kernel requires OpenMP-enabled M4RI; build with M4RI_OMP or explicitly set WILD_ALLOW_SERIAL");
#endif
        const char* thread_limit=std::getenv("OMP_NUM_THREADS");
        int requested_threads=thread_limit?std::stoi(thread_limit):0;
        std::cerr<<"ASSEMBLE "<<2*rows<<" x "<<2*N<<" binary matrix; "<<bytes/1073741824.
            <<" GiB; OpenMP="<<(__M4RI_HAVE_OPENMP?"enabled":"disabled")<<"; requested threads="<<requested_threads<<'\n';
        Packed A(mzd_init(2*rows,2*N));
        for(int block=0;block<int(points.size());++block)for(int c=0;c<N;++c) {
            auto mon=mons[c];for(auto derivative=mon;;derivative=(derivative-1)&mon) {
                auto found=slots.find(derivative);if(found!=slots.end()) {
                    Byte value=monomial_value(p,points[block],mon^derivative);
                    if(value)write_block(A.get(),block*int(derivatives.size())+found->second,c,value);
                }
                if(!derivative)break;
            }
        }
        double assembly=elapsed(start);auto tick=Clock::now();std::cerr<<"ELIMINATE\n";
        int rank2=mzd_echelonize_m4ri(A.get(),1,0);double elimination=elapsed(tick);tick=Clock::now();
        require(rank2%2==0,"expanded rank is not even");
        std::vector<int> pivots(rank2),free;std::vector<bool> used(2*N);int cursor=0;
        for(int r=0;r<rank2;++r){while(cursor<2*N && !mzd_read_bit(A.get(),r,cursor))++cursor;require(cursor<2*N,"invalid echelon output");pivots[r]=cursor;used[cursor++]=true;}
        for(int c=0;c<N;++c){require(used[2*c]==used[2*c+1],"expanded pivots do not respect GF(4) column pairs");if(!used[2*c])free.push_back(c);}
        int full_dimension=int(free.size());require(full_dimension>0,"holdout kernel is zero");
        // Distinct free GF(4) columns give independent kernel vectors. Spread the
        // selected columns across the complete basis with a fixed public shuffle.
        std::mt19937_64 rng(0);for(int i=int(free.size())-1;i>0;--i)std::swap(free[i],free[rng()%unsigned(i+1)]);
        free.resize(std::min(cap,full_dimension));int M=int(free.size());
        Matrix K(M,N);
        for(int r=0;r<M;++r) {
            int free_bit=2*free[r];K(r,free[r])=1;
            for(int a=0;a<rank2;++a)if(mzd_read_bit(A.get(),a,free_bit))K(r,pivots[a]/2)|=Byte(1u<<(pivots[a]%2));
        }
        A.reset();
        std::ofstream out(argv[2],std::ios::binary);out.write("WGK1",4);
        for(int x:{p.k,p.n,p.D(),p.held,p.d,N,M,full_dimension})write_u32(out,std::uint32_t(x));
        Vec packed((N+3)/4);
        for(int r=0;r<M;++r){std::fill(packed.begin(),packed.end(),0);for(int c=0;c<N;++c)packed[c/4]|=K(r,c)<<(2*(c%4));out.write(reinterpret_cast<char*>(packed.data()),packed.size());}
        require(bool(out),"kernel write failed");
        std::cout<<"{\"status\":\"computed\",\"heldout_positions\":1,\"gf4_rows\":"<<rows<<",\"gf4_columns\":"<<N
            <<",\"gf2_rows\":"<<2*rows<<",\"gf2_columns\":"<<2*N<<",\"rank_gf4\":"<<rank2/2
            <<",\"kernel_dimension_gf4\":"<<full_dimension<<",\"exported_polynomials\":"<<M<<",\"packed_matrix_bytes\":"<<bytes
            <<",\"openmp_enabled\":"<<(__M4RI_HAVE_OPENMP?"true":"false")<<",\"requested_threads\":"<<requested_threads
            <<",\"assembly_seconds\":"<<assembly<<",\"elimination_seconds\":"<<elimination<<",\"extraction_seconds\":"<<elapsed(tick)
            <<",\"total_seconds\":"<<elapsed(start)<<"}\n";
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
