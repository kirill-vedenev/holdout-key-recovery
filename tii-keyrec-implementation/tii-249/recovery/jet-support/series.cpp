// Public-only colex formal substitution over GF(2^m), m <= 8.
// Incremental corrections update the two newest coefficients algebraically;
// they do not approximate any convolution or discard any Hasse coefficient.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <functional>
#include <thread>
#include <vector>
using U64=uint64_t;
struct Cache {
    struct Level {std::vector<U64> masks;std::vector<uint32_t> parent;std::vector<uint8_t> last,values,dr,ds;};
    int k,d,q,stride,words,count,threads;U64 C[64][9]={};
    std::vector<uint8_t> arcs,mul;std::vector<U64> panel;std::vector<Level> levels;
    U64 rank(U64 a) {U64 r=0;int u=1;while(a){int b=__builtin_ctzll(a);r+=C[b][u++];a&=a-1;}return r;}
    template<class F> void parallel(size_t n,F f) {
        std::vector<std::thread> workers;
        for(int t=1;t<threads;++t)workers.emplace_back([&,t]{f(n*t/threads,n*(t+1)/threads);});
        f(0,n/threads);for(auto& w:workers)w.join();
    }
    Cache(int kk,int dd,int m,int modulus,int depth,int cc,int tt,const char* path):
        k(kk),d(dd),q(1<<m),stride(depth+1),words((cc+63)/64),count(cc),threads(tt),arcs(kk*(depth+1)),mul(q*q),levels(dd+1) {
        for(int n=0;n<=k;++n){C[n][0]=1;for(int r=1;r<=d;++r)C[n][r]=n?C[n-1][r-1]+C[n-1][r]:0;}
        for(int a=0;a<q;++a)for(int b=0;b<q;++b){int x=a,y=b,c=0;while(y){if(y&1)c^=x;y>>=1;x<<=1;if(x&q)x^=modulus;}mul[a*q+b]=c;}
        for(int j=0;j<=d;++j) {
            auto& l=levels[j];const U64 N=C[k][j];l.values.resize(N*stride);l.dr.resize(N);l.ds.resize(N);
            if(!j){l.masks={0};l.values[0]=1;continue;}
            U64 mask=(U64(1)<<j)-1;
            for(U64 a=0;a<N;++a){l.masks.push_back(mask);int b=63-__builtin_clzll(mask);l.last.push_back(b);l.parent.push_back(rank(mask^(U64(1)<<b)));
                U64 low=mask&(~mask+1),next=mask+low;mask=next|(((next^mask)>>2)/low);}
        }
        std::ifstream in(path,std::ios::binary|std::ios::ate);auto bytes=in.tellg();in.seekg(0);
        const U64 N=C[k][d];int source_words=U64(bytes)/(N*8);if(source_words<words || U64(bytes)!=N*source_words*8)throw 1;
        std::vector<U64> raw(N*source_words);in.read(reinterpret_cast<char*>(raw.data()),bytes);if(!in)throw 2;
        panel.resize(N*words);for(U64 a=0;a<N;++a)for(int w=0;w<words;++w)panel[a*words+w]=raw[a*source_words+w];
    }
    void set(int r,const uint8_t* jet){for(int a=0;a<k;++a)arcs[a*stride+r]=jet[a];}
    void eval(int r) {
        for(int j=1;j<=d;++j) {
            auto& l=levels[j];auto& prev=levels[j-1];
            parallel(l.masks.size(),[&](size_t start,size_t end){for(size_t a=start;a<end;++a) {
                auto* p=prev.values.data()+size_t(l.parent[a])*stride;auto* f=arcs.data()+l.last[a]*stride;uint8_t c=0;
                for(int u=0;u<=r;++u)c^=mul[p[u]*q+f[r-u]];
                l.values[a*stride+r]=c;
            }});
        }
    }
    void eval_pair(int r) {
        for(int j=1;j<=d;++j) {
            auto& l=levels[j];auto& prev=levels[j-1];
            parallel(l.masks.size(),[&](size_t start,size_t end){for(size_t a=start;a<end;++a) {
                auto* p=prev.values.data()+size_t(l.parent[a])*stride;auto* f=arcs.data()+l.last[a]*stride;uint8_t c=0,e=0;
                for(int u=0;u<=r;++u){const int t=p[u]*q;c^=mul[t+f[r-u]];e^=mul[t+f[r+1-u]];}
                e^=mul[p[r+1]*q+f[0]];
                l.values[a*stride+r]=c;l.values[a*stride+r+1]=e;
            }});
        }
    }
    // r >= 3 ensures no product of two corrections contributes at r or r+1.
    void patch_pair(int r,const uint8_t* odd,const uint8_t* even) {
        for(int j=1;j<=d;++j) {
            auto& l=levels[j];auto& prev=levels[j-1];
            parallel(l.masks.size(),[&](size_t start,size_t end){for(size_t a=start;a<end;++a) {
                size_t parent=l.parent[a];int last=l.last[a];auto* p=prev.values.data()+parent*stride;auto* f=arcs.data()+last*stride;
                uint8_t c=mul[prev.dr[parent]*q+f[0]]^mul[p[0]*q+odd[last]];
                uint8_t e=mul[prev.ds[parent]*q+f[0]]^mul[prev.dr[parent]*q+f[1]]^mul[p[1]*q+odd[last]]^mul[p[0]*q+even[last]];
                l.dr[a]=c;l.ds[a]=e;l.values[a*stride+r]^=c;l.values[a*stride+r+1]^=e;
            }});
        }
        set(r,odd);set(r+1,even);
    }
    void contract(int r,uint8_t* output) {
        std::vector<U64> buckets(q*words),bits(8*words);auto& l=levels[d];
        for(size_t a=0;a<l.masks.size();++a){uint8_t c=l.values[a*stride+r];U64* row=buckets.data()+c*words;
            for(int w=0;w<words;++w)row[w]^=panel[a*words+w];}
        for(int c=1;c<q;++c)for(int bit=0;bit<8;++bit)if((c>>bit)&1)for(int w=0;w<words;++w)bits[bit*words+w]^=buckets[c*words+w];
        for(int a=0;a<count;++a){uint8_t c=0;for(int bit=0;bit<8;++bit)c|=((bits[bit*words+a/64]>>(a%64))&1)<<bit;output[a]=c;}
    }
};
extern "C" void* cache_create(int k,int d,int m,int modulus,int depth,int count,int threads,const char* panel) {
    if(k<1||k>63||d<1||d>8||m<1||m>8||depth<4||count<1||threads<1)return nullptr;
    try{return new Cache(k,d,m,modulus,depth,count,threads,panel);}catch(...){return nullptr;}
}
extern "C" void cache_free(void* p){delete static_cast<Cache*>(p);}
extern "C" void cache_set(void* p,int r,const uint8_t* jet){static_cast<Cache*>(p)->set(r,jet);}
extern "C" void cache_eval(void* p,int r,uint8_t* output){auto& s=*static_cast<Cache*>(p);s.eval(r);s.contract(r,output);}
extern "C" void cache_pair(void* p,int r,uint8_t* output){auto& s=*static_cast<Cache*>(p);s.eval_pair(r);s.contract(r,output);s.contract(r+1,output+s.count);}
extern "C" void cache_patch(void* p,int r,const uint8_t* odd,const uint8_t* even,uint8_t* output){auto& s=*static_cast<Cache*>(p);s.patch_pair(r,odd,even);s.contract(r,output);s.contract(r+1,output+s.count);}
