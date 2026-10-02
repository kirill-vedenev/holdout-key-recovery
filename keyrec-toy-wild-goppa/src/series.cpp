// Exact truncated polynomial products and incremental coefficient cache.
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

