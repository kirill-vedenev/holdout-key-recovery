// Exact Hasse Taylor coefficients at a binary point, in the saved colex basis.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using U64=uint64_t;
namespace fs=std::filesystem;
int main(int argc,char** argv) { try {
    if(argc!=8) throw std::runtime_error("usage: extract_jets panel monomials k degree point count output-dir");
    const int k=std::stoi(argv[3]),d=std::stoi(argv[4]),count=std::stoi(argv[6]);
    const U64 p=std::stoull(argv[5]); const int words=(count+63)/64;
    if(k>63 || d>8 || p>>k) throw std::runtime_error("unsupported dimensions");
    U64 C[64][9]={};for(int n=0;n<=k;++n) {C[n][0]=1;for(int r=1;r<=d;++r)C[n][r]=n?C[n-1][r-1]+C[n-1][r]:0;}
    auto rank=[&](U64 a) {U64 r=0;int u=1;while(a){int b=__builtin_ctzll(a);r+=C[b][u++];a&=a-1;}return r;};
    const U64 N=C[k][d];
    const U64 panel_bytes=fs::file_size(argv[1]),map_bytes=fs::file_size(argv[2]);
    if(map_bytes!=N*8 || panel_bytes%(N*8)) throw std::runtime_error("invalid input size");
    const int stride=panel_bytes/(N*8); if(stride<words)throw std::runtime_error("insufficient panel words");
    std::vector<U64> panel(N*stride),mons(N);
    std::ifstream f(argv[1],std::ios::binary),g(argv[2],std::ios::binary);
    f.read(reinterpret_cast<char*>(panel.data()),panel_bytes);g.read(reinterpret_cast<char*>(mons.data()),map_bytes);
    if(!f||!g)throw std::runtime_error("input read failed");
    std::vector<std::vector<U64>> jets(d);
    for(int u=0;u<d;++u) jets[u].resize(C[k][u]*words);
    U64 updates=0;
    for(U64 c=0;c<N;++c) {
        const U64 a=mons[c];if(__builtin_popcountll(a)!=d || a>>k || rank(a)!=c)throw std::runtime_error("input map is not the declared colex basis");
        for(int w=words;w<stride;++w)if(panel[c*stride+w])throw std::runtime_error("nonzero unused panel column");
        if(count%64 && panel[c*stride+words-1]>>(count%64))throw std::runtime_error("nonzero padding coefficient");
        const U64 mandatory=a&~p,optional=a&p;
        for(U64 b=optional;;b=(b-1)&optional) {
            U64 derivative=mandatory|b;int u=__builtin_popcountll(derivative);
            if(u<d) {
                U64* dst=jets[u].data()+rank(derivative)*words;
                for(int w=0;w<words;++w)dst[w]^=panel[c*stride+w];++updates;
            }
            if(!b)break;
        }
    }
    fs::create_directories(argv[7]);
    for(int u=0;u<d;++u) {
        fs::path path=fs::path(argv[7])/("hasse-"+std::to_string(u)+".u64le");
        if(fs::exists(path))throw std::runtime_error("output already exists");
        std::ofstream out(path,std::ios::binary);out.write(reinterpret_cast<char*>(jets[u].data()),jets[u].size()*8);
        if(!out)throw std::runtime_error("output write failed");
    }
    std::cout<<"{\"k\":"<<k<<",\"degree\":"<<d<<",\"count\":"<<count<<",\"words\":"<<words<<",\"point\":"<<p<<",\"updates\":"<<updates<<",\"colex_map_checked\":true,\"padding_checked\":true}"<<std::endl;
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;} }
