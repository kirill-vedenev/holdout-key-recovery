#pragma once
#include "algebra.hpp"
#include <fstream>
#include <string>

namespace wild {
using namespace toy;
static constexpr Byte embedding[4]={0,1,58,59};
static constexpr Byte mul4[4][4]={{0,0,0,0},{0,1,2,3},{0,2,3,1},{0,3,1,2}};
struct Public {
    int m=3,t=0,n=0,k=0,modulus=67,held=0,d=0,s=0;
    std::vector<std::uint64_t> points;
    int D() const { return n-4*t-1; }
    Byte digit(int a,int j) const { return (points[j]>>(2*a))&3; }
    Vec point(int j) const { Vec v(k);for(int a=0;a<k;++a)v[a]=embedding[digit(a,j)];return v; }
    Matrix generator() const { Matrix y(k,n);for(int j=0;j<n;++j)for(int a=0;a<k;++a)y(a,j)=embedding[digit(a,j)];return y; }
    int unit_row(int j) const {
        int row=-1;for(int a=0;a<k;++a)if(Byte x=digit(a,j)) {
            if(x!=1 || row>=0)return -1;row=a;
        }
        return row;
    }
    std::vector<int> retained_nonunits() const {
        std::vector<int> out;std::vector<bool> units(k);
        for(int j=0;j<n;++j)if(j!=held){int a=unit_row(j);if(a>=0)units[a]=true;else out.push_back(j);}
        require(std::all_of(units.begin(),units.end(),[](bool x){return x;}),"retained columns must include an information set of unit vectors");return out;
    }
};
inline Public read_public(const std::string& path) {
    std::ifstream in(path);std::string magic;Public p;
    require(bool(in>>magic>>p.m>>p.t>>p.n>>p.k>>p.modulus>>p.held>>p.d>>p.s) && magic=="WGPK1","invalid public header");
    require(p.m==3 && p.modulus==67 && p.t>=2 && p.n<=64 && p.n>4*p.t && p.k>=5 && p.k<=31 &&
            p.d>=3 && p.d<=6 && p.s==p.d-1 && p.held>=0 && p.held<p.n,"unsupported public parameters");
    p.points.resize(p.n);for(auto& x:p.points)require(bool(in>>x) && !(x>>(2*p.k)),"invalid public column");
    std::string extra;require(!(in>>extra),"trailing public data");Field f(6,67);
    require(rank(f,p.generator())==p.k,"public generator rank");p.retained_nonunits();
    require(p.s*(p.n-1)>p.d*p.D(),"ordinary interpolation bound is not satisfied");return p;
}
inline Byte monomial_value(const Public& p,int point,std::uint64_t mask) {
    Byte value=1;while(mask){int a=__builtin_ctzll(mask);value=mul4[value][p.digit(a,point)];mask&=mask-1;}return value;
}
inline std::vector<int> derivative_orders(int d,int s) {
    std::vector<int> orders;
    for(int u=s-1;u>=0;--u){bool implied=false;for(int v:orders)if(((v-u)&~(d-u))==0)implied=true;if(!implied)orders.push_back(u);}
    std::reverse(orders.begin(),orders.end());return orders;
}
inline void write_u32(std::ostream& out,std::uint32_t x) {
    for(int b=0;b<4;++b)out.put(char((x>>(8*b))&255));
}
inline std::uint32_t read_u32(std::istream& in) {
    std::uint32_t x=0;for(int b=0;b<4;++b){int a=in.get();require(a!=EOF,"truncated kernel header");x|=std::uint32_t(a)<<(8*b);}return x;
}
} // namespace wild
