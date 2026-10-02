#pragma once
#include "../third_party/saarinen/colex.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <csignal>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <sstream>
#include <thread>
#include <vector>

namespace tii {
using u64 = std::uint64_t;
using u32 = std::uint32_t;
namespace fs = std::filesystem;
constexpr u32 NONE = std::numeric_limits<u32>::max();
inline volatile std::sig_atomic_t interrupted = 0;
inline void on_signal(int) { interrupted = 1; }
inline double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
inline u64 choose(int n, int d) { return saarinen::HOST_CHOOSE.get(n,d); }
inline u64 rank_mask(u64 v) { return saarinen::colex_rank_host(v); }
inline u64 unrank(u64 r,int d,int k) { return saarinen::colex_unrank_host(r,d,k); }
inline void require(bool ok,const std::string& why) { if (!ok) throw std::runtime_error(why); }
inline std::string escape(const std::string& s) {
    std::ostringstream o;
    for (unsigned char c:s) {
        if (c=='"' || c=='\\') o << '\\' << c;
        else if (c<32) o << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
        else o << c;
    }
    return o.str();
}
inline void put32(std::ostream& f,u32 v) { for(int b=0;b<4;++b) f.put(char(v>>(8*b))); }
inline void put64(std::ostream& f,u64 v) { for(int b=0;b<8;++b) f.put(char(v>>(8*b))); }
inline u64 get_number(std::istream& f,int bytes) {
    u64 v=0;
    for(int b=0;b<bytes;++b) { int c=f.get(); require(c!=EOF,"truncated binary file"); v|=u64(c)<<(8*b); }
    return v;
}
inline void finish_file(std::ofstream& f,const fs::path& temporary,const fs::path& final) {
    f.flush(); require(bool(f),"output write failed: "+temporary.string()); f.close();
    fs::rename(temporary,final);
}
inline void json_file(const fs::path& path,const std::string& contents) {
    fs::path temp=path.string()+".tmp";
    std::ofstream f(temp); require(bool(f),"cannot create "+temp.string());
    f << contents << '\n'; finish_file(f,temp,path);
}
struct Point { int original; u64 mask; int multiplicity; };
struct Row { u64 point, derivative; };
struct Request {
    std::string id;
    int k=0,d=0,minimum=0,heldout_original=0;
    u64 heldout_mask=0,N=0,words=0;
    std::vector<Point> points;
    std::vector<Row> rows;
};
inline std::vector<int> orders(int d,int s) {
    std::vector<int> keep;
    for(int u=s-1;u>=0;--u) {
        bool implied=false;
        for(int v:keep) if (((v-u)&~(d-u))==0) implied=true;
        if(!implied) keep.push_back(u);
    }
    std::sort(keep.begin(),keep.end()); return keep;
}
inline Request read_request(const fs::path& path) {
    std::ifstream in(path); require(bool(in),"cannot open request");
    std::string magic; in>>magic; require(magic=="TII_HOLDOUT_V1","bad request magic");
    auto key=[&](const char* expected) { std::string s; in>>s; require(s==expected,std::string("expected ")+expected); };
    Request r; key("instance_id"); in>>r.id;
    require(r.id.size()==64 && r.id.find_first_not_of("0123456789abcdef")==std::string::npos,"bad instance identity");
    key("k"); in>>r.k; key("degree"); in>>r.d;
    require(r.k>=2 && r.k<=63 && r.d>=2 && r.d<=8 && r.d<=r.k,"unsupported k/degree");
    key("minimum"); in>>r.minimum; require(r.minimum>0,"minimum must be positive");
    key("heldout_original"); in>>r.heldout_original;
    key("heldout_mask"); in>>r.heldout_mask;
    require(r.heldout_original>=0 && !(r.heldout_mask>>r.k),"bad held-out position");
    int count; key("points"); in>>count; require(count>=0 && count<=10000,"bad point count");
    std::vector<int> seen;
    for(int i=0;i<count;++i) {
        Point p; in>>p.original>>p.mask>>p.multiplicity;
        require(bool(in) && p.original>=-1 && p.original!=r.heldout_original && !(p.mask>>r.k)
                && p.multiplicity>0 && p.multiplicity<r.d,"invalid constraint point");
        require(std::find(seen.begin(),seen.end(),p.original)==seen.end(),"duplicate original point label");
        seen.push_back(p.original); r.points.push_back(p);
    }
    std::string extra; require(!(in>>extra),"trailing request tokens");
    r.N=choose(r.k,r.d); require(r.N>0 && r.N<NONE,"monomial count exceeds supported u32 indexing");
    r.words=(r.N+63)/64;
    for(const Point& p:r.points) {
        // If fewer than d-s+1 coordinates are nonzero, every relevant derivative is identically zero.
        if(__builtin_popcountll(p.mask)<r.d-p.multiplicity+1) continue;
        for(int u:orders(r.d,p.multiplicity))
            for(u64 a=0;a<choose(r.k,u);++a) r.rows.push_back({p.mask,unrank(a,u,r.k)});
    }
    require(r.rows.size()<NONE,"constraint count exceeds supported u32 indexing");
    return r;
}
template<class F> inline void combinations(u64 allowed,int count,F visit) {
    if(count<0 || __builtin_popcountll(allowed)<count) return;
    if(count==0) { visit(u64(0)); return; }
    int pos[63],n=0; while(allowed) { pos[n++]=__builtin_ctzll(allowed); allowed&=allowed-1; }
    u64 v=(u64(1)<<count)-1,limit=u64(1)<<n;
    while(v<limit) {
        u64 a=v,mask=0;
        while(a) { mask|=u64(1)<<pos[__builtin_ctzll(a)]; a&=a-1; }
        visit(mask);
        u64 low=v&(~v+1),next=v+low;
        if(next>=limit) break;
        v=next|(((next^v)>>2)/low);
    }
}
inline std::vector<u64> assemble_cpu(const Request& r) {
    std::vector<u64> a(r.rows.size()*r.words);
    for(std::size_t i=0;i<r.rows.size();++i) {
        auto row=r.rows[i]; int extra=r.d-__builtin_popcountll(row.derivative);
        combinations(row.point&~row.derivative,extra,[&](u64 x) {
            u64 c=rank_mask(x|row.derivative); a[i*r.words+c/64]|=u64(1)<<(c%64);
        });
    }
    return a;
}
struct Kernel {
    u32 k,d,N,count;
    std::vector<unsigned char> data;
    std::size_t stride() const { return (u64(N)+7)/8; }
    bool bit(u32 row,u32 col) const { return (data[row*stride()+col/8]>>(col%8))&1; }
};
inline Kernel read_kernel(const fs::path& path) {
    std::ifstream in(path,std::ios::binary); require(bool(in),"cannot open kernel file");
    char magic[4]; in.read(magic,4); require(in && std::memcmp(magic,"THK1",4)==0,"expected THK1 colex kernel");
    Kernel k; k.k=get_number(in,4); k.d=get_number(in,4); k.N=get_number(in,4); k.count=get_number(in,4);
    require(k.k<=63 && k.d>=2 && k.d<=8 && k.d<=k.k && k.N==choose(k.k,k.d)
            && k.count>0 && k.count<=k.N,"invalid kernel dimensions");
    u64 bytes=k.stride()*u64(k.count);
    require(fs::file_size(path)==20+bytes,"kernel file length mismatch");
    require(bytes<=u64(8)<<30,"kernel reader memory cap exceeded");
    k.data.resize(bytes); in.read(reinterpret_cast<char*>(k.data.data()),bytes); require(bool(in),"kernel read failed");
    if(k.N%8) for(u32 i=0;i<k.count;++i)
        require((k.data[(i+1)*k.stride()-1]>>(k.N%8))==0,"nonzero kernel padding bits");
    return k;
}
inline void write_kernel(const fs::path& path,const Request& r,u32 count,const std::vector<u64>& columns) {
    require(!fs::exists(path),"kernel output already exists");
    const u64 sw=(u64(count)+63)/64; fs::path tmp=path.string()+".tmp";
    std::ofstream f(tmp,std::ios::binary); f.write("THK1",4);
    put32(f,r.k); put32(f,r.d); put32(f,r.N); put32(f,count);
    std::vector<unsigned char> row((r.N+7)/8);
    for(u32 v=0;v<count;++v) {
        std::fill(row.begin(),row.end(),0);
        for(u64 c=0;c<r.N;++c) row[c/8]|=((columns[c*sw+v/64]>>(v%64))&1)<<(c%8);
        f.write(reinterpret_cast<char*>(row.data()),row.size());
    }
    finish_file(f,tmp,path);
}
inline u32 independent_rank(const Kernel& k) {
    const std::size_t words=(u64(k.N)+63)/64;
    std::vector<u64> basis(std::size_t(k.count)*words),v(words);
    std::vector<int> owner(k.N,-1); u32 rank=0;
    for(u32 r=0;r<k.count;++r) {
        std::fill(v.begin(),v.end(),0);
        for(std::size_t j=0;j<k.stride();++j) v[j/8]|=u64(k.data[r*k.stride()+j])<<(8*(j%8));
        for(std::size_t w=0;w<words;++w) while(v[w]) {
            u32 p=w*64+__builtin_ctzll(v[w]);
            if(owner[p]<0) {
                owner[p]=rank; std::copy(v.begin(),v.end(),basis.begin()+std::size_t(rank)*words); ++rank;
                goto next_vector;
            }
            const u64* b=basis.data()+std::size_t(owner[p])*words;
            for(std::size_t j=w;j<words;++j) v[j]^=b[j];
        }
        next_vector:;
    }
    return rank;
}
inline std::vector<u64> kernel_columns(const Kernel& k) {
    u64 sw=(u64(k.count)+63)/64; std::vector<u64> columns(u64(k.N)*sw);
    for(u32 v=0;v<k.count;++v) for(u32 c=0;c<k.N;++c)
        if(k.bit(v,c)) columns[u64(c)*sw+v/64]|=u64(1)<<(v%64);
    return columns;
}
inline std::string verify(const Request& r,const Kernel& k,int threads) {
    require(k.k==u32(r.k) && k.d==u32(r.d) && k.N==r.N,"kernel/request dimension mismatch");
    require(k.count>=u32(r.minimum),"too few kernel elements");
    double start=now(); u32 rank=independent_rank(k);
    require(rank==k.count,"kernel elements are not independent");
    auto columns=kernel_columns(k); const u64 sw=(u64(k.count)+63)/64;
    std::vector<u64> mons(r.N); for(u64 c=0;c<r.N;++c) mons[c]=unrank(c,r.d,r.k);
    std::atomic<std::size_t> next{0}; std::atomic<bool> failed{false};
    std::atomic<u64> checked{0},implicit{0};
    std::exception_ptr error; std::mutex error_mutex;
    auto worker=[&] {
      try {
        for(;;) {
            std::size_t j=next.fetch_add(1); if(j>=r.points.size()) break;
            auto p=r.points[j]; std::vector<u64> offsets(p.multiplicity); u64 rows=0;
            for(int u=0;u<p.multiplicity;++u) { offsets[u]=rows; rows+=choose(r.k,u); }
            if(__builtin_popcountll(p.mask)<r.d-p.multiplicity+1) { implicit+=rows; continue; }
            std::vector<u64> residual(rows*sw);
            // Independent assembly orientation: visit every monomial and all its derivative submasks.
            for(u64 c=0;c<r.N;++c) {
                u64 a=mons[c];
                for(u64 s=a;;s=(s-1)&a) {
                    int u=__builtin_popcountll(s);
                    if(u<p.multiplicity && !((a^s)&~p.mask)) {
                        u64 dest=(offsets[u]+rank_mask(s))*sw;
                        for(u64 w=0;w<sw;++w) residual[dest+w]^=columns[c*sw+w];
                    }
                    if(!s) break;
                }
            }
            if(std::any_of(residual.begin(),residual.end(),[](u64 x){return x!=0;})) failed=true;
            checked+=rows;
        }
      } catch(...) { std::lock_guard<std::mutex> lock(error_mutex); if(!error) error=std::current_exception(); }
    };
    std::vector<std::thread> pool;
    for(int i=0;i<std::max(1,std::min<int>(threads,r.points.size()));++i) pool.emplace_back(worker);
    for(auto& t:pool) t.join(); if(error) std::rethrow_exception(error);
    require(!failed,"original Hasse constraints failed");
    // The held-out value is diagnostic: it is not silently imposed as an extra equation.
    std::vector<u64> value(sw);
    for(u64 c=0;c<r.N;++c) if(!(mons[c]&~r.heldout_mask))
        for(u64 w=0;w<sw;++w) value[w]^=columns[c*sw+w];
    u32 held_nonzero=0; for(u64 x:value) held_nonzero+=__builtin_popcountll(x);
    std::ostringstream out; out << std::setprecision(12)
        << "{\"schema\":\"tii-holdout-verification-v1\",\"instance_id\":\""<<r.id
        <<"\",\"single_position\":true,\"heldout_original\":"<<r.heldout_original
        <<",\"verified_elements\":"<<k.count<<",\"independent_rank\":"<<rank
        <<",\"minimum_elements\":"<<r.minimum<<",\"all_original_orders_verified\":true"
        <<",\"explicit_constraint_rows_checked\":"<<checked.load()
        <<",\"identically_zero_constraint_rows_checked\":"<<implicit.load()
        <<",\"heldout_nonzero_values\":"<<held_nonzero
        <<",\"curve_identity_certified\":false,\"seconds\":"<<now()-start<<"}";
    return out.str();
}
inline void export_whk(const Kernel& k,const fs::path& output) {
    require(!fs::exists(output),"export output already exists");
    std::vector<u32> order; order.reserve(k.N);
    std::function<void(int,int,u64)> visit=[&](int start,int left,u64 mask) {
        if(!left) { order.push_back(rank_mask(mask)); return; }
        for(int a=start;a<=int(k.k)-left;++a) visit(a+1,left-1,mask|(u64(1)<<a));
    }; visit(0,k.d,0);
    fs::path tmp=output.string()+".tmp"; std::ofstream f(tmp,std::ios::binary);
    f.write("WHK1",4); put32(f,k.k);put32(f,k.d);put32(f,k.N);put32(f,k.count);
    std::vector<unsigned char> row(k.stride());
    for(u32 v=0;v<k.count;++v) {
        std::fill(row.begin(),row.end(),0);
        for(u32 c=0;c<k.N;++c) if(k.bit(v,order[c])) row[c/8]|=1<<(c%8);
        f.write(reinterpret_cast<char*>(row.data()),row.size());
    }
    finish_file(f,tmp,output);
}
} // namespace tii
