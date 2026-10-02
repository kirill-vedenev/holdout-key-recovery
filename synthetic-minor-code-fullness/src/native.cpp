// C++/M4RIE backend.
//
// minor_*: minor code of a public matrix and supplied tangent representatives
// (used for the GRS, alternant and wild Goppa families, generated in Sage).
// goppa_trial: one complete binary Goppa trial (key, shortening, tangent planes,
// minor code, certificate and recovery).
//
// Recovery sees only the public matrix, the minor matrix and D. The secret
// support enters only the generation of the tangent planes and the checks.
#include <m4rie/m4rie.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using u16 = std::uint16_t;
using Clock = std::chrono::steady_clock;
struct MatFree { void operator()(mzed_t* p) const { if (p) mzed_free(p); } };
using Mat = std::unique_ptr<mzed_t,MatFree>;
struct BitFree { void operator()(mzd_t* p) const { if (p) mzd_free(p); } };
using Bits = std::unique_ptr<mzd_t,BitFree>;
static thread_local std::string last_error;

static double since(Clock::time_point& start) {
    auto now=Clock::now(); double s=std::chrono::duration<double>(now-start).count(); start=now; return s;
}

struct Field {
    gf2e* ff;
    unsigned q;
    int m;
    std::vector<unsigned> log;
    std::vector<u16> exp;
    explicit Field(unsigned modulus): ff(gf2e_init(modulus)),q(1u<<ff->degree),m(ff->degree),log(q),exp(2*(q-1)) {
        std::vector<unsigned> factors;
        unsigned left=q-1;
        for(unsigned p=2;p*p<=left;++p) if(left%p==0) {
            factors.push_back(p); while(left%p==0) left/=p;
        }
        if(left>1) factors.push_back(left);
        auto power=[&](unsigned a,unsigned e) {
            unsigned r=1;
            while(e) { if(e&1) r=gf2e_mul(ff,r,a); a=gf2e_mul(ff,a,a); e>>=1; }
            return r;
        };
        unsigned primitive=1;
        for(;primitive<q;++primitive) {
            bool good=true;
            for(auto p:factors) if(power(primitive,(q-1)/p)==1) { good=false; break; }
            if(good) break;
        }
        unsigned x=1;
        for(unsigned i=0;i<q-1;++i) { exp[i]=exp[i+q-1]=x; log[x]=i; x=gf2e_mul(ff,x,primitive); }
    }
    ~Field() { gf2e_free(ff); }
    Field(const Field&)=delete;
    Field& operator=(const Field&)=delete;
    u16 mul(u16 a,u16 b) const { return a&&b ? exp[log[a]+log[b]] : 0; }
    u16 inv(u16 a) const { if(!a) throw std::runtime_error("division by zero"); return exp[q-1-log[a]]; }
    u16 div(u16 a,u16 b) const { return mul(a,inv(b)); }
    u16 sqrt(u16 a) const { return a ? exp[(std::uint64_t(log[a])*(q/2))%(q-1)] : 0; }
};

static Mat matrix(const Field& f,int rows,int cols) { return Mat(mzed_init(f.ff,rows,cols)); }
static Mat load(const Field& f,int rows,int cols,const u16* data) {
    auto a=matrix(f,rows,cols);
    for(int i=0;i<rows;++i) for(int j=0;j<cols;++j) mzed_write_elem(a.get(),i,j,data[std::size_t(i)*cols+j]);
    return a;
}
static std::vector<int> pivots(const mzed_t* a,int rank) {
    std::vector<int> out;
    for(int i=0;i<rank;++i) {
        int j=0; while(j<a->ncols&&!mzed_read_elem(a,i,j)) ++j;
        if(j==a->ncols) throw std::runtime_error("invalid reduced matrix");
        out.push_back(j);
    }
    return out;
}
static Mat first_rows(const Field& f,const mzed_t* a,int count) {
    auto b=matrix(f,count,a->ncols);
    for(int i=0;i<count;++i) mzed_copy_row(b.get(),i,a,i);
    return b;
}
// The input is already reduced: never echelonize it a second time for a kernel.
static Mat kernel(const Field& f,const mzed_t* a,int rank) {
    auto p=pivots(a,rank); std::vector<bool> used(a->ncols);
    for(int j:p) used[j]=true;
    auto b=matrix(f,a->ncols-rank,a->ncols);
    int row=0;
    for(int j=0;j<a->ncols;++j) if(!used[j]) {
        mzed_write_elem(b.get(),row,j,1);
        for(int i=0;i<rank;++i) mzed_write_elem(b.get(),row,p[i],
            f.div(mzed_read_elem(a,i,j),mzed_read_elem(a,i,p[i])));
        ++row;
    }
    return b;
}
static std::vector<u16> derivatives(const Field& f,const std::vector<u16>& a) {
    std::vector<u16> d(a.size(),1);
    for(std::size_t i=0;i<a.size();++i) for(std::size_t j=0;j<i;++j) {
        auto x=a[i]^a[j]; if(!x) throw std::runtime_error("duplicate support");
        d[i]=f.mul(d[i],x); d[j]=f.mul(d[j],x);
    }
    return d;
}
// Return the transpose of a dual-GRS parity check, ready for multiplication.
static Mat parity_transpose(const Field& f,const std::vector<u16>& a,int dim,const std::vector<u16>& v) {
    int n=a.size(); if(dim<1||dim>=n||int(v.size())!=n) throw std::runtime_error("invalid GRS dimensions");
    auto dp=derivatives(f,a); auto h=matrix(f,n,n-dim);
    for(int j=0;j<n;++j) {
        u16 value=f.inv(f.mul(v[j],dp[j]));
        for(int r=0;r<n-dim;++r) { mzed_write_elem(h.get(),j,r,value); value=f.mul(value,a[j]); }
    }
    return h;
}
// Row space of g inside GRS_dim(a, v).
static bool contains(const Field& f,const mzed_t* g,const std::vector<u16>& a,int dim,const std::vector<u16>& v) {
    auto h=parity_transpose(f,a,dim,v);
    Mat product(mzed_mul(nullptr,g,h.get()));
    for(int i=0;i<product->nrows;++i) for(int j=0;j<product->ncols;++j)
        if(mzed_read_elem(product.get(),i,j)) return false;
    return true;
}
// A reduced matrix r equals the reduced generator of GRS_k(a, v), k = rank of r,
// on the same pivots: entry (i, j) is v_j L_i(a_j) / v_{p_i} for the Lagrange
// basis L_i on the pivot nodes. Equality with this matrix is equivalent to
// containment of rowspace(r) in GRS_k(a, v).
static bool equals_grs(const Field& f,const mzed_t* r,const std::vector<u16>& a,const std::vector<u16>& v) {
    int k=r->nrows,n=r->ncols;
    if(k<1||k>=n) return false;
    auto p=pivots(r,k);std::vector<bool> used(n);
    std::vector<u16> nodes; for(int j:p) { used[j]=true; nodes.push_back(a[j]); }
    auto dp=derivatives(f,nodes);
    for(int i=0;i<k;++i) for(int l=0;l<k;++l) if(mzed_read_elem(r,i,p[l])!=word(i==l)) return false;
    std::vector<u16> scale(k);
    for(int i=0;i<k;++i) scale[i]=f.inv(f.mul(dp[i],v[p[i]]));
    for(int j=0;j<n;++j) if(!used[j]) {
        u16 product=v[j];for(auto x:nodes) product=f.mul(product,a[j]^x);
        for(int i=0;i<k;++i)
            if(mzed_read_elem(r,i,j)!=f.mul(f.div(product,a[j]^nodes[i]),scale[i])) return false;
    }
    return true;
}
// Product of a binary matrix y (k x n) with a field matrix w (n x c, row major),
// one GF(2) product per bit of the field elements.
static std::vector<u16> binary_times(const Field& f,const mzd_t* y,const std::vector<u16>& w,int c) {
    int k=y->nrows,n=y->ncols;
    std::vector<u16> out(std::size_t(k)*c);
    Bits plane(mzd_init(n,c));
    for(int b=0;b<f.m;++b) {
        mzd_set_ui(plane.get(),0);
        for(int j=0;j<n;++j) for(int r=0;r<c;++r) if((w[std::size_t(j)*c+r]>>b)&1) mzd_write_bit(plane.get(),j,r,1);
        Bits product(mzd_mul(nullptr,y,plane.get(),0));
        for(int a=0;a<k;++a) for(int r=0;r<c;++r) if(mzd_read_bit(product.get(),a,r)) out[std::size_t(a)*c+r]|=u16(1u<<b);
    }
    return out;
}
// Binary rows y inside GRS_dim(a, v).
static bool contains_binary(const Field& f,const mzd_t* y,const std::vector<u16>& a,int dim,const std::vector<u16>& v) {
    int n=a.size(),c=n-dim; if(dim<1||dim>=n||int(v.size())!=n) throw std::runtime_error("invalid GRS dimensions");
    auto dp=derivatives(f,a); std::vector<u16> w(std::size_t(n)*c);
    for(int j=0;j<n;++j) {
        u16 value=f.inv(f.mul(v[j],dp[j]));
        for(int r=0;r<c;++r) { w[std::size_t(j)*c+r]=value; value=f.mul(value,a[j]); }
    }
    auto product=binary_times(f,y,w,c);
    return std::all_of(product.begin(),product.end(),[](u16 x) { return x==0; });
}

// Sidelnikov-Shestakov on a reduced generator: fix two coordinates to 0 and 1
// and try every value of a third one.
static std::vector<u16> recover_support(const Field& f,const mzed_t* r) {
    int k=r->nrows,n=r->ncols;
    if(k<3||n-k<2) return {};
    auto order=pivots(r,k); std::vector<bool> used(n);
    for(int j:order) used[j]=true;
    for(int j=0;j<n;++j) if(!used[j]) order.push_back(j);
    auto get=[&](int i,int j) { return u16(mzed_read_elem(r,i,order[j])); };
    for(unsigned xk=2;xk<f.q;++xk) {
        std::vector<u16> x(n); std::vector<bool> seen(f.q);
        x[0]=0;x[1]=1;x[k]=xk;seen[0]=seen[1]=seen[xk]=true;
        bool good=true;
        for(int j=k+1;j<n&&good;++j) {
            u16 denominator=f.mul(get(0,k),get(1,j));
            if(!denominator) { good=false; break; }
            u16 gamma=f.div(f.mul(get(0,j),get(1,k)),denominator);
            denominator=x[k]^f.mul(gamma,x[k]^1);
            if(!denominator) { good=false; break; }
            x[j]=f.div(x[k],denominator);
            if(seen[x[j]]) good=false; else seen[x[j]]=true;
        }
        for(int i=2;i<k&&good;++i) {
            u16 denominator=f.mul(get(i,k+1),get(0,k));
            if(!denominator) { good=false; break; }
            u16 gamma=f.div(f.mul(get(i,k),get(0,k+1)),denominator);
            denominator=f.mul(gamma,x[k+1])^x[k];
            if(!denominator) { good=false; break; }
            x[i]=f.div(f.mul(x[k],x[k+1])^f.mul(gamma,f.mul(x[k],x[k+1])),denominator);
            if(seen[x[i]]) good=false; else seen[x[i]]=true;
        }
        if(good) { std::vector<u16> out(n); for(int i=0;i<n;++i) out[order[i]]=x[i]; return out; }
    }
    return {};
}

// Multiplier of a full GRS code with support a, read from the Cauchy entries of
// its reduced generator r and checked against every row; empty if none.
static std::vector<u16> full_multiplier(const Field& f,const mzed_t* r,const std::vector<u16>& a) {
    int k=r->nrows,n=r->ncols;
    if(k<1||k>=n) return {};
    auto p=pivots(r,k);std::vector<bool> used(n);
    std::vector<u16> nodes; for(int j:p) { used[j]=true; nodes.push_back(a[j]); }
    auto dp=derivatives(f,nodes);std::vector<u16> v(n);
    int anchor=-1;u16 anchor_product=0;
    for(int j=0;j<n;++j) if(!used[j]) {
        u16 product=1;for(auto x:nodes) product=f.mul(product,a[j]^x);
        v[j]=f.div(f.mul(mzed_read_elem(r,0,j),f.mul(a[j]^nodes[0],dp[0])),product);
        if(anchor<0) { anchor=j;anchor_product=product; }
    }
    for(int i=0;i<k;++i) {
        u16 denominator=f.mul(a[anchor]^nodes[i],f.mul(dp[i],mzed_read_elem(r,i,anchor)));
        if(!denominator) return {};
        v[p[i]]=f.div(f.mul(v[anchor],anchor_product),denominator);
    }
    if(std::find(v.begin(),v.end(),0)!=v.end()||!equals_grs(f,r,a,v)) return {};
    return v;
}

// A multiplier lambda with rowspace(y) inside GRS_dim(a, lambda): solve
// y diag(z) H^T = 0 on random row combinations, then check every row.
static std::vector<u16> embedding_multiplier(const Field& f,const mzed_t* y,const std::vector<u16>& a,int dim,std::uint64_t seed) {
    int k=y->nrows,n=y->ncols,hrows=n-dim;
    auto ht=parity_transpose(f,a,dim,std::vector<u16>(n,1));
    std::mt19937_64 rng(seed);
    int used=std::min(k,std::max(8,(n+1+hrows-1)/hrows));
    while(true) {
        Mat gp;
        if(used==k) gp=Mat(mzed_copy(nullptr,y));
        else {
            auto mix=matrix(f,used,k);
            for(int i=0;i<used;++i) for(int j=0;j<k;++j) mzed_write_elem(mix.get(),i,j,rng()&(f.q-1));
            gp=Mat(mzed_mul(nullptr,mix.get(),y));
        }
        auto eq=matrix(f,used*hrows,n);
        for(int h=0;h<hrows;++h) for(int i=0;i<used;++i) for(int j=0;j<n;++j)
            mzed_write_elem(eq.get(),h*used+i,j,f.mul(mzed_read_elem(gp.get(),i,j),mzed_read_elem(ht.get(),j,h)));
        int rank=mzed_echelonize(eq.get(),1),nullity=n-rank;
        if(nullity>1&&used<k) { used=std::min(k,2*used); continue; }
        if(!nullity) return {};
        auto ker=kernel(f,eq.get(),rank);
        for(int attempt=0;attempt<nullity+200;++attempt) {
            std::vector<u16> b(n);
            if(attempt<nullity) for(int j=0;j<n;++j) b[j]=mzed_read_elem(ker.get(),attempt,j);
            else for(int i=0;i<nullity;++i) {
                u16 c=rng()&(f.q-1);
                for(int j=0;j<n;++j) b[j]^=f.mul(c,mzed_read_elem(ker.get(),i,j));
            }
            if(std::find(b.begin(),b.end(),0)!=b.end()) continue;
            for(auto& value:b) value=f.inv(value);
            if(contains(f,y,a,dim,b)) return b;
        }
        return {};
    }
}

// Minor matrix T of sampled pairs, reduced; t keeps its rank rows.
struct Minor {
    Field f;
    Mat y,t;
    int rank,rows;
    double build_seconds,rank_seconds;
    Minor(unsigned modulus,int k,int n,const u16* public_y,const u16* tangent,int count,const std::int32_t* pairs):f(modulus),y(load(f,k,n,public_y)),rows(count) {
        auto start=Clock::now();t=matrix(f,count,n);
        for(int r=0;r<count;++r) {
            int a=pairs[2*r],b=pairs[2*r+1];
            if(a<0||a>=b||b>=k) throw std::runtime_error("invalid minor pair");
            for(int j=0;j<n;++j) {
                u16 value=f.mul(public_y[std::size_t(a)*n+j],tangent[std::size_t(b)*n+j])^
                          f.mul(public_y[std::size_t(b)*n+j],tangent[std::size_t(a)*n+j]);
                mzed_write_elem(t.get(),r,j,f.sqrt(value));
            }
        }
        build_seconds=since(start);
        rank=mzed_echelonize(t.get(),1);
        rank_seconds=since(start);
        t=first_rows(f,t.get(),rank);
    }
};

// 0: no support, 1: support and minor multiplier, 2: also a public embedding.
static int recover(Minor& m,int degree,std::uint64_t seed,u16* support,u16* multiplier,double* seconds) {
    auto start=Clock::now();int n=m.t->ncols;
    Mat g=m.rank>n-m.rank ? kernel(m.f,m.t.get(),m.rank) : Mat(mzed_copy(nullptr,m.t.get()));
    mzed_echelonize(g.get(),1);
    auto a=recover_support(m.f,g.get());
    seconds[0]=since(start);
    if(a.empty()) return 0;
    auto v=full_multiplier(m.f,g.get(),a);
    seconds[1]=since(start);
    if(v.empty()) return 0;
    auto lambda=embedding_multiplier(m.f,m.y.get(),a,degree+1,seed);
    seconds[2]=since(start);
    if(lambda.empty()) return 1;
    std::copy(a.begin(),a.end(),support);std::copy(lambda.begin(),lambda.end(),multiplier);
    return 2;
}

// ---------------------------------------------------------------------------
// Binary Goppa trials

struct Random {
    std::mt19937_64 g;
    explicit Random(std::uint64_t seed): g(seed) {}
    std::uint64_t below(std::uint64_t n) {   // uniform in [0, n)
        std::uint64_t limit=std::numeric_limits<std::uint64_t>::max();
        limit-=limit%n;
        std::uint64_t x; do x=g(); while(x>=limit);
        return x%n;
    }
    u16 element(unsigned q) { return u16(g()&(q-1)); }
    u16 nonzero(unsigned q) { u16 x; do x=element(q); while(!x); return x; }
};

using Poly = std::vector<u16>;   // coefficients, constant term first

static void trim(Poly& a) { while(!a.empty()&&!a.back()) a.pop_back(); }
static u16 evaluate(const Field& f,const Poly& a,u16 x) {
    u16 y=0; for(auto c=a.rbegin();c!=a.rend();++c) y=f.mul(y,x)^*c; return y;
}
static Poly multiply(const Field& f,const Poly& a,const Poly& b) {
    Poly c(a.size()+b.size()-1);
    for(std::size_t i=0;i<a.size();++i) for(std::size_t j=0;j<b.size();++j) c[i+j]^=f.mul(a[i],b[j]);
    return c;
}
static void reduce(const Field& f,Poly& a,const Poly& g) {   // g monic
    int t=int(g.size())-1;
    for(int d=int(a.size())-1;d>=t;--d) {
        u16 c=a[d]; if(!c) continue;
        for(int i=0;i<=t;++i) a[d-t+i]^=f.mul(c,g[i]);
    }
    if(int(a.size())>t) a.resize(t);
    trim(a);
}
static int gcd_degree(const Field& f,Poly a,Poly b) {
    trim(a); trim(b);
    while(!b.empty()) {
        u16 lead=f.inv(b.back());
        for(auto& c:b) c=f.mul(c,lead);
        reduce(f,a,b); std::swap(a,b);
    }
    return int(a.size())-1;
}
// Ben-Or: g of degree t is irreducible iff gcd(g, x^(q^i) - x) = 1 for i <= t/2.
static bool irreducible(const Field& f,const Poly& g) {
    int t=int(g.size())-1;
    Poly h={0,1};
    for(int i=1;i<=t/2;++i) {
        for(int s=0;s<f.m;++s) {
            Poly square(2*h.size());
            for(std::size_t j=0;j<h.size();++j) square[2*j]=f.mul(h[j],h[j]);
            reduce(f,square,g); h=square;
        }
        Poly d=h; d.resize(std::max<std::size_t>(d.size(),2)); d[1]^=1;
        if(gcd_degree(f,g,d)>0) return false;
    }
    return true;
}
static Poly random_irreducible(const Field& f,Random& rng,int t,int& attempts) {
    while(true) {
        ++attempts;
        Poly g(t+1); for(int i=0;i<t;++i) g[i]=rng.element(f.q); g[t]=1;
        if(t==1||(g[0]&&irreducible(f,g))) return g;
    }
}
static Poly goppa_polynomial(const Field& f,Random& rng,int t,int kind,int& attempts) {
    if(kind==0) return random_irreducible(f,rng,t,attempts);
    if(kind==1) {                                     // t distinct roots
        std::set<u16> roots; Poly g={1};
        while(int(roots.size())<t) {
            u16 r=rng.element(f.q);
            if(roots.insert(r).second) g=multiply(f,g,Poly{r,1});
        }
        return g;
    }
    std::vector<int> degrees;                         // irreducible factors of distinct degrees
    if(t==5) degrees={2,3}; else if(t==6) degrees={2,4}; else if(t==7) degrees={3,4};
    else if(t==8) degrees={3,5}; else if(t==12) degrees={5,7};
    else throw std::runtime_error("no separable split for this degree");
    Poly g={1};
    for(int d:degrees) g=multiply(f,g,random_irreducible(f,rng,d,attempts));
    return g;
}
static bool squarefree(const Field& f,const Poly& g) {
    Poly d(g.size()>1 ? g.size()-1 : 1);
    for(std::size_t i=1;i<g.size();i+=2) d[i-1]=g[i];
    trim(d);
    return !d.empty()&&gcd_degree(f,g,d)==0;
}

enum { O_N, O_K, O_K0, O_D, O_BOUND, O_HRANK, O_AMBIENT, O_NOMINAL, O_SQUAREFREE,
       O_DEFECTS, O_RANK, O_ROWS, O_AVAILABLE, O_COMPLETE, O_ETA, O_FULL, O_RECOVERY,
       O_ATTEMPTS, O_COUNT };
enum { T_KEY, T_CODE, T_AMBIENT, T_TANGENTS, T_MINOR_BUILD, T_MINOR_RANK, T_CERTIFICATE,
       T_SUPPORT, T_MINOR_MULTIPLIER, T_EMBEDDING, T_COUNT };

static std::vector<std::pair<int,int>> sample_pairs(Random& rng,int k,long count) {
    std::vector<std::pair<int,int>> out;
    long available=long(k)*(k-1)/2;
    if(count>=available) {
        for(int a=0;a<k;++a) for(int b=a+1;b<k;++b) out.emplace_back(a,b);
        return out;
    }
    std::set<std::pair<int,int>> chosen;
    while(long(chosen.size())<count) {
        int a=int(rng.below(k)),b=int(rng.below(k-1));
        if(b>=a) ++b;
        chosen.emplace(std::min(a,b),std::max(a,b));
    }
    return std::vector<std::pair<int,int>>(chosen.begin(),chosen.end());
}

static void goppa_run(unsigned modulus,int n0,int t,int ell,int kind,const std::uint64_t* seeds,int want_recovery,
                      int* out,double* seconds,u16* alpha0_out,int* keep_out,u16* g_out,std::uint8_t* y_out,
                      int* defect_out,u16* support_out,u16* multiplier_out,u16* tangent_out) {
    Field f(modulus);
    std::fill(out,out+O_COUNT,0); std::fill(seconds,seconds+T_COUNT,0.0);
    out[O_RECOVERY]=-1;
    auto start=Clock::now();

    // Key: Goppa polynomial and n0 distinct support elements outside its roots.
    Random key(seeds[0]);
    int attempts=0;
    Poly g=goppa_polynomial(f,key,t,kind,attempts);
    out[O_ATTEMPTS]=attempts;
    out[O_SQUAREFREE]=squarefree(f,g);
    std::vector<u16> admissible;
    for(unsigned x=0;x<f.q;++x) if(evaluate(f,g,u16(x))) admissible.push_back(u16(x));
    if(int(admissible.size())<n0) throw std::runtime_error("not enough support elements");
    for(std::size_t i=admissible.size()-1;i>0;--i) std::swap(admissible[i],admissible[key.below(i+1)]);
    std::vector<u16> alpha0(admissible.begin(),admissible.begin()+n0);
    std::copy(alpha0.begin(),alpha0.end(),alpha0_out); std::copy(g.begin(),g.end(),g_out);
    seconds[T_KEY]=since(start);

    // Code: binary expansion of alpha_j^i / G(alpha_j), 0 <= i < t, in reduced form.
    int m=f.m;
    Bits h(mzd_init(m*t,n0));
    for(int j=0;j<n0;++j) {
        u16 value=f.inv(evaluate(f,g,alpha0[j]));
        for(int i=0;i<t;++i) {
            for(int b=0;b<m;++b) if((value>>b)&1) mzd_write_bit(h.get(),i*m+b,j,1);
            value=f.mul(value,alpha0[j]);
        }
    }
    int hrank=mzd_echelonize(h.get(),1);
    std::vector<int> pivot; std::vector<bool> is_pivot(n0);
    for(int r=0,j=0;r<hrank;++r) {
        while(!mzd_read_bit(h.get(),r,j)) ++j;
        pivot.push_back(j); is_pivot[j]=true;
    }
    // Shorten at the first ell non-pivot positions (an information set of the code).
    std::vector<bool> dropped(n0); int left=ell;
    for(int j=0;j<n0&&left>0;++j) if(!is_pivot[j]) { dropped[j]=true; --left; }
    if(left) throw std::runtime_error("shortening exceeds the dimension");
    std::vector<int> keep, position(n0,-1);
    for(int j=0;j<n0;++j) if(!dropped[j]) { position[j]=int(keep.size()); keep.push_back(j); }
    int n=int(keep.size()),k=n0-hrank-ell;
    Bits yb(mzd_init(k,n));
    for(int j=0,row=0;j<n0;++j) if(!dropped[j]&&!is_pivot[j]) {
        mzd_write_bit(yb.get(),row,position[j],1);
        for(int r=0;r<hrank;++r) if(mzd_read_bit(h.get(),r,j)) mzd_write_bit(yb.get(),row,position[pivot[r]],1);
        ++row;
    }
    mzd_echelonize(yb.get(),1);
    std::vector<u16> alpha(n),y(std::size_t(k)*n);
    for(int j=0;j<n;++j) { alpha[j]=alpha0[keep[j]]; keep_out[j]=keep[j]; }
    for(int a=0;a<k;++a) for(int j=0;j<n;++j) y[std::size_t(a)*n+j]=y_out[std::size_t(a)*n+j]=mzd_read_bit(yb.get(),a,j);
    int D=n-2*t-1,bound=D-t;
    out[O_N]=n; out[O_K]=k; out[O_K0]=n0-hrank; out[O_D]=D; out[O_BOUND]=bound; out[O_HRANK]=hrank;
    seconds[T_CODE]=since(start);

    // Ambient GRS_{D+1}(alpha, lambda), lambda_j = G(alpha_j)^2 / Pi'(alpha_j).
    auto dp=derivatives(f,alpha);
    std::vector<u16> galpha(n),lambda(n),nu(n);
    for(int j=0;j<n;++j) {
        galpha[j]=evaluate(f,g,alpha[j]);
        lambda[j]=f.div(f.mul(galpha[j],galpha[j]),dp[j]);
        nu[j]=f.inv(f.mul(galpha[j],galpha[j]));
    }
    out[O_AMBIENT]=contains_binary(f,yb.get(),alpha,D+1,lambda);
    out[O_NOMINAL]=!contains_binary(f,yb.get(),alpha,D,lambda);
    seconds[T_AMBIENT]=since(start);

    // Tangent representatives yhat_j = sum_{l != j} nu_l y_l / (alpha_j - alpha_l),
    // then yhat_j <- a_j yhat_j + b_j y_j with random a_j != 0 and b_j.
    std::vector<u16> cauchy(std::size_t(n)*n);
    for(int l=0;l<n;++l) for(int j=0;j<n;++j) if(l!=j) cauchy[std::size_t(l)*n+j]=f.div(nu[l],alpha[j]^alpha[l]);
    auto tangent=binary_times(f,yb.get(),cauchy,n);
    std::vector<u16>().swap(cauchy);
    Random representatives(seeds[1]);
    std::vector<u16> scale(n),shift(n);
    for(int j=0;j<n;++j) scale[j]=representatives.nonzero(f.q);
    for(int j=0;j<n;++j) shift[j]=representatives.element(f.q);
    std::vector<u16> yhat(std::size_t(k)*n);
    for(int a=0;a<k;++a) for(int j=0;j<n;++j)
        yhat[std::size_t(a)*n+j]=f.mul(tangent[std::size_t(a)*n+j],scale[j])^(y[std::size_t(a)*n+j] ? shift[j] : 0);
    if(tangent_out) std::copy(yhat.begin(),yhat.end(),tangent_out);
    int defects=0;
    for(int j=0;j<n;++j) {
        int a0=0; while(a0<k&&!y[std::size_t(a0)*n+j]) ++a0;
        bool plane=false;
        if(a0<k) {
            u16 c=yhat[std::size_t(a0)*n+j];
            for(int a=0;a<k&&!plane;++a) plane=yhat[std::size_t(a)*n+j]!=(y[std::size_t(a)*n+j] ? c : 0);
        }
        if(!plane) defect_out[defects++]=j;
    }
    out[O_DEFECTS]=defects;
    seconds[T_TANGENTS]=since(start);

    // Minor code: sample pairs, enlarge until the bound D - t is attained or all
    // pairs are used. Every sampled row is checked against GRS_{D-t}(alpha, eta),
    // eta_j = G(alpha_j)^2 sqrt(a_j) / Pi'(alpha_j).
    Random sampler(seeds[2]);
    long available=long(k)*(k-1)/2;
    long budget=std::min<long>(available,bound+64);
    std::vector<u16> eta(n);
    for(int j=0;j<n;++j) eta[j]=f.div(f.mul(f.mul(galpha[j],galpha[j]),f.sqrt(scale[j])),dp[j]);
    std::unique_ptr<Minor> minor;
    while(true) {
        auto pairs=sample_pairs(sampler,k,budget);
        std::vector<std::int32_t> flat;
        for(auto& p:pairs) { flat.push_back(p.first); flat.push_back(p.second); }
        minor.reset(new Minor(modulus,k,n,y.data(),yhat.data(),int(pairs.size()),flat.data()));
        seconds[T_MINOR_BUILD]+=minor->build_seconds; seconds[T_MINOR_RANK]+=minor->rank_seconds;
        start=Clock::now();
        out[O_ETA]=minor->rank==bound ? equals_grs(minor->f,minor->t.get(),alpha,eta)
                  : minor->rank==0||contains(minor->f,minor->t.get(),alpha,bound,eta);
        seconds[T_CERTIFICATE]+=since(start);
        if(!out[O_ETA]) throw std::runtime_error("sampled minor rows violate the predicted containment");
        if(minor->rank==bound||budget==available) break;
        budget=std::min(available,2*budget);
    }
    out[O_RANK]=minor->rank; out[O_ROWS]=minor->rows;
    out[O_AVAILABLE]=int(std::min<long>(available,std::numeric_limits<int>::max()));
    out[O_COMPLETE]=1;
    if(minor->rank==bound) out[O_FULL]=1;
    else if(minor->rank>0&&!defects) {
        start=Clock::now();
        out[O_FULL]=!full_multiplier(minor->f,minor->t.get(),alpha).empty();
        seconds[T_CERTIFICATE]+=since(start);
    }

    // Recovery from public data: the minor matrix, y and D.
    if(want_recovery&&out[O_FULL])
        out[O_RECOVERY]=recover(*minor,D,seeds[3],support_out,multiplier_out,seconds+T_SUPPORT);
}

extern "C" {
const char* minor_error() { return last_error.c_str(); }
void* minor_create(unsigned modulus,int k,int n,const u16* y,const u16* tangent,int count,const std::int32_t* pairs) {
    try { return new Minor(modulus,k,n,y,tangent,count,pairs); }
    catch(const std::exception& e) { last_error=e.what();return nullptr; }
}
void minor_free(void* handle) { delete static_cast<Minor*>(handle); }
int minor_rank(void* handle) { return static_cast<Minor*>(handle)->rank; }
void minor_times(void* handle,double* out) {
    auto& m=*static_cast<Minor*>(handle);out[0]=m.build_seconds;out[1]=m.rank_seconds;
}
// Verification does not retain its secret inputs or change the public state.
int minor_verify(void* handle,const u16* support,const u16* multiplier,int dim) {
    try {
        auto& m=*static_cast<Minor*>(handle);int n=m.t->ncols;
        return contains(m.f,m.t.get(),std::vector<u16>(support,support+n),dim,std::vector<u16>(multiplier,multiplier+n));
    } catch(const std::exception& e) { last_error=e.what();return -1; }
}
int minor_recover(void* handle,int degree,std::uint64_t seed,u16* support,u16* multiplier,double* seconds) {
    try { return recover(*static_cast<Minor*>(handle),degree,seed,support,multiplier,seconds); }
    catch(const std::exception& e) { last_error=e.what();return -1; }
}
int goppa_trial(unsigned modulus,int n0,int t,int ell,int kind,const std::uint64_t* seeds,int want_recovery,
                int* out,double* seconds,u16* alpha0,int* keep,u16* g,std::uint8_t* y,int* defects,
                u16* support,u16* multiplier,u16* tangent) {
    try {
        goppa_run(modulus,n0,t,ell,kind,seeds,want_recovery,out,seconds,alpha0,keep,g,y,defects,support,multiplier,tangent);
        return 0;
    } catch(const std::exception& e) { last_error=e.what();return -1; }
}
}
