#include "algebra.hpp"
#include "series_api.hpp"
#include <m4ri/m4ri.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>

using namespace toy;
namespace fs=std::filesystem;
using Clock=std::chrono::steady_clock;
static double elapsed(Clock::time_point t) { return std::chrono::duration<double>(Clock::now()-t).count(); }
struct BinaryFree { void operator()(mzd_t* p) const { if(p)mzd_free(p); } };
using Binary=std::unique_ptr<mzd_t,BinaryFree>;
static Binary binary(int r,int c) { require(r>0 && c>0,"empty packed matrix");return Binary(mzd_init(r,c)); }
static void array_json(std::ostream& out,const Vec& v) {
    out<<'[';for(std::size_t i=0;i<v.size();++i){if(i)out<<',';out<<int(v[i]);}out<<']';
}
struct Report {
    std::map<std::string,double> numbers;std::map<std::string,bool> checks;
    void write(std::ostream& out) const {
        out<<'{';bool comma=false;
        for(const auto& kv:numbers){if(comma)out<<',';comma=true;out<<'"'<<kv.first<<"\":"<<kv.second;}
        for(const auto& kv:checks){if(comma)out<<',';comma=true;out<<'"'<<kv.first<<"\":"<<(kv.second?"true":"false");}out<<'}';
    }
};
struct Public {
    int m,t,n,k,modulus,held,d,s;std::vector<std::uint64_t> points;
    int D() const { return n-2*t-1; }
    Vec point(int j) const { Vec v(k);for(int a=0;a<k;++a)v[a]=(points[j]>>a)&1;return v; }
    Matrix generator() const { Matrix y(k,n);for(int j=0;j<n;++j)for(int a=0;a<k;++a)y(a,j)=(points[j]>>a)&1;return y; }
    std::vector<std::uint64_t> rows() const {
        std::vector<std::uint64_t> a(k);for(int i=0;i<k;++i)for(int j=0;j<n;++j)a[i]|=((points[j]>>i)&1)<<j;return a;
    }
};
static Public read_public(const fs::path& path) {
    std::ifstream in(path);std::string magic;Public p{};
    require(bool(in>>magic) && magic=="BGPK1","invalid public format");
    require(bool(in>>p.m>>p.t>>p.n>>p.k>>p.modulus>>p.held>>p.d>>p.s),"public header");
    require(p.m==6 && p.t==6 && p.n==64 && p.k==28 && p.d==5 && p.s==4,
            "this implementation targets the unshortened (6,6,64,28) degree-five toy");
    require(p.held>=0 && p.held<p.n,"held-out index");
    p.points.resize(p.n);for(auto& x:p.points)require(bool(in>>x) && !(x>>p.k),"invalid public column");
    std::string extra;require(!(in>>extra),"trailing public data");
    require(int(binary_rref(p.rows(),p.n).size())==p.k,"public generator rank");
    Field f(p.m,p.modulus);(void)f;return p;
}
static void write_public(const fs::path& path,const Public& p) {
    std::ofstream out(path);out<<"BGPK1\n"<<p.m<<' '<<p.t<<' '<<p.n<<' '<<p.k<<' '
        <<p.modulus<<' '<<p.held<<' '<<p.d<<' '<<p.s<<'\n';
    for(auto x:p.points)out<<x<<' ';out<<'\n';require(bool(out),"write public failed");
}
static void prepare(const Public& p,const fs::path& output) {
    std::vector<std::uint64_t> points;std::uint64_t units=0;
    for(int j=0;j<p.n;++j)if(j!=p.held) {
        auto v=p.points[j];if(__builtin_popcountll(v)==1)units|=v;else points.push_back(v);
    }
    require(units==((std::uint64_t(1)<<p.k)-1),"retained information set must contain every unit column");
    std::ofstream out(output);out<<p.k<<' '<<p.d<<' '<<points.size()<<'\n';
    for(auto point:points)out<<point<<' '<<p.s<<'\n';require(bool(out),"write operator failed");
}
static void save_matrix(const fs::path& path,const Matrix& a) {
    std::ofstream out(path);out<<"GF64M1\n"<<a.rows<<' '<<a.cols<<'\n';
    for(int r=0;r<a.rows;++r){for(int c=0;c<a.cols;++c)out<<int(a(r,c))<<' ';out<<'\n';}
    require(bool(out),"write matrix failed");
}
static std::uint32_t read_u32(std::istream& in) {
    unsigned char b[4];require(bool(in.read(reinterpret_cast<char*>(b),4)),"truncated kernel header");
    return std::uint32_t(b[0])|(std::uint32_t(b[1])<<8)|(std::uint32_t(b[2])<<16)|(std::uint32_t(b[3])<<24);
}
static Binary read_kernel(const Public& p,const fs::path& path) {
    std::ifstream in(path,std::ios::binary);char magic[4];
    require(bool(in.read(magic,4)) && std::string(magic,4)=="WHK1","invalid holdout kernel");
    auto k=read_u32(in),d=read_u32(in),N=read_u32(in),M=read_u32(in);
    require(k==unsigned(p.k) && d==unsigned(p.d) && N==masks(p.k,p.d).size(),"kernel parameter mismatch");
    require(M>0 && M<=N,"invalid kernel dimension");
    auto K=binary(M,N);Vec packed((N+7)/8);
    for(unsigned r=0;r<M;++r) {
        require(bool(in.read(reinterpret_cast<char*>(packed.data()),packed.size())),"truncated kernel");
        for(unsigned c=0;c<N;++c)if((packed[c/8]>>(c%8))&1)mzd_write_bit(K.get(),r,c,1);
    }
    char extra;require(!in.get(extra),"trailing kernel data");
    Binary copy(mzd_copy(nullptr,K.get()));require(mzd_echelonize_m4ri(copy.get(),0,0)==int(M),"dependent kernel vectors");
    return K;
}
static void verify_hasse(const Public& p,const mzd_t* K) {
    auto mons=masks(p.k,p.d);std::unordered_map<std::uint64_t,int> index;
    for(int c=0;c<int(mons.size());++c)index[mons[c]]=c;
    std::vector<std::uint64_t> derivatives;
    for(int u=0;u<p.s;++u){auto v=masks(p.k,u);derivatives.insert(derivatives.end(),v.begin(),v.end());}
    Binary KT(mzd_transpose(nullptr,K));
    for(int j=0;j<p.n;++j)if(j!=p.held && __builtin_popcountll(p.points[j])!=1) {
        auto A=binary(int(derivatives.size()),int(mons.size()));
        for(int r=0;r<int(derivatives.size());++r) {
            auto s=derivatives[r];int degree=p.d-__builtin_popcountll(s);
            for(auto extra:masks(p.k,degree,p.points[j]&~s))mzd_write_bit(A.get(),r,index.at(s|extra),1);
        }
        Binary residual(mzd_mul(nullptr,A.get(),KT.get(),0));require(mzd_is_zero(residual.get()),"a returned polynomial violates a retained Hasse condition");
    }
}
static Matrix apply_binary(const Field& f,const mzd_t* K,const Matrix& values) {
    require(K->ncols==values.rows && values.cols>0,"binary product dimensions");
    auto B=binary(values.rows,values.cols*f.m);
    for(int r=0;r<values.rows;++r)for(int c=0;c<values.cols;++c) {
        Byte x=values(r,c);for(int b=0;b<f.m;++b)if((x>>b)&1)mzd_write_bit(B.get(),r,c*f.m+b,1);
    }
    Binary C(mzd_mul(nullptr,K,B.get(),0));Matrix result(K->nrows,values.cols);
    for(int r=0;r<result.rows;++r)for(int c=0;c<result.cols;++c)
        for(int b=0;b<f.m;++b)result(r,c)|=mzd_read_bit(C.get(),r,c*f.m+b)<<b;
    return result;
}
struct Local { Matrix J,Q,B;std::vector<std::pair<int,int>> pairs;std::vector<Vec> branches;Report report; };
static Matrix polar(const Field& f,const Local& local,const Vec& v) {
    Matrix P(local.Q.cols,int(v.size()));
    for(int r=0;r<int(local.pairs.size());++r){auto [a,b]=local.pairs[r];P(r,a)=v[b];P(r,b)=v[a];}
    return multiply(f,local.Q,P);
}
static Local recover_branches(const Public& p,const Field& f,const mzd_t* K) {
    Local local;for(int a=0;a<p.k;++a)for(int b=a+1;b<p.k;++b)local.pairs.emplace_back(a,b);
    int count=1+p.k+int(local.pairs.size()),width=(K->ncols+7)/8;Vec bytes(std::size_t(count)*width);
    require(local_rows(p.k,p.d,p.points[p.held],bytes.data())==0,"Taylor-map construction");
    auto A=binary(count,K->ncols);
    for(int r=0;r<count;++r)for(int c=0;c<K->ncols;++c)
        if((bytes[std::size_t(r)*width+c/8]>>(c%8))&1)mzd_write_bit(A.get(),r,c,1);
    Binary AT(mzd_transpose(nullptr,A.get())),maps(mzd_mul(nullptr,K,AT.get(),0));
    local.J=Matrix(K->nrows,p.k);local.Q=Matrix(K->nrows,int(local.pairs.size()));
    for(int r=0;r<K->nrows;++r) {
        require(!mzd_read_bit(maps.get(),r,0),"holdout polynomial does not vanish at the seed");
        for(int c=0;c<p.k;++c)local.J(r,c)=mzd_read_bit(maps.get(),r,1+c);
        for(int c=0;c<local.Q.cols;++c)local.Q(r,c)=mzd_read_bit(maps.get(),r,1+p.k+c);
    }
    Matrix W=nullspace(f,local.J);require(W.rows==p.m+1,"gradient nullity is not m+1");
    std::vector<Vec> basis={p.point(p.held)};
    for(int r=0;r<W.rows;++r) {
        auto trial=basis;trial.push_back(W.row(r));if(rank(f,from_rows(trial,p.k))>int(basis.size()))basis.push_back(W.row(r));
    }
    require(int(basis.size())==p.m+1,"gradient quotient basis");basis.erase(basis.begin());local.B=from_rows(basis,p.k);
    std::vector<std::pair<int,int>> qm;for(int a=0;a<p.m;++a)for(int b=a;b<p.m;++b)qm.emplace_back(a,b);
    Matrix substitution(local.Q.cols,int(qm.size()));
    for(int r=0;r<substitution.rows;++r)for(int c=0;c<substitution.cols;++c) {
        auto [a,b]=local.pairs[r];auto [u,v]=qm[c];
        substitution(r,c)=f.mul(local.B(u,a),local.B(v,b));
        if(u!=v)substitution(r,c)^=f.mul(local.B(v,a),local.B(u,b));
    }
    Matrix restrictions=multiply(f,local.Q,substitution);
    require(rank(f,restrictions)==p.m*(p.m-1)/2,"quadratic separation rank criterion");
    Matrix quotient=nullspace(f,restrictions);require(quotient.rows==p.m,"quadratic quotient dimension");
    std::vector<Matrix> multiplication;
    for(int a=0;a<p.m;++a) {
        Matrix M(p.m,p.m);
        for(int j=0;j<p.m;++j) {
            auto target=std::make_pair(std::min(a,j),std::max(a,j));
            int c=int(std::find(qm.begin(),qm.end(),target)-qm.begin());
            for(int r=0;r<p.m;++r)M(r,j)=quotient(r,c);
        }
        multiplication.push_back(M);
    }
    Matrix A1=multiplication[0];require(rank(f,A1)==p.m,"chosen projective chart misses a branch");
    for(int mask=2;mask<(1<<p.m) && local.branches.empty();++mask) {
        Matrix A2(p.m,p.m);for(int a=0;a<p.m;++a)if((mask>>a)&1)
            for(std::size_t c=0;c<A2.data.size();++c)A2.data[c]^=multiplication[a].data[c];
        std::vector<Vec> found;
        for(int theta=0;theta<f.q;++theta) {
            Matrix pencil=A2;for(std::size_t c=0;c<pencil.data.size();++c)pencil.data[c]^=f.mul(Byte(theta),A1.data[c]);
            Matrix left=nullspace(f,pencil.transpose());if(left.rows!=1)continue;
            Vec coordinates=multiply(f,A1.transpose(),left.row(0));found.push_back(multiply(f,local.B.transpose(),coordinates));
        }
        if(int(found.size())==p.m)local.branches=found;
    }
    require(int(local.branches.size())==p.m,"could not separate all Frobenius branches");
    for(const Vec& v:local.branches) {
        auto value=multiply(f,local.J,v);require(std::all_of(value.begin(),value.end(),[](Byte x){return !x;}),"branch gradient residual");
        Vec products;for(auto [a,b]:local.pairs)products.push_back(f.mul(v[a],v[b]));
        auto residual=multiply(f,local.Q,products);require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"branch quadratic residual");
        require(rank(f,stack(local.J,polar(f,local,v)))==p.k-2,"binary continuation rank criterion on a branch");
    }
    local.report.numbers={{"gradient_rank",double(p.k-p.m-1)},{"gradient_nullity",double(p.m+1)},
        {"quadratic_rank",double(p.m*(p.m-1)/2)},{"branches",double(local.branches.size())},{"continuation_rank_all_branches",double(p.k-2)}};
    return local;
}

struct Jet { Matrix coefficients;Report report; };
static Jet continue_jet(const Public& p,const Field& f,const mzd_t* K,const Local& local,int depth) {
    require(depth>=2 && depth%2==0,"normalized continuation requests an even final order");
    auto started=Clock::now();Matrix A=stack(local.J,polar(f,local,local.branches[0]));FixedSolver solve(f,A);
    Matrix jets(depth+1,p.k);jets.set_row(0,p.point(p.held));jets.set_row(1,local.branches[0]);
    Vec second=local.branches[0];for(auto& x:second)x=f.sq(x);jets.set_row(2,second);
    std::unique_ptr<void,decltype(&series_cache_free)> cache(series_cache_create(p.k,p.d,p.m,p.modulus,depth),series_cache_free);
    require(bool(cache),"series cache allocation");Matrix values(K->ncols,1);
    auto coefficient=[&](int r) {
        require(series_cache_eval(cache.get(),r,values.data.data())==0,"series coefficient evaluation");
        return apply_binary(f,K,values).column(0);
    };
    for(int r=0;r<=2;++r) {
        auto v=jets.row(r);require(series_cache_set(cache.get(),r,v.data())==0,"set jet coefficient");
        auto residual=coefficient(r);require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"initial jet substitution");
    }
    int first_full=-1;
    for(int r=3;r<depth;r+=2) {
        Vec even=jets.row((r+1)/2);for(auto& x:even)x=f.sq(x);
        Vec rhs=coefficient(r),next=coefficient(r+1),Jnext=multiply(f,local.J,even);
        for(int j=0;j<int(next.size());++j)next[j]^=Jnext[j];rhs.insert(rhs.end(),next.begin(),next.end());
        Vec odd=solve.solve(rhs);jets.set_row(r,odd);jets.set_row(r+1,even);
        require(series_cache_set(cache.get(),r,odd.data())==0 && series_cache_set(cache.get(),r+1,even.data())==0,"update jet");
        for(int u=r;u<=r+1;++u) {
            auto residual=coefficient(u);require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"continued jet substitution");
            if(first_full<0 && rank(f,jets.first_rows(u+1))==p.k)first_full=u;
        }
    }
    double recurrence_seconds=elapsed(started);auto tick=Clock::now();Matrix arcs=jets.transpose(),fresh(K->ncols,depth+1);
    require(series_products(p.k,p.d,p.m,p.modulus,depth,arcs.data.data(),fresh.data.data())==0,"fresh full product computation");
    require(apply_binary(f,K,fresh).zero(),"independent full jet substitution failed");
    for(int r=1;2*r<=depth;++r)for(int a=0;a<p.k;++a)require(jets(2*r,a)==f.sq(jets(r,a)),"binary normalization failed");
    Report report;report.numbers={{"depth",double(depth)},{"continuation_rank",double(solve.columns.size())},
        {"first_full_coefficient_span_order",double(first_full)},{"recurrence_seconds",recurrence_seconds},
        {"fresh_substitution_seconds",elapsed(tick)}};
    report.checks={{"normalized",true},{"full_substitution",true},{"all_polynomials_used",true}};
    return {jets,report};
}

using Pairs=std::vector<std::pair<int,int>>;
static Pairs quadric_monomials(int k) {
    Pairs pairs;for(int a=0;a<k;++a)for(int b=a;b<k;++b)pairs.emplace_back(a,b);return pairs;
}
static Vec quadratic_values(const Field& f,const Vec& point,const Pairs& pairs) {
    Vec v;for(auto [a,b]:pairs)v.push_back(f.mul(point[a],point[b]));return v;
}
static Matrix contacts(const Field& f,const Matrix& jets,const Pairs& pairs,int order) {
    require(jets.rows>order,"insufficient contact precision");Matrix a(order+1,int(pairs.size()));
    for(int r=0;r<=order;++r)for(int c=0;c<a.cols;++c) {
        auto [u,v]=pairs[c];for(int e=0;e<=r;++e)a(r,c)^=f.mul(jets(e,u),jets(r-e,v));
    }
    return a;
}
struct Bootstrap { Matrix forms,tangents;Report report; };
static Bootstrap bootstrap(const Public& p,const Field& f,const Matrix& jets) {
    auto start=Clock::now();int depth=std::max(0,2*p.D()-p.n+1);Pairs pairs=quadric_monomials(p.k);
    Matrix A(p.n+depth,int(pairs.size()));
    for(int j=0;j<p.n;++j)A.set_row(j,quadratic_values(f,p.point(j),pairs));
    Matrix C=contacts(f,jets,pairs,depth);for(int r=1;r<=depth;++r)A.set_row(p.n+r-1,C.row(r));
    Matrix H=nullspace(f,A),tangents(p.k,p.n);
    require(multiply(f,A,H.transpose()).zero(),"bootstrap linear residual");
    for(int j=0;j<p.n;++j) {
        Vec point=p.point(j);Matrix J(H.rows,p.k);
        for(int r=0;r<H.rows;++r)for(int c=0;c<H.cols;++c) {
            auto [a,b]=pairs[c];J(r,a)^=f.mul(H(r,c),point[b]);J(r,b)^=f.mul(H(r,c),point[a]);
        }
        Matrix W=nullspace(f,J);require(W.rows==2,"a bootstrapped tangent kernel is not two-dimensional");
        bool found=false;
        for(int r=0;r<W.rows;++r)if(rank(f,from_rows({point,W.row(r)},p.k))==2) {
            for(int a=0;a<p.k;++a)tangents(a,j)=W(r,a);found=true;break;
        }
        require(found,"could not choose a tangent representative");
    }
    Report report;report.numbers={{"jet_order",double(depth)},{"quadratic_unknowns",double(pairs.size())},
        {"constraint_rank",double(A.cols-H.rows)},{"ideal_dimension",double(H.rows)},
        {"gradient_rank_all_positions",double(p.k-2)},{"positions",double(p.n)},{"seconds",elapsed(start)}};
    report.checks={{"ordinary_bootstrap_bound",p.n+depth>2*p.D()}};
    return {H,tangents,report};
}
static Matrix grs(const Field& f,const Vec& support,int dimension,const Vec& weights) {
    require(support.size()==weights.size(),"GRS weights length");Matrix A(dimension,int(support.size()));
    for(int j=0;j<A.cols;++j){Byte x=weights[j];for(int r=0;r<dimension;++r){A(r,j)=x;x=f.mul(x,support[j]);}}return A;
}
static Matrix parity(const Field& f,const Vec& support,int dimension,const Vec& weights) {
    Vec dp=support_derivatives(f,support),dual(support.size());
    for(int j=0;j<int(support.size());++j)dual[j]=f.inv(f.mul(weights[j],dp[j]));
    return grs(f,support,int(support.size())-dimension,dual);
}
static Vec recover_multiplier(const Field& f,const Matrix& Y,const Vec& support,int dimension) {
    Matrix H=parity(f,support,dimension,Vec(support.size(),1)),A(Y.rows*H.rows,Y.cols);
    for(int a=0;a<Y.rows;++a)for(int b=0;b<H.rows;++b)for(int j=0;j<Y.cols;++j)A(a*H.rows+b,j)=f.mul(Y(a,j),H(b,j));
    Matrix K=nullspace(f,A);require(K.rows==1,"multiplier space is not one-dimensional");
    Vec lambda=K.row(0);for(auto& x:lambda)x=f.inv(x);
    require(multiply(f,Y,parity(f,support,dimension,lambda).transpose()).zero(),"public GRS containment failed");return lambda;
}
static Vec recover_grs_support(const Field& f,const Matrix& input) {
    Matrix G=row_basis(f,input);int n=G.cols;
    if(n-G.rows<G.rows)G=nullspace(f,G);
    auto pivots=rref(f,G);int k=G.rows;require(k>=2 && k<=n-2,"Sidelnikov-Shestakov dimensions");
    std::vector<int> order=pivots;for(int j=0;j<n;++j)if(std::find(pivots.begin(),pivots.end(),j)==pivots.end())order.push_back(j);
    auto get=[&](int i,int j){return G(i,order[j]);};
    for(int i=0;i<k;++i)for(int j=k;j<n;++j)require(get(i,j)!=0,"zero entry in systematic GRS block");
    for(int trial=2;trial<f.q;++trial) {
        Vec a(n);a[1]=1;a[k]=Byte(trial);bool good=true;
        for(int j=k+1;j<n && good;++j) {
            Byte ratio=f.div(f.mul(get(0,j),get(1,k)),f.mul(get(1,j),get(0,k)));
            Byte denominator=a[k]^f.mul(ratio,a[k]^1);if(!denominator){good=false;break;}
            a[j]=f.div(a[k],denominator);
        }
        for(int i=2;i<k && good;++i) {
            Byte ratio=f.div(f.mul(get(i,k),get(0,k+1)),f.mul(get(0,k),get(i,k+1)));
            Byte denominator=f.mul(ratio,a[k+1])^a[k];if(!denominator){good=false;break;}
            a[i]=f.div(f.mul(f.mul(a[k],a[k+1]),ratio^1),denominator);
        }
        if(!good)continue;
        Vec sorted=a;std::sort(sorted.begin(),sorted.end());if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())continue;
        Vec support(n);for(int j=0;j<n;++j)support[order[j]]=a[j];
        try {
            Vec lambda=recover_multiplier(f,G,support,k);
            require(same_span(f,G,grs(f,support,k,lambda)),"GRS equality");return support;
        } catch(const std::exception&) { /* A trial with no full GRS description is discarded. */ }
    }
    throw std::runtime_error("Sidelnikov-Shestakov failed to find a finite support");
}
struct Key { Vec support,lambda,g;Report report; };
static std::vector<std::uint64_t> expand_binary(const Field& f,const Matrix& H) {
    require(H.cols<=64,"binary check width");std::vector<std::uint64_t> rows(H.rows*f.m);
    for(int r=0;r<H.rows;++r)for(int j=0;j<H.cols;++j)for(int b=0;b<f.m;++b)
        rows[r*f.m+b]|=std::uint64_t((H(r,j)>>b)&1)<<j;
    return binary_rref(rows,H.cols);
}
static Report verify_key(const Public& p,const Field& f,const Key& key) {
    require(int(key.support.size())==p.n && int(key.lambda.size())==p.n,"key lengths");
    for(Byte a:key.support)require(a<f.q,"support field encoding");
    for(Byte a:key.lambda)require(a>0 && a<f.q,"multiplier field encoding");
    auto target=binary_rref(binary_kernel(p.rows(),p.n),p.n);
    Matrix H=parity(f,key.support,p.D()+1,key.lambda);
    require(multiply(f,p.generator(),H.transpose()).zero(),"recovered key does not contain the public code");
    auto actual=expand_binary(f,H);require(actual==target,"recovered alternant code differs from public code");
    require(int(key.g.size())==p.t+1 && key.g.back()==1 && irreducible(f,key.g),"recovered Goppa polynomial is not monic irreducible degree six");
    Vec beta(p.n);for(int j=0;j<p.n;++j)beta[j]=f.inv(evaluate(f,key.g,key.support[j]));
    require(expand_binary(f,grs(f,key.support,p.t,beta))==target,"recovered Goppa code differs from public code");
    Vec dp=support_derivatives(f,key.support);Byte scale=0;
    for(int j=0;j<p.n;++j) {
        Byte value=f.div(f.mul(key.lambda[j],dp[j]),f.sq(evaluate(f,key.g,key.support[j])));
        if(!j)scale=value;require(value==scale,"Goppa multiplier square identity");
    }
    Report report;report.numbers={{"support_count",double(p.n)},{"binary_parity_rank",double(target.size())},{"goppa_degree",double(p.t)}};
    report.checks={{"distinct_finite_support",true},{"exact_public_alternant_equality",true},
        {"exact_public_goppa_equality",true},{"goppa_irreducible",true},{"multiplier_square_identity",true}};return report;
}
static void complete_key(const Public& p,const Field& f,Key& key) {
    key.lambda=recover_multiplier(f,p.generator(),key.support,p.D()+1);
    Vec dp=support_derivatives(f,key.support),values(p.n);
    for(int j=0;j<p.n;++j)values[j]=f.root(f.mul(key.lambda[j],dp[j]));
    key.g=interpolate(f,key.support,values);require(int(key.g.size())==p.t+1,"Goppa polynomial interpolation degree");
    key.g=scale(f,key.g,f.inv(key.g.back()));Report verified=verify_key(p,f,key);
    key.report.numbers.insert(verified.numbers.begin(),verified.numbers.end());
    key.report.checks.insert(verified.checks.begin(),verified.checks.end());
}
static Key finish_minor(const Public& p,const Field& f,const Bootstrap& bootstrap,Matrix& minors) {
    auto start=Clock::now();minors=Matrix(p.k*(p.k-1)/2,p.n);int row=0;Matrix Y=p.generator();
    for(int a=0;a<p.k;++a)for(int b=a+1;b<p.k;++b) {
        for(int j=0;j<p.n;++j)minors(row,j)=f.root(f.mul(Y(a,j),bootstrap.tangents(b,j))^f.mul(Y(b,j),bootstrap.tangents(a,j)));++row;
    }
    int dimension=rank(f,minors);require(dimension==p.D()-p.t,"minor code is not full at the Goppa dimension bound");
    Key key;key.support=recover_grs_support(f,minors);
    Vec lambda=recover_multiplier(f,minors,key.support,dimension);
    require(same_span(f,minors,grs(f,key.support,dimension,lambda)),"complete minor code is not the recovered GRS code");
    complete_key(p,f,key);key.report.numbers["minor_rank"]=dimension;key.report.numbers["minor_rows"]=minors.rows;
    key.report.numbers["seconds"]=elapsed(start);key.report.checks["exact_minor_grs_equality"]=true;return key;
}
static Key finish_direct(const Public& p,const Field& f,const Matrix& jets,Matrix& sections) {
    auto start=Clock::now();Pairs pairs=quadric_monomials(p.k);int depth=2*p.D()-1;
    // Only the prefix through 2D-1 enters this construction, even when order 2D is available.
    Matrix C=contacts(f,jets.first_rows(depth+1),pairs,depth),WR=nullspace(f,C);
    Vec target=quadratic_values(f,p.point((p.held+1)%p.n),pairs),R;
    for(int r=0;r<WR.rows;++r)if(Byte value=dot(f,WR.row(r),target)){R=WR.row(r);for(auto& x:R)x=f.div(x,value);break;}
    require(!R.empty(),"direct numerator section is absent");
    int order=2*p.D()-2;Matrix WS=nullspace(f,C.first_rows(order));Vec S;
    for(int r=0;r<WS.rows;++r)if(Byte value=dot(f,WS.row(r),C.row(order))){S=WS.row(r);for(auto& x:S)x=f.div(x,value);break;}
    require(!S.empty(),"direct denominator section is absent");sections=from_rows({R,S},int(pairs.size()));
    Matrix ordinary=nullspace(f,C.first_rows(depth));
    for(int r=0;r<ordinary.rows;++r)require(dot(f,ordinary.row(r),C.row(depth))==0,"unexpected adjacent-order section in a binary Goppa instance");
    std::vector<int> projective(p.n);std::vector<bool> seen(f.q+1);
    for(int j=0;j<p.n;++j) {
        if(j==p.held)projective[j]=0;
        else {
            Vec v=quadratic_values(f,p.point(j),pairs);Byte a=dot(f,R,v),b=dot(f,S,v);
            require(a!=0,"direct numerator vanishes at another support point");projective[j]=b?int(f.root(f.div(a,b))):f.q;
        }
        require(!seen[projective[j]],"duplicate direct support");seen[projective[j]]=true;
    }
    int pole=-1;if(seen[f.q])for(int a=0;a<f.q;++a)if(!seen[a]){pole=a;break;}
    Key key;key.support.resize(p.n);
    for(int j=0;j<p.n;++j)key.support[j]=pole<0?Byte(projective[j]):
        (projective[j]==f.q?Byte(0):f.inv(Byte(projective[j]^pole)));
    complete_key(p,f,key);key.report.numbers["jet_order_used"]=depth;
    key.report.numbers["numerator_space_dimension"]=WR.rows;key.report.numbers["denominator_space_dimension"]=WS.rows;
    key.report.numbers["seconds"]=elapsed(start);key.report.checks["square_root_quotient"]=true;
    key.report.checks["adjacent_order_denominator_absent"]=true;return key;
}
static int support_equivalence(const Field& f,const Vec& original,const Vec& recovered,int conjugates) {
    require(original.size()==recovered.size() && original.size()>=3,"support comparison sizes");
    for(int sigma=0;sigma<conjugates;++sigma) {
        Vec support=original;for(auto& a:support)a=f.power(a,1u<<sigma);
        Matrix A(3,4);for(int j=0;j<3;++j){A(j,0)=support[j];A(j,1)=1;A(j,2)=f.mul(recovered[j],support[j]);A(j,3)=recovered[j];}
        Matrix K=nullspace(f,A);if(K.rows!=1)continue;Vec v=K.row(0);
        if(!(f.mul(v[0],v[3])^f.mul(v[1],v[2])))continue;bool good=true;
        for(int j=0;j<int(support.size());++j) {
            Byte denominator=f.mul(v[2],support[j])^v[3];
            if(!denominator || f.div(f.mul(v[0],support[j])^v[1],denominator)!=recovered[j]){good=false;break;}
        }
        if(good)return sigma;
    }
    throw std::runtime_error("supports are not globally Frobenius/Mobius equivalent");
}
static void save_key(const Public& p,const Key& key,const fs::path& basename) {
    std::ofstream text(basename.string()+".txt");text<<"BGKEY1\n"<<p.m<<' '<<p.t<<' '<<p.n<<' '<<p.modulus<<'\n';
    for(const Vec* v:{&key.support,&key.lambda,&key.g}){for(Byte a:*v)text<<int(a)<<' ';text<<'\n';}
    require(bool(text),"write recovered key failed");
    std::ofstream json(basename.string()+".json");json<<"{\"m\":"<<p.m<<",\"t\":"<<p.t<<",\"n\":"<<p.n<<",\"modulus\":"<<p.modulus<<",\"support\":";
    array_json(json,key.support);json<<",\"grs_multiplier\":";array_json(json,key.lambda);json<<",\"goppa_polynomial\":";array_json(json,key.g);
    json<<",\"verification\":";key.report.write(json);json<<"}\n";require(bool(json),"write key JSON failed");
}
static Key read_key(const Public& p,const fs::path& path) {
    std::ifstream in(path);std::string magic;int m,t,n,mod;
    require(bool(in>>magic>>m>>t>>n>>mod) && magic=="BGKEY1" && m==p.m && t==p.t && n==p.n && mod==p.modulus,"invalid recovered key header");
    Key key;key.support.resize(n);key.lambda.resize(n);key.g.resize(t+1);
    for(Vec* v:{&key.support,&key.lambda,&key.g})for(auto& x:*v){int a;require(bool(in>>a) && a>=0 && a<(1<<m),"invalid key field element");x=Byte(a);}
    std::string extra;require(!(in>>extra),"trailing key data");return key;
}

static Report certify_geometry(const Public& p,const Field& f,const Key& key,const Local& local,
                                const Matrix& jets,const Bootstrap* boot,const fs::path& output) {
    auto start=Clock::now();std::vector<Poly> curve;Matrix Y=p.generator();
    for(int a=0;a<p.k;++a) {
        Vec values(p.n);for(int j=0;j<p.n;++j)values[j]=f.div(Y(a,j),key.lambda[j]);
        Poly poly=interpolate(f,key.support,values);require(int(poly.size())<=p.D()+1,"reconstructed curve degree");curve.push_back(poly);
    }
    std::vector<Poly> P;for(const auto& poly:curve)P.push_back(scale(f,compose(f,poly,Poly{key.support[p.held],1},p.D()+1),key.lambda[p.held]));
    Vec point=p.point(p.held);std::vector<Vec> true_first;
    for(int sigma=0;sigma<p.m;++sigma){Vec v;for(const auto& poly:P)v.push_back(f.power(coefficient(poly,1),1u<<sigma));true_first.push_back(v);}
    for(const auto& v:local.branches) {
        bool found=false;for(const auto& w:true_first)if(same_span(f,from_rows({point,v},p.k),from_rows({point,w},p.k)))found=true;
        require(found,"a recovered branch differs from every reconstructed-curve branch");
    }
    for(const auto& w:true_first) {
        bool found=false;for(const auto& v:local.branches)if(same_span(f,from_rows({point,v},p.k),from_rows({point,w},p.k)))found=true;
        require(found,"a reconstructed-curve branch was not recovered");
    }
    int sigma=-1;
    for(int s=0;s<p.m;++s)if(same_span(f,from_rows({point,jets.row(1)},p.k),from_rows({point,true_first[s]},p.k))){sigma=s;break;}
    require(sigma>=0,"seed branch alignment failed");
    for(auto& poly:P)for(auto& c:poly)c=f.power(c,1u<<sigma);
    FixedSolver solve(f,from_rows({point,true_first[sigma]},p.k).transpose());Poly scalar={1},parameter;
    for(int r=1;r<jets.rows;++r) {
        Vec rhs=jets.row(r);
        for(int a=0;a<p.k;++a)rhs[a]^=coefficient(product(f,scalar,compose(f,P[a],parameter,r+1),r+1),r);
        Vec delta=solve.solve(rhs);scalar.resize(r+1);parameter.resize(r+1);scalar[r]=delta[0];parameter[r]=delta[1];
    }
    require(scalar[0]==1 && parameter[0]==0 && parameter[1]!=0,"invalid local scalar/parameter series");
    for(int a=0;a<p.k;++a) {
        Poly exact=product(f,scalar,compose(f,P[a],parameter,jets.rows),jets.rows);
        for(int r=0;r<jets.rows;++r)require(coefficient(exact,r)==jets(r,a),"coherent jet audit failed");
    }
    Matrix parametrization(p.k,p.D()+1);for(int a=0;a<p.k;++a)for(int r=0;r<=p.D();++r)parametrization(a,r)=coefficient(curve[a],r);
    save_matrix(output/"reconstructed-curve.txt",parametrization);
    scalar.resize(jets.rows);parameter.resize(jets.rows);save_matrix(output/"coherence-series.txt",from_rows({scalar,parameter},jets.rows));
    if(boot) {
        std::vector<Poly> selected=curve;for(auto& poly:selected)for(auto& c:poly)c=f.power(c,1u<<sigma);
        Pairs pairs=quadric_monomials(p.k);Matrix coefficients(2*p.D()+1,int(pairs.size()));
        for(int c=0;c<int(pairs.size());++c) {
            auto [a,b]=pairs[c];Poly poly=product(f,selected[a],selected[b]);
            for(int r=0;r<coefficients.rows;++r)coefficients(r,c)=coefficient(poly,r);
        }
        require(multiply(f,coefficients,boot->forms.transpose()).zero(),"a bootstrapped quadric is not a global curve identity");
        require(coefficients.cols-rank(f,coefficients)==boot->forms.rows,"bootstrapped quadrics do not span the full ideal");
        for(int j=0;j<p.n;++j) {
            Byte a=f.power(key.support[j],1u<<sigma);Vec tangent;
            for(const auto& poly:selected) {
                Poly derivative(poly.size());for(int r=1;r<int(poly.size());r+=2)derivative[r-1]=poly[r];
                tangent.push_back(evaluate(f,derivative,a));
            }
            require(same_span(f,from_rows({p.point(j),tangent},p.k),from_rows({p.point(j),boot->tangents.column(j)},p.k)),"tangent alignment audit failed");
        }
    }
    // All retained Hasse conditions and full jet substitution were checked earlier.
    // Coherence now makes their root count an exact identity certificate on this curve.
    int zeros=p.s*(p.n-1)+jets.rows;
    require(zeros>p.d*p.D(),"insufficient final curve-identity certificate");
    Report report;report.numbers={{"curve_degree_bound",double(p.D())},{"coherent_coefficients",double(jets.rows*p.k)},
        {"frobenius_alignment",double(sigma)},{"certified_root_count",double(zeros)},
        {"holdout_composition_degree_bound",double(p.d*p.D())},{"seconds",elapsed(start)}};
    report.checks={{"coherent_jet",true},{"all_branches_verified",true},{"holdout_curve_identities_certified",true},
        {"uses_original_secret",false}};
    if(boot){report.checks["full_quadratic_ideal_verified"]=true;report.checks["all_tangents_aligned"]=true;}
    return report;
}
static void generate(std::uint64_t seed,const fs::path& output) {
    require(!fs::exists(output),"generation output already exists");Field f(6,67);std::mt19937_64 rng(seed);
    Poly g(7);g[6]=1;
    do { for(int i=0;i<6;++i)g[i]=Byte(rng()&63); } while(!irreducible(f,g));
    Vec support(64);std::iota(support.begin(),support.end(),Byte(0));
    for(int i=63;i>0;--i)std::swap(support[i],support[rng()%unsigned(i+1)]);
    Vec beta(64);for(int j=0;j<64;++j)beta[j]=f.inv(evaluate(f,g,support[j]));
    auto check=expand_binary(f,grs(f,support,6,beta));
    auto rows=binary_rref(binary_kernel(check,64),64);require(rows.size()==28,"generated code dimension differs from 28");
    Vec squared=beta;for(auto& x:squared)x=f.sq(x);
    require(expand_binary(f,grs(f,support,12,squared))==check,"binary square-free Goppa identity failed");
    Public p{6,6,64,28,67,0,5,4,{}};p.points.resize(64);
    for(int a=0;a<28;++a)for(int j=0;j<64;++j)p.points[j]|=((rows[a]>>j)&1)<<a;
    while(p.held<64 && __builtin_popcountll(p.points[p.held])==1)++p.held;
    require(p.held<64,"no held-out column");fs::create_directories(output/"audit");write_public(output/"public.txt",p);
    std::ofstream out(output/"audit"/"secret.txt");out<<"BGSK1\n6 6 64 67\n";
    for(Byte x:support)out<<int(x)<<' ';out<<'\n';for(Byte x:g)out<<int(x)<<' ';out<<'\n';
    std::ofstream info(output/"generation.json");info<<"{\"seed\":"<<seed<<",\"generator\":\"C++ mt19937_64\",\"field_modulus\":67,\"binary_rank\":36,\"k\":28}\n";
    std::cout<<"Generated public [64,28] key; secret audit data are in a separate directory.\n";
}
static void audit_secret(const Public& p,const Field& f,const fs::path& path,const Key& key) {
    std::ifstream in(path);std::string magic;int m,t,n,mod;
    require(bool(in>>magic>>m>>t>>n>>mod) && magic=="BGSK1" && m==p.m && t==p.t && n==p.n && mod==p.modulus,"invalid secret audit header");
    Vec support(n),g(t+1);for(Vec* v:{&support,&g})for(auto& x:*v){int a;require(bool(in>>a) && a>=0 && a<f.q,"secret audit encoding");x=Byte(a);}
    require(irreducible(f,g),"original polynomial is not irreducible");
    Vec beta(n);for(int j=0;j<n;++j)beta[j]=f.inv(evaluate(f,g,support[j]));
    require(expand_binary(f,grs(f,support,t,beta))==binary_rref(binary_kernel(p.rows(),p.n),p.n),"original audit key does not describe public input");
    verify_key(p,f,key);int sigma=support_equivalence(f,support,key.support,p.m);
    std::cout<<"{\"status\":\"passed\",\"global_support_equivalence\":true,\"frobenius_power\":"<<sigma<<"}\n";
}
static void recover(const Public& p,const fs::path& kernel_path,const fs::path& output,const std::string& route) {
    require(route=="both" || route=="minor" || route=="direct","route must be both, minor, or direct");
    require(!fs::exists(output),"recovery output already exists; use a fresh directory");fs::create_directories(output);
    auto started=Clock::now(),tick=started;Field f(p.m,p.modulus);
    std::cout<<"Read and verify the public holdout kernel\n"<<std::flush;
    Binary K=read_kernel(p,kernel_path);verify_hasse(p,K.get());double kernel_seconds=elapsed(tick);
    std::cout<<"Separate the six branches and check every continuation rank\n"<<std::flush;
    tick=Clock::now();Local local=recover_branches(p,f,K.get());local.report.numbers["seconds"]=elapsed(tick);
    save_matrix(output/"branches.txt",from_rows(local.branches,p.k));
    int depth=route=="minor"?2*((std::max(0,2*p.D()-p.n+1)+1)/2):2*p.D();
    std::cout<<"Recover a binary-normalized jet through order "<<depth<<"\n"<<std::flush;
    Jet jet=continue_jet(p,f,K.get(),local,depth);save_matrix(output/"jet.txt",jet.coefficients);
    Bootstrap boot;Key minor,direct;Matrix minors,sections;
    if(route!="direct") {
        std::cout<<"Bootstrap one branch, align tangents, and recover the minor-code key\n"<<std::flush;
        boot=bootstrap(p,f,jet.coefficients);save_matrix(output/"quadrics.txt",boot.forms);save_matrix(output/"tangents.txt",boot.tangents);
        minor=finish_minor(p,f,boot,minors);save_matrix(output/"minor-matrix.txt",minors);save_key(p,minor,output/"minor-key");
    }
    if(route!="minor") {
        std::cout<<"Recover the direct square-root quotient key\n"<<std::flush;
        direct=finish_direct(p,f,jet.coefficients,sections);save_matrix(output/"direct-sections.txt",sections);save_key(p,direct,output/"direct-key");
    }
    if(route=="both")require(support_equivalence(f,minor.support,direct.support,1)==0,"the two routes do not share a Mobius frame");
    const Key& recovered=route=="minor"?minor:direct;
    std::cout<<"Certify jet coherence and global curve identities from the recovered public key\n"<<std::flush;
    Report geometry=certify_geometry(p,f,recovered,local,jet.coefficients,route=="direct"?nullptr:&boot,output);
    Key bad=recovered;bad.lambda[0]=f.mul(bad.lambda[0],2);bool rejected=false;
    try { verify_key(p,f,bad); } catch(const std::exception&) { rejected=true; }
    require(rejected,"negative verification control accepted a corrupted multiplier");
    std::ofstream out(output/"recovery.json");out<<"{\"status\":\"passed\",\"route\":\""<<route<<"\",\"m\":"<<p.m
        <<",\"t\":"<<p.t<<",\"n\":"<<p.n<<",\"k\":"<<p.k<<",\"D\":"<<p.D()<<",\"heldout\":"<<p.held
        <<",\"modulus\":"<<p.modulus<<",\"kernel_dimension\":"<<K->nrows<<",\"monomials\":"<<K->ncols
        <<",\"kernel_input_verification_seconds\":"<<kernel_seconds<<",\"initial_forced_zeros\":"<<p.s*(p.n-1)
        <<",\"initial_composition_degree_bound\":"<<p.d*p.D()<<",\"initial_ordinary_zero_count_suffices\":false"
        <<",\"public_inputs_only\":true,\"full_hasse_verified\":true,\"kernel_independence_verified\":true,\"branches\":";
    local.report.write(out);out<<",\"jet\":";jet.report.write(out);out<<",\"geometry\":";geometry.write(out);
    if(route!="direct"){out<<",\"bootstrap\":";boot.report.write(out);out<<",\"minor\":";minor.report.write(out);}
    if(route!="minor"){out<<",\"direct\":";direct.report.write(out);}
    out<<",\"tampered_multiplier_rejected\":true";
    if(route=="both")out<<",\"routes_mobius_equivalent\":true";
    out<<",\"total_seconds_excluding_kernel_acquisition\":"<<elapsed(started)<<"}\n";
    require(bool(out),"report write failed");std::cout<<"PASS: recovered all 64 support positions and an equivalent degree-six Goppa key.\n";
}

static void self_test() {
    Field f(6,67);
    for(int a=0;a<64;++a) {
        require(f.sq(f.root(Byte(a)))==a,"square root self-test");
        for(int b=0;b<64;++b)for(int c=0;c<64;++c)
            require(f.mul(Byte(a),Byte(b^c))==(f.mul(Byte(a),Byte(b))^f.mul(Byte(a),Byte(c))),"field distributivity");
    }
    std::mt19937_64 rng(17);Matrix A(7,11);for(auto& x:A.data)x=Byte(rng()&63);
    Matrix N=nullspace(f,A);require(multiply(f,A,N.transpose()).zero() && N.rows==A.cols-rank(f,A),"nullspace self-test");
    Vec x(11);for(auto& c:x)c=Byte(rng()&63);FixedSolver solve(f,A);
    require(multiply(f,A,solve.solve(multiply(f,A,x)))==multiply(f,A,x),"linear-solve self-test");
    Vec support={0,1,2,3,4,5,6,7},values;Poly poly={7,4,21,3};
    for(Byte a:support)values.push_back(evaluate(f,poly,a));require(interpolate(f,support,values)==poly,"interpolation self-test");
    require(irreducible(f,Poly{19,42,61,23,18,33,1}),"known irreducible polynomial test");
    require(!irreducible(f,product(f,Poly{1,1},Poly{2,1})),"reducible polynomial test");
    constexpr int k=7,d=4,depth=12;Matrix jets(depth+1,k);for(auto& a:jets.data)a=Byte(rng()&63);
    auto mons=masks(k,d);Matrix arcs=jets.transpose(),whole(int(mons.size()),depth+1);
    require(series_products(k,d,6,67,depth,arcs.data.data(),whole.data.data())==0,"native products self-test");
    std::unique_ptr<void,decltype(&series_cache_free)> cache(series_cache_create(k,d,6,67,depth),series_cache_free);
    Vec output(mons.size());
    for(int r=0;r<=depth;++r) {
        Vec jet=jets.row(r);require(series_cache_set(cache.get(),r,jet.data())==0,"cache set self-test");
        require(series_cache_eval(cache.get(),r,output.data())==0,"cache eval self-test");
        for(int c=0;c<int(mons.size());++c)require(output[c]==whole(c,r),"incremental product disagrees with fresh products");
    }
    for(int c=0;c<int(mons.size());++c) {
        Poly expected={1};for(int a=0;a<k;++a)if((mons[c]>>a)&1)expected=product(f,expected,arcs.row(a),depth+1);
        for(int r=0;r<=depth;++r)require(coefficient(expected,r)==whole(c,r),"native series product disagrees with independent multiplication");
    }
    // Independent literal-Hasse check of the native local-map generator.
    auto ms=masks(7,3);int width=(int(ms.size())+7)/8;Vec maps_data((1+7+21)*width);
    require(local_rows(7,3,85,maps_data.data())==0,"local map test");int row=0;
    for(int order=0;order<=2;++order)for(auto deriv:masks(7,order)) {
        for(int c=0;c<int(ms.size());++c) {
            bool expected=(deriv&ms[c])==deriv && !((ms[c]^deriv)&~std::uint64_t(85));
            require(bool((maps_data[row*width+c/8]>>(c%8))&1)==expected,"local Taylor map mismatch");
        }
        ++row;
    }
    std::cout<<"PASS: field, matrices, interpolation, irreducibility, series cache, fresh products, and Hasse maps.\n";
}
int main(int argc,char** argv) {
    try {
        require(argc>=2,"commands: self-test | generate SEED DIR | prepare PUBLIC OPERATOR | recover PUBLIC KERNEL DIR [both|minor|direct] | verify PUBLIC KEY | audit PUBLIC SECRET KEY");
        std::string command=argv[1];
        if(command=="self-test"){require(argc==2,"self-test arguments");self_test();}
        else if(command=="generate"){require(argc==4,"generate SEED DIR");generate(std::stoull(argv[2]),argv[3]);}
        else if(command=="prepare"){require(argc==4,"prepare PUBLIC OPERATOR");prepare(read_public(argv[2]),argv[3]);}
        else if(command=="recover"){require(argc==5 || argc==6,"recover PUBLIC KERNEL DIR [both|minor|direct]");recover(read_public(argv[2]),argv[3],argv[4],argc==6?argv[5]:"both");}
        else if(command=="verify"){require(argc==4,"verify PUBLIC KEY");Public p=read_public(argv[2]);Field f(p.m,p.modulus);auto report=verify_key(p,f,read_key(p,argv[3]));report.write(std::cout);std::cout<<'\n';}
        else if(command=="audit"){require(argc==5,"audit PUBLIC SECRET KEY");Public p=read_public(argv[2]);Field f(p.m,p.modulus);audit_secret(p,f,argv[3],read_key(p,argv[4]));}
        else throw std::runtime_error("unknown command");
        return 0;
    } catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
