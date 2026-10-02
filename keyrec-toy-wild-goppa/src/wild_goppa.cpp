#include "public.hpp"
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

using namespace toy;using namespace wild;
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
static void save_matrix(const fs::path& path,const Matrix& a) {
    std::ofstream out(path);out<<"GF64M1\n"<<a.rows<<' '<<a.cols<<'\n';
    for(int r=0;r<a.rows;++r){for(int c=0;c<a.cols;++c)out<<int(a(r,c))<<' ';out<<'\n';}
    require(bool(out),"write matrix failed");
}

struct Kernel {
    int rows,cols,full_dimension;
    Matrix digits;
    Binary low,high,vectors;
    Kernel(int r,int c,int full):rows(r),cols(c),full_dimension(full),digits(r,c),
        low(binary(r,c)),high(binary(r,c)),vectors(binary(2*c,r)) {}
};
static Kernel read_kernel(const Public& p,const fs::path& path) {
    std::ifstream in(path,std::ios::binary);char magic[4];
    require(bool(in.read(magic,4)) && std::string(magic,4)=="WGK1","invalid kernel magic");
    int k=int(read_u32(in)),n=int(read_u32(in)),D=int(read_u32(in)),held=int(read_u32(in)),d=int(read_u32(in));
    int N=int(read_u32(in)),M=int(read_u32(in)),full=int(read_u32(in));
    require(k==p.k && n==p.n && D==p.D() && held==p.held && d==p.d && N==int(masks(p.k,p.d).size()),"kernel/public mismatch");
    require(M>0 && M<=4096 && M<=full && full<=N,"kernel dimensions");Kernel K(M,N,full);Vec packed((N+3)/4);
    auto expanded=binary(2*M,2*N);
    for(int r=0;r<M;++r) {
        require(bool(in.read(reinterpret_cast<char*>(packed.data()),packed.size())),"truncated kernel");
        for(int c=0;c<N;++c) {
            Byte a=(packed[c/4]>>(2*(c%4)))&3;K.digits(r,c)=a;
            for(int bit=0;bit<2;++bit)if((a>>bit)&1){mzd_write_bit(K.vectors.get(),2*c+bit,r,1);mzd_write_bit(expanded.get(),2*r,2*c+bit,1);}
            Byte b=mul4[a][2];for(int bit=0;bit<2;++bit)if((b>>bit)&1)mzd_write_bit(expanded.get(),2*r+1,2*c+bit,1);
            if(a&1)mzd_write_bit(K.low.get(),r,c,1);if(a&2)mzd_write_bit(K.high.get(),r,c,1);
        }
    }
    char extra;require(!in.get(extra),"trailing kernel bytes");
    require(mzd_echelonize_m4ri(expanded.get(),0,0)==2*M,"exported polynomials are dependent over GF(4)");return K;
}
static Matrix apply_binary(const Field& f,const mzd_t* K,const Matrix& values) {
    require(K->ncols==values.rows && values.cols>0,"binary product dimensions");auto B=binary(values.rows,values.cols*f.m);
    for(int r=0;r<values.rows;++r)for(int c=0;c<values.cols;++c)for(int bit=0;bit<f.m;++bit)
        if((values(r,c)>>bit)&1)mzd_write_bit(B.get(),r,c*f.m+bit,1);
    Binary C(mzd_mul(nullptr,K,B.get(),0));Matrix result(K->nrows,values.cols);
    for(int r=0;r<result.rows;++r)for(int c=0;c<result.cols;++c)for(int bit=0;bit<f.m;++bit)
        result(r,c)|=mzd_read_bit(C.get(),r,c*f.m+bit)<<bit;
    return result;
}
static Matrix apply(const Field& f,const Kernel& K,const Matrix& values) {
    Matrix low=apply_binary(f,K.low.get(),values),high=apply_binary(f,K.high.get(),values);
    for(std::size_t i=0;i<low.data.size();++i)low.data[i]^=f.mul(embedding[2],high.data[i]);return low;
}
static void write_block(mzd_t* A,int row,int col,Byte a) {
    Byte b=mul4[a][2];
    if(a&1)mzd_write_bit(A,2*row,2*col,1);if(a&2)mzd_write_bit(A,2*row+1,2*col,1);
    if(b&1)mzd_write_bit(A,2*row,2*col+1,1);if(b&2)mzd_write_bit(A,2*row+1,2*col+1,1);
}
static void verify_hasse(const Public& p,const Kernel& K) {
    auto mons=masks(p.k,p.d);std::unordered_map<std::uint64_t,int> index;
    for(int c=0;c<int(mons.size());++c)index[mons[c]]=c;
    std::vector<std::uint64_t> derivatives;for(int r=0;r<p.s;++r){auto v=masks(p.k,r);derivatives.insert(derivatives.end(),v.begin(),v.end());}
    for(int j:p.retained_nonunits()) {
        std::uint64_t active=0;for(int a=0;a<p.k;++a)if(p.digit(a,j))active|=std::uint64_t(1)<<a;
        auto A=binary(2*int(derivatives.size()),2*K.cols);
        // Assemble by derivative index and extensions, independently of the
        // monomial-submask enumeration used by the holdout program.
        for(int r=0;r<int(derivatives.size());++r) {
            auto derivative=derivatives[r];int degree=p.d-__builtin_popcountll(derivative);
            for(auto extra:masks(p.k,degree,active&~derivative))write_block(A.get(),r,index.at(derivative|extra),monomial_value(p,j,extra));
        }
        Binary residual(mzd_mul(nullptr,A.get(),K.vectors.get(),0));
        require(mzd_is_zero(residual.get()),"kernel violates a retained Hasse condition");
    }
}
struct Local { Matrix J,Q,B,stationary;std::vector<std::pair<int,int>> pairs;std::vector<Vec> branches;Report report; };
static Matrix polar(const Field& f,const Local& local,const Vec& v) {
    Matrix P(local.Q.cols,int(v.size()));for(int r=0;r<int(local.pairs.size());++r){auto [a,b]=local.pairs[r];P(r,a)=v[b];P(r,b)=v[a];}
    return multiply(f,local.Q,P);
}
static Local recover_branches(const Public& p,const Field& f,const Kernel& K) {
    Local local;for(int a=0;a<p.k;++a)for(int b=a+1;b<p.k;++b)local.pairs.emplace_back(a,b);
    std::vector<std::uint64_t> indices={0};for(int a=0;a<p.k;++a)indices.push_back(std::uint64_t(1)<<a);
    for(auto [a,b]:local.pairs)indices.push_back((std::uint64_t(1)<<a)|(std::uint64_t(1)<<b));
    std::unordered_map<std::uint64_t,int> slots;for(int c=0;c<int(indices.size());++c)slots[indices[c]]=c;
    auto mons=masks(p.k,p.d);Matrix functionals(K.cols,int(indices.size()));
    for(int c=0;c<K.cols;++c)for(auto derivative=mons[c];;derivative=(derivative-1)&mons[c]) {
        auto found=slots.find(derivative);if(found!=slots.end())functionals(c,found->second)=embedding[monomial_value(p,p.held,mons[c]^derivative)];
        if(!derivative)break;
    }
    Matrix maps=apply(f,K,functionals);local.J=Matrix(K.rows,p.k);local.Q=Matrix(K.rows,int(local.pairs.size()));
    for(int r=0;r<K.rows;++r) {
        require(!maps(r,0),"a holdout polynomial is nonzero at the seed");
        for(int c=0;c<p.k;++c)local.J(r,c)=maps(r,c+1);
        for(int c=0;c<local.Q.cols;++c)local.Q(r,c)=maps(r,c+1+p.k);
    }
    Matrix W=nullspace(f,local.J);require(W.rows==p.m+1,"gradient kernel does not have dimension m+1");
    std::vector<Vec> basis={p.point(p.held)};
    for(int r=0;r<W.rows;++r){auto trial=basis;trial.push_back(W.row(r));if(rank(f,from_rows(trial,p.k))>int(basis.size()))basis.push_back(W.row(r));}
    require(int(basis.size())==p.m+1,"gradient quotient basis");basis.erase(basis.begin());local.B=from_rows(basis,p.k);
    local.stationary=nullspace(f,local.J.transpose());
    std::vector<std::pair<int,int>> qm;for(int a=0;a<p.m;++a)for(int b=a;b<p.m;++b)qm.emplace_back(a,b);
    Matrix substitution(local.Q.cols,int(qm.size()));
    for(int r=0;r<substitution.rows;++r)for(int c=0;c<substitution.cols;++c) {
        auto [a,b]=local.pairs[r];auto [u,v]=qm[c];substitution(r,c)=f.mul(local.B(u,a),local.B(v,b));
        if(u!=v)substitution(r,c)^=f.mul(local.B(v,a),local.B(u,b));
    }
    Matrix stationary_quadrics=multiply(f,local.stationary,local.Q);
    Matrix restrictions=multiply(f,stationary_quadrics,substitution);
    require(rank(f,restrictions)==p.m*(p.m-1)/2,"stationary quadratic rank criterion fails");
    Matrix quotient=nullspace(f,restrictions);require(quotient.rows==p.m,"quadratic quotient dimension");
    std::vector<Matrix> multiplication;
    for(int a=0;a<p.m;++a) {
        Matrix M(p.m,p.m);for(int j=0;j<p.m;++j) {
            auto pair=std::make_pair(std::min(a,j),std::max(a,j));int c=int(std::find(qm.begin(),qm.end(),pair)-qm.begin());
            for(int r=0;r<p.m;++r)M(r,j)=quotient(r,c);
        }multiplication.push_back(M);
    }
    Matrix A1=multiplication[0];require(rank(f,A1)==p.m,"first projective chart misses a branch");
    for(int mask=2;mask<(1<<p.m) && local.branches.empty();++mask) {
        Matrix A2(p.m,p.m);for(int a=0;a<p.m;++a)if((mask>>a)&1)
            for(std::size_t c=0;c<A2.data.size();++c)A2.data[c]^=multiplication[a].data[c];
        std::vector<Vec> found;
        for(int theta=0;theta<f.q;++theta) {
            Matrix pencil=A2;for(std::size_t c=0;c<pencil.data.size();++c)pencil.data[c]^=f.mul(Byte(theta),A1.data[c]);
            Matrix left=nullspace(f,pencil.transpose());if(left.rows!=1)continue;
            found.push_back(multiply(f,local.B.transpose(),multiply(f,A1.transpose(),left.row(0))));
        }
        if(int(found.size())==p.m)local.branches=found;
    }
    require(int(local.branches.size())==p.m,"could not recover all three branches");
    for(const Vec& v:local.branches) {
        Vec residual=multiply(f,local.J,v);require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"branch gradient residual");
        Vec products;for(auto [a,b]:local.pairs)products.push_back(f.mul(v[a],v[b]));
        residual=multiply(f,stationary_quadrics,products);require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"stationary quadratic residual");
        Matrix A=stack(local.J,multiply(f,local.stationary,polar(f,local,v)));
        require(rank(f,A)==p.k-2,"stationary continuation does not have the two-dimensional gauge kernel");
    }
    local.report.numbers={{"gradient_rank",double(p.k-p.m-1)},{"gradient_nullity",double(p.m+1)},
        {"stationary_quadratic_rank",3},{"branches",3},{"continuation_rank_all_branches",double(p.k-2)}};
    return local;
}
struct Jet {Matrix coefficients;Report report;};
static Jet continue_jet(const Public& p,const Field& f,const Kernel& K,const Local& local,int depth) {
    require(depth>=2,"jet depth must be at least two");auto started=Clock::now();
    Matrix A=stack(local.J,multiply(f,local.stationary,polar(f,local,local.branches[0])));FixedSolver solver(f,A);
    Matrix jets(depth+1,p.k);jets.set_row(0,p.point(p.held));jets.set_row(1,local.branches[0]);
    std::unique_ptr<void,decltype(&series_cache_free)> cache(series_cache_create(p.k,p.d,6,p.modulus,depth+1),series_cache_free);
    require(bool(cache),"series cache allocation");Matrix values(K.cols,1);
    auto coefficient=[&](int r){require(series_cache_eval(cache.get(),r,values.data.data())==0,"series coefficient evaluation");return apply(f,K,values).column(0);};
    for(int r=0;r<2;++r){Vec v=jets.row(r);require(series_cache_set(cache.get(),r,v.data())==0,"set initial jet");auto residual=coefficient(r);require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"initial jet residual");}
    int first_full=-1;
    for(int r=2;r<=depth;++r) {
        Vec rhs=coefficient(r),next=multiply(f,local.stationary,coefficient(r+1));rhs.insert(rhs.end(),next.begin(),next.end());
        Vec v=solver.solve(rhs);jets.set_row(r,v);require(series_cache_set(cache.get(),r,v.data())==0,"set recovered jet");
        Vec residual=coefficient(r);require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"regular recurrence residual");
        residual=multiply(f,local.stationary,coefficient(r+1));require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"stationary next-order residual");
        if(first_full<0 && rank(f,jets.first_rows(r+1))==p.k)first_full=r;
    }
    double recurrence=elapsed(started);auto tick=Clock::now();Matrix arcs=jets.transpose(),fresh(K.cols,depth+1);
    require(series_products(p.k,p.d,6,p.modulus,depth,arcs.data.data(),fresh.data.data())==0,"fresh formal products");
    require(apply(f,K,fresh).zero(),"fresh full substitution failed");
    Report report;report.numbers={{"depth",double(depth)},{"continuation_rank",double(solver.columns.size())},
        {"first_full_coefficient_span_order",double(first_full)},{"recurrence_seconds",recurrence},{"fresh_substitution_seconds",elapsed(tick)}};
    report.checks={{"general_stationary_continuation",true},{"full_substitution",true},{"all_exported_polynomials_used",true}};
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

static const std::vector<Vec>& components(const Field& f) {
    static std::vector<Vec> table;
    if(table.empty()) {
        table.assign(64,Vec(3));std::vector<bool> seen(64);
        for(int a=0;a<4;++a)for(int b=0;b<4;++b)for(int c=0;c<4;++c) {
            Byte x=embedding[a]^f.mul(embedding[b],2)^f.mul(embedding[c],4);
            require(!seen[x],"invalid GF(64)/GF(4) basis");seen[x]=true;
            table[x]={embedding[a],embedding[b],embedding[c]};
        }
    }
    return table;
}
static Matrix expand_subfield(const Field& f,const Matrix& H) {
    const auto& basis=components(f);Matrix out(3*H.rows,H.cols);
    for(int r=0;r<H.rows;++r)for(int j=0;j<H.cols;++j)for(int a=0;a<3;++a)out(3*r+a,j)=basis[H(r,j)][a];
    return out;
}
static std::vector<Poly> curve_from_key(const Public& p,const Field& f,const Key& key) {
    Matrix Y=p.generator();std::vector<Poly> curve;
    for(int a=0;a<p.k;++a) {
        Vec values(p.n);for(int j=0;j<p.n;++j)values[j]=f.div(Y(a,j),key.lambda[j]);
        Poly poly=interpolate(f,key.support,values);require(int(poly.size())<=p.D()+1,"reconstructed curve degree");curve.push_back(poly);
    }
    return curve;
}
static bool same_projective_point(const Field& f,const Vec& a,const Vec& b) {
    int pivot=-1;for(int j=0;j<int(b.size());++j)if(b[j]){pivot=j;break;}
    if(pivot<0 || !a[pivot])return false;Byte scale=f.div(a[pivot],b[pivot]);
    for(int j=0;j<int(b.size());++j)if(a[j]!=f.mul(scale,b[j]))return false;return true;
}
static Report verify_key(const Public& p,const Field& f,const Key& key) {
    require(int(key.support.size())==p.n && int(key.lambda.size())==p.n,"key lengths");
    for(Byte a:key.support)require(a<64,"support encoding");for(Byte a:key.lambda)require(a>0 && a<64,"multiplier encoding");
    Matrix Y=p.generator(),H=expand_subfield(f,parity(f,key.support,p.D()+1,key.lambda));
    require(rank(f,H)==p.n-p.k && multiply(f,Y,H.transpose()).zero(),"recovered alternant code differs from public GF(4) code");
    require(int(key.g.size())==p.t+1 && key.g.back()==1 && irreducible(f,key.g),"recovered radical is not monic irreducible of the prescribed degree");
    Vec values(p.n);for(int j=0;j<p.n;++j){values[j]=evaluate(f,key.g,key.support[j]);require(values[j]!=0,"radical vanishes on support");}
    for(int exponent:{3,4}) {
        Vec beta(p.n);for(int j=0;j<p.n;++j)beta[j]=f.inv(f.power(values[j],exponent));
        Matrix check=expand_subfield(f,grs(f,key.support,exponent*p.t,beta));
        require(rank(f,check)==p.n-p.k && multiply(f,Y,check.transpose()).zero(),"recovered wild Goppa code differs from the public code");
    }
    Vec derivative=support_derivatives(f,key.support);Byte common=0;
    for(int j=0;j<p.n;++j) {
        Byte value=f.div(f.mul(key.lambda[j],derivative[j]),f.power(values[j],4));
        if(j==0)common=value;require(value==common,"wild Goppa multiplier identity");
    }
    Report report;report.numbers={{"support_count",double(p.n)},{"subfield_parity_rank",double(p.n-p.k)},{"radical_degree",double(p.t)}};
    report.checks={{"distinct_finite_support",true},{"exact_public_alternant_equality",true},
        {"exact_public_goppa_cube_equality",true},{"exact_public_goppa_fourth_power_equality",true},
        {"radical_irreducible",true},{"multiplier_fourth_power_identity",true}};return report;
}
static void complete_key(const Public& p,const Field& f,Key& key) {
    key.lambda=recover_multiplier(f,p.generator(),key.support,p.D()+1);
    auto curve=curve_from_key(p,f,key);Vec infinity(p.k);
    for(int j=0;j<p.n;++j){Vec point=p.point(j);for(int a=0;a<p.k;++a)infinity[a]^=point[a];}
    std::vector<int> matches;
    for(int x=0;x<=64;++x) {
        Vec point;for(const auto& poly:curve)point.push_back(x==64?coefficient(poly,p.D()):evaluate(f,poly,Byte(x)));
        if(same_projective_point(f,point,infinity))matches.push_back(x);
    }
    require(matches.size()==1,"public infinity does not identify a unique parameter");
    int pole=matches[0];
    if(pole!=64) {
        for(auto& a:key.support)a=f.inv(a^Byte(pole));
        key.lambda=recover_multiplier(f,p.generator(),key.support,p.D()+1);
    }
    Vec dp=support_derivatives(f,key.support),values(p.n);
    for(int j=0;j<p.n;++j)values[j]=f.power(f.mul(key.lambda[j],dp[j]),16);
    key.g=interpolate(f,key.support,values);require(int(key.g.size())==p.t+1,"radical interpolation degree");
    key.g=scale(f,key.g,f.inv(key.g.back()));Report verified=verify_key(p,f,key);
    key.report.numbers.insert(verified.numbers.begin(),verified.numbers.end());key.report.checks.insert(verified.checks.begin(),verified.checks.end());
    key.report.checks["public_infinity_identified"]=true;
}
static Key finish_minor(const Public& p,const Field& f,const Bootstrap& boot,Matrix& minors) {
    auto start=Clock::now();minors=Matrix(p.k*(p.k-1)/2,p.n);int row=0;Matrix Y=p.generator();
    for(int a=0;a<p.k;++a)for(int b=a+1;b<p.k;++b) {
        for(int j=0;j<p.n;++j)minors(row,j)=f.root(f.mul(Y(a,j),boot.tangents(b,j))^f.mul(Y(b,j),boot.tangents(a,j)));++row;
    }
    int dimension=rank(f,minors);require(dimension==p.D(),"minor code does not attain its GRS dimension bound");
    Key key;key.support=recover_grs_support(f,minors);Vec lambda=recover_multiplier(f,minors,key.support,dimension);
    require(same_span(f,minors,grs(f,key.support,dimension,lambda)),"minor-GRS row-space equality failed");
    complete_key(p,f,key);key.report.numbers["minor_rank"]=dimension;key.report.numbers["minor_rows"]=minors.rows;
    key.report.numbers["seconds"]=elapsed(start);key.report.checks["exact_minor_grs_equality"]=true;return key;
}
static Key finish_direct(const Public& p,const Field& f,const Matrix& jets,Matrix& sections) {
    auto start=Clock::now();Pairs pairs=quadric_monomials(p.k);int L=2*p.D(),depth=L-1;
    Matrix C=contacts(f,jets.first_rows(depth+1),pairs,depth),WR=nullspace(f,C);
    Vec target=quadratic_values(f,p.point((p.held+1)%p.n),pairs),R;
    for(int r=0;r<WR.rows;++r)if(Byte a=dot(f,WR.row(r),target)){R=WR.row(r);for(auto& x:R)x=f.div(x,a);break;}
    require(!R.empty(),"direct numerator section is absent");
    Matrix WS=nullspace(f,C.first_rows(L-1));Vec S;
    for(int r=0;r<WS.rows;++r)if(Byte a=dot(f,WS.row(r),C.row(L-1))){S=WS.row(r);for(auto& x:S)x=f.div(x,a);break;}
    require(!S.empty(),"adjacent-order denominator section is absent");sections=from_rows({R,S},int(pairs.size()));
    std::vector<int> projective(p.n);std::vector<bool> seen(65);
    for(int j=0;j<p.n;++j) {
        if(j==p.held)projective[j]=0;
        else {Vec value=quadratic_values(f,p.point(j),pairs);Byte a=dot(f,R,value),b=dot(f,S,value);require(a!=0,"numerator vanishes at another public point");projective[j]=b?int(f.div(a,b)):64;}
        require(!seen[projective[j]],"duplicate direct support");seen[projective[j]]=true;
    }
    int pole=-1;if(seen[64])for(int a=0;a<64;++a)if(!seen[a]){pole=a;break;}
    Key key;key.support.resize(p.n);for(int j=0;j<p.n;++j)key.support[j]=pole<0?Byte(projective[j]):(projective[j]==64?Byte(0):f.inv(Byte(projective[j]^pole)));
    complete_key(p,f,key);key.report.numbers["jet_order_used"]=depth;key.report.numbers["numerator_contact_order"]=L;
    key.report.numbers["denominator_contact_order"]=L-1;key.report.numbers["seconds"]=elapsed(start);
    key.report.checks["ordinary_evaluation_quotient"]=true;return key;
}
static int support_equivalence(const Field& f,const Vec& original,const Vec& recovered,int conjugates) {
    require(original.size()==recovered.size() && original.size()>=3,"support comparison dimensions");
    for(int sigma=0;sigma<conjugates;++sigma) {
        Vec a=original;for(auto& x:a)x=f.power(x,1u<<(2*sigma));Matrix A(3,4);
        for(int j=0;j<3;++j){A(j,0)=a[j];A(j,1)=1;A(j,2)=f.mul(recovered[j],a[j]);A(j,3)=recovered[j];}
        Matrix K=nullspace(f,A);if(K.rows!=1)continue;Vec v=K.row(0);
        if(!(f.mul(v[0],v[3])^f.mul(v[1],v[2])))continue;bool good=true;
        for(int j=0;j<int(a.size());++j){Byte denominator=f.mul(v[2],a[j])^v[3];if(!denominator || f.div(f.mul(v[0],a[j])^v[1],denominator)!=recovered[j]){good=false;break;}}
        if(good)return sigma;
    }
    throw std::runtime_error("supports are not globally Frobenius/Mobius equivalent");
}
static void save_key(const Public& p,const Key& key,const fs::path& basename) {
    std::ofstream text(basename.string()+".txt");text<<"WGKEY1\n"<<p.m<<' '<<p.t<<' '<<p.n<<' '<<p.modulus<<'\n';
    for(const Vec* v:{&key.support,&key.lambda,&key.g}){for(Byte x:*v)text<<int(x)<<' ';text<<'\n';}require(bool(text),"write key failed");
    std::ofstream json(basename.string()+".json");json<<"{\"q\":4,\"m\":"<<p.m<<",\"t\":"<<p.t<<",\"n\":"<<p.n<<",\"k\":"<<p.k<<",\"modulus\":"<<p.modulus<<",\"support\":";
    array_json(json,key.support);json<<",\"grs_multiplier\":";array_json(json,key.lambda);json<<",\"goppa_radical\":";array_json(json,key.g);
    json<<",\"goppa_exponent\":3,\"verification\":";key.report.write(json);json<<"}\n";require(bool(json),"write key JSON failed");
}
static Key read_key(const Public& p,const fs::path& path) {
    std::ifstream in(path);std::string magic;int m,t,n,mod;
    require(bool(in>>magic>>m>>t>>n>>mod) && magic=="WGKEY1" && m==p.m && t==p.t && n==p.n && mod==p.modulus,"invalid key header");
    Key key;key.support.resize(n);key.lambda.resize(n);key.g.resize(t+1);
    for(Vec* v:{&key.support,&key.lambda,&key.g})for(auto& x:*v){int a;require(bool(in>>a) && a>=0 && a<64,"invalid field encoding");x=Byte(a);}
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
    for(int sigma=0;sigma<p.m;++sigma){Vec v;for(const auto& poly:P)v.push_back(f.power(coefficient(poly,1),1u<<(2*sigma)));true_first.push_back(v);}
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
    for(auto& poly:P)for(auto& c:poly)c=f.power(c,1u<<(2*sigma));
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
        std::vector<Poly> selected=curve;for(auto& poly:selected)for(auto& c:poly)c=f.power(c,1u<<(2*sigma));
        Pairs pairs=quadric_monomials(p.k);Matrix coefficients(2*p.D()+1,int(pairs.size()));
        for(int c=0;c<int(pairs.size());++c) {
            auto [a,b]=pairs[c];Poly poly=product(f,selected[a],selected[b]);
            for(int r=0;r<coefficients.rows;++r)coefficients(r,c)=coefficient(poly,r);
        }
        require(multiply(f,coefficients,boot->forms.transpose()).zero(),"a bootstrapped quadric is not a global curve identity");
        require(coefficients.cols-rank(f,coefficients)==boot->forms.rows,"bootstrapped quadrics do not span the full ideal");
        for(int j=0;j<p.n;++j) {
            Byte a=f.power(key.support[j],1u<<(2*sigma));Vec tangent;
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

static void write_public(const Public& p,const fs::path& path) {
    std::ofstream out(path);out<<"WGPK1\n"<<p.m<<' '<<p.t<<' '<<p.n<<' '<<p.k<<' '<<p.modulus<<' '<<p.held<<' '<<p.d<<' '<<p.s<<'\n';
    for(auto value:p.points)out<<value<<' ';out<<'\n';require(bool(out),"write public file failed");
}
static void generate(int t,int n,int degree,std::uint64_t seed,const fs::path& output) {
    require(!fs::exists(output),"generation output exists");require(t>=2 && n<=64 && n-9*t>=5 && n-9*t<=31 && degree>=3 && degree<=6,"generator parameters");
    Field f(6,67);std::mt19937_64 rng(seed);Poly g(t+1);g[t]=1;
    do{for(int i=0;i<t;++i)g[i]=Byte(rng()&63);}while(!irreducible(f,g));
    Vec support(64);std::iota(support.begin(),support.end(),Byte(0));for(int i=63;i>0;--i)std::swap(support[i],support[rng()%unsigned(i+1)]);support.resize(n);
    Vec beta(n);for(int j=0;j<n;++j)beta[j]=f.inv(f.power(evaluate(f,g,support[j]),3));
    Matrix check=expand_subfield(f,grs(f,support,3*t,beta));Matrix Y=row_basis(f,nullspace(f,check));
    require(Y.rows==n-9*t,"generated dimension differs from the nominal dimension");
    for(int j=0;j<n;++j)beta[j]=f.inv(f.power(evaluate(f,g,support[j]),4));
    require(same_span(f,check,expand_subfield(f,grs(f,support,4*t,beta))),"wild Goppa identity failed");
    Public p;p.t=t;p.n=n;p.k=Y.rows;p.d=degree;p.s=degree-1;p.points.resize(n);
    for(int j=0;j<n;++j)for(int a=0;a<p.k;++a) {
        int digit=int(std::find(embedding,embedding+4,Y(a,j))-embedding);require(digit<4,"public value is outside GF(4)");p.points[j]|=std::uint64_t(digit)<<(2*a);
    }
    while(p.held<n && p.unit_row(p.held)>=0)++p.held;require(p.held<n && p.s*(n-1)>p.d*p.D(),"generation parameters do not satisfy the interpolation bound");
    fs::create_directories(output/"audit");write_public(p,output/"public.txt");
    std::ofstream sk(output/"audit/secret.txt");sk<<"WGSK1\n3 "<<t<<' '<<n<<" 67\n";
    for(Byte a:support)sk<<int(a)<<' ';sk<<'\n';for(Byte a:g)sk<<int(a)<<' ';sk<<'\n';
    std::ofstream info(output/"generation.json");info<<"{\"seed\":"<<seed<<",\"q\":4,\"m\":3,\"t\":"<<t<<",\"n\":"<<n<<",\"k\":"<<p.k<<",\"degree\":"<<degree<<"}\n";
    std::cout<<"Generated a public ["<<n<<','<<p.k<<"] wild Goppa instance over GF(4).\n";
}
static void audit_secret(const Public& p,const Field& f,const fs::path& path,const Key& key) {
    std::ifstream in(path);std::string magic;int m,t,n,mod;
    require(bool(in>>magic>>m>>t>>n>>mod) && magic=="WGSK1" && m==p.m && t==p.t && n==p.n && mod==p.modulus,"invalid audit key header");
    Vec support(n),g(t+1);for(Vec* v:{&support,&g})for(auto& x:*v){int a;require(bool(in>>a) && a>=0 && a<64,"audit key encoding");x=Byte(a);}
    require(irreducible(f,g),"audit radical is reducible");Vec beta(n);for(int j=0;j<n;++j)beta[j]=f.inv(f.power(evaluate(f,g,support[j]),3));
    Matrix check=expand_subfield(f,grs(f,support,3*t,beta));
    require(rank(f,check)==n-p.k && multiply(f,p.generator(),check.transpose()).zero(),"audit key does not describe public input");
    verify_key(p,f,key);int sigma=support_equivalence(f,support,key.support,p.m);
    std::cout<<"{\"status\":\"passed\",\"global_support_equivalence\":true,\"frobenius_power\":"<<sigma<<"}\n";
}
static void recover(const Public& p,const fs::path& kernel_path,const fs::path& output,const std::string& route) {
    require(route=="both" || route=="minor" || route=="direct","route must be both, minor, or direct");
    require(!fs::exists(output),"recovery output exists; choose a fresh directory");fs::create_directories(output);
    auto started=Clock::now(),tick=started;Field f(6,67);
    std::cout<<"Verify one public holdout panel over GF(4)\n"<<std::flush;
    Kernel K=read_kernel(p,kernel_path);verify_hasse(p,K);double verification=elapsed(tick);
    std::cout<<"Separate three branches using stationary quadrics\n"<<std::flush;
    tick=Clock::now();Local local=recover_branches(p,f,K);local.report.numbers["seconds"]=elapsed(tick);save_matrix(output/"branches.txt",from_rows(local.branches,p.k));
    int depth=route=="minor"?std::max(2,2*p.D()-p.n+1):2*p.D()-1;
    std::cout<<"Continue a coherent jet through order "<<depth<<"\n"<<std::flush;
    Jet jet=continue_jet(p,f,K,local,depth);save_matrix(output/"jet.txt",jet.coefficients);
    Bootstrap boot;Key minor,direct;Matrix minors,sections;
    if(route!="direct") {
        std::cout<<"Bootstrap the branch and finish through the minor code\n"<<std::flush;
        boot=bootstrap(p,f,jet.coefficients);save_matrix(output/"quadrics.txt",boot.forms);save_matrix(output/"tangents.txt",boot.tangents);
        minor=finish_minor(p,f,boot,minors);save_matrix(output/"minor-matrix.txt",minors);save_key(p,minor,output/"minor-key");
    }
    if(route!="minor") {
        std::cout<<"Finish by direct evaluation ratios\n"<<std::flush;
        direct=finish_direct(p,f,jet.coefficients,sections);save_matrix(output/"direct-sections.txt",sections);save_key(p,direct,output/"direct-key");
    }
    if(route=="both")support_equivalence(f,minor.support,direct.support,1);
    const Key& result=route=="minor"?minor:direct;
    std::cout<<"Certify coherence and the global branch equations\n"<<std::flush;
    Report geometry=certify_geometry(p,f,result,local,jet.coefficients,route=="direct"?nullptr:&boot,output);
    Key bad=result;bad.lambda[0]=f.mul(bad.lambda[0],2);bool rejected=false;
    try{verify_key(p,f,bad);}catch(const std::exception&){rejected=true;}require(rejected,"corrupted multiplier accepted");
    std::ofstream out(output/"recovery.json");out<<"{\"status\":\"passed\",\"route\":\""<<route<<"\",\"q\":4,\"m\":3,\"t\":"<<p.t<<",\"n\":"<<p.n<<",\"k\":"<<p.k<<",\"D\":"<<p.D()
        <<",\"heldout\":"<<p.held<<",\"degree\":"<<p.d<<",\"multiplicity\":"<<p.s<<",\"holdout_computations\":1,\"modulus\":67"
        <<",\"kernel_dimension_gf4\":"<<K.full_dimension<<",\"exported_polynomials\":"<<K.rows<<",\"monomials\":"<<K.cols
        <<",\"kernel_verification_seconds\":"<<verification<<",\"forced_zeros\":"<<p.s*(p.n-1)<<",\"composition_degree_bound\":"<<p.d*p.D()
        <<",\"ordinary_interpolation_bound\":true,\"public_inputs_only\":true,\"all_hasse_orders_verified\":true,\"panel_independence_verified\":true,\"branches\":";
    local.report.write(out);out<<",\"jet\":";jet.report.write(out);out<<",\"geometry\":";geometry.write(out);
    if(route!="direct"){out<<",\"bootstrap\":";boot.report.write(out);out<<",\"minor\":";minor.report.write(out);}
    if(route!="minor"){out<<",\"direct\":";direct.report.write(out);}
    if(route=="both")out<<",\"routes_mobius_equivalent\":true";
    out<<",\"tampered_multiplier_rejected\":true,\"total_seconds_excluding_holdout\":"<<elapsed(started)<<"}\n";require(bool(out),"report write failed");
    std::cout<<"PASS: recovered all "<<p.n<<" public coordinates and an equivalent wild Goppa key.\n";
}
static void self_test() {
    Field f(6,67);std::mt19937_64 rng(19);
    for(int a=0;a<4;++a)for(int b=0;b<4;++b)require(f.mul(embedding[a],embedding[b])==embedding[mul4[a][b]],"subfield multiplication");
    const auto& basis=components(f);
    for(int x=0;x<64;++x)require((basis[x][0]^f.mul(basis[x][1],2)^f.mul(basis[x][2],4))==x,"subfield decomposition");
    Kernel K(7,11,7);Matrix embedded(7,11),values(11,9);
    for(int r=0;r<7;++r)for(int c=0;c<11;++c){Byte a=Byte(rng()&3);embedded(r,c)=embedding[a];if(a&1)mzd_write_bit(K.low.get(),r,c,1);if(a&2)mzd_write_bit(K.high.get(),r,c,1);}
    for(auto& x:values.data)x=Byte(rng()&63);require(apply(f,K,values).data==multiply(f,embedded,values).data,"packed GF(4) by GF(64) multiplication");
    Matrix A(7,11);for(auto& x:A.data)x=Byte(rng()&63);Matrix N=nullspace(f,A);
    require(multiply(f,A,N.transpose()).zero() && N.rows==A.cols-rank(f,A),"kernel residual");
    Vec x(11);for(auto& c:x)c=Byte(rng()&63);FixedSolver solver(f,A);require(multiply(f,A,solver.solve(multiply(f,A,x)))==multiply(f,A,x),"fixed solve");
    Vec support={0,1,2,3,4,5,6,7},evaluations;Poly polynomial={7,4,21,3};for(Byte a:support)evaluations.push_back(evaluate(f,polynomial,a));
    require(interpolate(f,support,evaluations)==polynomial,"interpolation");
    constexpr int k=7,d=4,depth=12;Matrix arcs(k,depth+2);auto mons=masks(k,d);
    std::unique_ptr<void,decltype(&series_cache_free)> cache(series_cache_create(k,d,6,67,depth+1),series_cache_free);Vec output(mons.size());
    for(int r=0;r<=depth;++r) {
        Vec jet(k);for(int a=0;a<k;++a){jet[a]=Byte(rng()&63);arcs(a,r)=jet[a];}
        require(series_cache_set(cache.get(),r,jet.data())==0,"cache set");
        for(int order=r;order<=r+1;++order) {
            require(series_cache_eval(cache.get(),order,output.data())==0,"cache evaluation");
            for(int c=0;c<int(mons.size());++c) {
                Poly expected={1};for(int a=0;a<k;++a)if((mons[c]>>a)&1)expected=product(f,expected,arcs.row(a),order+1);
                require(coefficient(expected,order)==output[c],"incremental lookahead differs from independent polynomial multiplication");
            }
        }
    }
    Matrix fresh(int(mons.size()),depth+2);require(series_products(k,d,6,67,depth+1,arcs.data.data(),fresh.data.data())==0,"fresh products");
    for(int c=0;c<int(mons.size());++c){Poly expected={1};for(int a=0;a<k;++a)if((mons[c]>>a)&1)expected=product(f,expected,arcs.row(a),depth+2);for(int r=0;r<depth+2;++r)require(coefficient(expected,r)==fresh(c,r),"fresh product residual");}
    std::cout<<"PASS: subfield arithmetic, packed products, linear algebra, interpolation, and incremental lookahead.\n";
}
int main(int argc,char** argv) {
    try {
        require(argc>=2,"commands: self-test | generate T N DEG SEED DIR | recover PUBLIC KERNEL DIR [both|minor|direct] | verify PUBLIC KEY | audit PUBLIC SECRET KEY");std::string command=argv[1];
        if(command=="self-test"){require(argc==2,"self-test arguments");self_test();}
        else if(command=="generate"){require(argc==7,"generate T N DEG SEED DIR");generate(std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]),std::stoull(argv[5]),argv[6]);}
        else if(command=="recover"){require(argc==5 || argc==6,"recover PUBLIC KERNEL DIR [both|minor|direct]");recover(read_public(argv[2]),argv[3],argv[4],argc==6?argv[5]:"both");}
        else if(command=="verify"){require(argc==4,"verify PUBLIC KEY");Public p=read_public(argv[2]);Field f(6,67);auto r=verify_key(p,f,read_key(p,argv[3]));r.write(std::cout);std::cout<<'\n';}
        else if(command=="audit"){require(argc==5,"audit PUBLIC SECRET KEY");Public p=read_public(argv[2]);Field f(6,67);audit_secret(p,f,argv[3],read_key(p,argv[4]));}
        else throw std::runtime_error("unknown command");return 0;
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
