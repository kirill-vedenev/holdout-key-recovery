#pragma once
#include <algorithm>
#include <cstdint>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace toy {
using Byte = std::uint8_t;
using Vec = std::vector<Byte>;
using Poly = Vec; // Low-degree coefficient first; zero is the empty vector.
inline void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

struct Field {
    int m, q, modulus;
    Vec table, inverse, square_root;
    Field(int degree, int polynomial):m(degree),q(1<<degree),modulus(polynomial),
        table(q*q),inverse(q),square_root(q) {
        require(m>=1 && m<=8 && (modulus>>m)==1,"invalid binary field");
        for(int a=0;a<q;++a) for(int b=0;b<q;++b) {
            int x=a,y=b,c=0;
            while(y) { if(y&1)c^=x; y>>=1; x<<=1; if(x&q)x^=modulus; }
            table[a*q+b]=Byte(c);
        }
        for(int a=1;a<q;++a) {
            inverse[a]=power(Byte(a),q-2);
            require(mul(Byte(a),inverse[a])==1,"field modulus is reducible");
        }
        for(int a=0;a<q;++a) square_root[a]=power(Byte(a),q/2);
    }
    Byte mul(Byte a,Byte b) const { return table[int(a)*q+b]; }
    Byte power(Byte a,unsigned n) const {
        Byte r=1;
        while(n) { if(n&1)r=mul(r,a); a=mul(a,a); n>>=1; }
        return r;
    }
    Byte inv(Byte a) const { require(a!=0,"division by zero"); return inverse[a]; }
    Byte div(Byte a,Byte b) const { return mul(a,inv(b)); }
    Byte sq(Byte a) const { return mul(a,a); }
    Byte root(Byte a) const { return square_root[a]; }
};

struct Matrix {
    int rows=0, cols=0;
    Vec data;
    Matrix()=default;
    Matrix(int r,int c):rows(r),cols(c),data(std::size_t(r)*c) {}
    Byte& operator()(int r,int c) { return data[std::size_t(r)*cols+c]; }
    Byte operator()(int r,int c) const { return data[std::size_t(r)*cols+c]; }
    Vec row(int r) const { return Vec(data.begin()+std::size_t(r)*cols,data.begin()+std::size_t(r+1)*cols); }
    Vec column(int c) const { Vec v(rows); for(int r=0;r<rows;++r)v[r]=(*this)(r,c); return v; }
    void set_row(int r,const Vec& v) {
        require(int(v.size())==cols,"row width mismatch");
        std::copy(v.begin(),v.end(),data.begin()+std::size_t(r)*cols);
    }
    Matrix transpose() const {
        Matrix b(cols,rows);
        for(int r=0;r<rows;++r)for(int c=0;c<cols;++c)b(c,r)=(*this)(r,c);
        return b;
    }
    Matrix first_rows(int n) const {
        require(n>=0 && n<=rows,"row range");
        Matrix b(n,cols);std::copy_n(data.begin(),b.data.size(),b.data.begin());return b;
    }
    bool zero() const { return std::all_of(data.begin(),data.end(),[](Byte a){return a==0;}); }
};
inline Matrix from_rows(const std::vector<Vec>& rows,int cols) {
    Matrix a(int(rows.size()),cols);for(int r=0;r<a.rows;++r)a.set_row(r,rows[r]);return a;
}
inline Matrix stack(const Matrix& a,const Matrix& b) {
    require(a.cols==b.cols,"stack width mismatch");Matrix c(a.rows+b.rows,a.cols);
    std::copy(a.data.begin(),a.data.end(),c.data.begin());
    std::copy(b.data.begin(),b.data.end(),c.data.begin()+a.data.size());return c;
}
inline Matrix multiply(const Field& f,const Matrix& a,const Matrix& b) {
    require(a.cols==b.rows,"product dimensions");Matrix c(a.rows,b.cols);
    for(int r=0;r<a.rows;++r)for(int u=0;u<a.cols;++u)if(Byte x=a(r,u))
        for(int v=0;v<b.cols;++v)c(r,v)^=f.mul(x,b(u,v));
    return c;
}
inline Vec multiply(const Field& f,const Matrix& a,const Vec& b) {
    require(a.cols==int(b.size()),"matrix-vector dimensions");Vec c(a.rows);
    for(int r=0;r<a.rows;++r)for(int u=0;u<a.cols;++u)c[r]^=f.mul(a(r,u),b[u]);return c;
}
inline Byte dot(const Field& f,const Vec& a,const Vec& b) {
    require(a.size()==b.size(),"dot dimensions");Byte c=0;
    for(std::size_t i=0;i<a.size();++i)c^=f.mul(a[i],b[i]);return c;
}
inline std::vector<int> rref(const Field& f,Matrix& a) {
    std::vector<int> pivots;int row=0;
    for(int c=0;c<a.cols && row<a.rows;++c) {
        int p=row;while(p<a.rows && !a(p,c))++p;
        if(p==a.rows)continue;
        if(p!=row)for(int j=0;j<a.cols;++j)std::swap(a(row,j),a(p,j));
        Byte scale=f.inv(a(row,c));
        for(int j=c;j<a.cols;++j)a(row,j)=f.mul(a(row,j),scale);
        for(int r=0;r<a.rows;++r)if(r!=row)if(Byte x=a(r,c))
            for(int j=c;j<a.cols;++j)a(r,j)^=f.mul(x,a(row,j));
        pivots.push_back(c);++row;
    }
    return pivots;
}
inline int rank(const Field& f,Matrix a) { return int(rref(f,a).size()); }
inline Matrix row_basis(const Field& f,Matrix a) { int r=int(rref(f,a).size());return a.first_rows(r); }
inline Matrix nullspace(const Field& f,Matrix a) {
    auto pivots=rref(f,a);std::vector<bool> used(a.cols);
    for(int c:pivots)used[c]=true;
    Matrix k(a.cols-int(pivots.size()),a.cols);int row=0;
    for(int c=0;c<a.cols;++c)if(!used[c]) {
        k(row,c)=1;for(int r=0;r<int(pivots.size());++r)k(row,pivots[r])=a(r,c);++row;
    }
    return k;
}
inline Matrix inverse(const Field& f,const Matrix& a) {
    require(a.rows==a.cols,"inverse of nonsquare matrix");int n=a.rows;Matrix b(n,2*n);
    for(int r=0;r<n;++r) { for(int c=0;c<n;++c)b(r,c)=a(r,c);b(r,n+r)=1; }
    auto p=rref(f,b);require(int(p.size())==n && (n==0 || p.back()==n-1),"singular matrix");
    Matrix result(n,n);for(int r=0;r<n;++r)for(int c=0;c<n;++c)result(r,c)=b(r,n+c);return result;
}
inline bool same_span(const Field& f,const Matrix& a,const Matrix& b) {
    int ra=rank(f,a);return ra==rank(f,b) && ra==rank(f,stack(a,b));
}

struct FixedSolver {
    const Field& f;Matrix a,inv;std::vector<int> columns,rows;
    FixedSolver(const Field& field,const Matrix& input):f(field),a(input) {
        Matrix copy=a;columns=rref(f,copy);
        Matrix sub(a.rows,int(columns.size()));
        for(int r=0;r<a.rows;++r)for(int c=0;c<sub.cols;++c)sub(r,c)=a(r,columns[c]);
        Matrix tr=sub.transpose();rows=rref(f,tr);
        Matrix square(int(rows.size()),int(columns.size()));
        for(int r=0;r<square.rows;++r)for(int c=0;c<square.cols;++c)square(r,c)=sub(rows[r],c);
        inv=inverse(f,square);
    }
    Vec solve(const Vec& rhs) const {
        require(int(rhs.size())==a.rows,"right-hand-side height");Vec selected(rows.size());
        for(int r=0;r<int(rows.size());++r)selected[r]=rhs[rows[r]];
        Vec small=multiply(f,inv,selected),answer(a.cols);
        for(int c=0;c<int(columns.size());++c)answer[columns[c]]=small[c];
        require(multiply(f,a,answer)==rhs,"inconsistent linear system");return answer;
    }
};

inline Poly trim(Poly a) { while(!a.empty() && !a.back())a.pop_back();return a; }
inline Poly add(Poly a,const Poly& b) {
    a.resize(std::max(a.size(),b.size()));for(std::size_t i=0;i<b.size();++i)a[i]^=b[i];return trim(a);
}
inline Poly scale(const Field& f,Poly a,Byte c) { for(auto& x:a)x=f.mul(x,c);return trim(a); }
inline Poly product(const Field& f,const Poly& a,const Poly& b,int precision=-1) {
    if(a.empty() || b.empty())return {};
    int length=int(a.size()+b.size()-1);if(precision>=0)length=std::min(length,precision);
    Poly c(length);
    for(int i=0;i<int(a.size()) && i<length;++i)if(a[i])
        for(int j=0;j<int(b.size()) && i+j<length;++j)c[i+j]^=f.mul(a[i],b[j]);
    return trim(c);
}
inline Byte evaluate(const Field& f,const Poly& a,Byte x) {
    Byte v=0;for(auto i=a.rbegin();i!=a.rend();++i)v=f.mul(v,x)^*i;return v;
}
inline Poly remainder(const Field& f,Poly a,const Poly& b) {
    a=trim(a);require(!b.empty() && b.back(),"zero polynomial divisor");
    while(a.size()>=b.size()) {
        Byte c=f.div(a.back(),b.back());int shift=int(a.size()-b.size());
        for(int i=0;i<int(b.size());++i)a[shift+i]^=f.mul(c,b[i]);a=trim(a);
    }
    return a;
}
inline Poly gcd(const Field& f,Poly a,Poly b) {
    a=trim(a);b=trim(b);while(!b.empty()){Poly r=remainder(f,a,b);a=b;b=r;}
    return a.empty()?a:scale(f,a,f.inv(a.back()));
}
inline Poly square_mod(const Field& f,const Poly& a,const Poly& g) {
    Poly b(2*a.size());for(int i=0;i<int(a.size());++i)b[2*i]=f.sq(a[i]);return remainder(f,trim(b),g);
}
inline bool irreducible(const Field& f,const Poly& g) {
    if(g.size()<2 || !g.back())return false;
    int degree=int(g.size())-1;Poly x={0,1},h=remainder(f,x,g);
    std::vector<int> checks;
    for(int p=2,left=degree;p<=degree;++p)if(left%p==0) {
        checks.push_back(degree/p);while(left%p==0)left/=p;
    }
    for(int r=1;r<=degree;++r) {
        for(int bit=0;bit<f.m;++bit)h=square_mod(f,h,g);
        if(std::find(checks.begin(),checks.end(),r)!=checks.end() && gcd(f,add(h,x),g)!=Poly{1})return false;
    }
    return add(h,remainder(f,x,g)).empty();
}
inline Vec support_derivatives(const Field& f,const Vec& support) {
    Vec dp(support.size(),1);
    for(int i=0;i<int(support.size());++i)for(int j=0;j<i;++j) {
        Byte x=support[i]^support[j];require(x!=0,"duplicate support");
        dp[i]=f.mul(dp[i],x);dp[j]=f.mul(dp[j],x);
    }
    return dp;
}
inline Poly interpolate(const Field& f,const Vec& support,const Vec& values) {
    require(support.size()==values.size(),"interpolation dimensions");
    Poly pi={1};for(Byte a:support)pi=product(f,pi,Poly{a,1});
    Vec dp=support_derivatives(f,support);Poly result(support.size());int n=int(support.size());
    for(int j=0;j<n;++j) {
        Poly quotient(n);quotient[n-1]=pi[n];
        for(int i=n-1;i>0;--i)quotient[i-1]=pi[i]^f.mul(support[j],quotient[i]);
        require((pi[0]^f.mul(support[j],quotient[0]))==0,"interpolation division");
        Byte c=f.div(values[j],dp[j]);for(int i=0;i<n;++i)result[i]^=f.mul(c,quotient[i]);
    }
    return trim(result);
}
inline Poly compose(const Field& f,const Poly& a,const Poly& b,int precision) {
    Poly out;for(auto i=a.rbegin();i!=a.rend();++i)out=add(product(f,out,b,precision),Poly{*i});return out;
}
inline Byte coefficient(const Poly& p,int r) { return r<int(p.size())?p[r]:0; }

inline std::vector<std::uint64_t> binary_rref(std::vector<std::uint64_t> a,int width) {
    int row=0;for(int c=0;c<width && row<int(a.size());++c) {
        int p=row;while(p<int(a.size()) && !((a[p]>>c)&1))++p;if(p==int(a.size()))continue;
        std::swap(a[p],a[row]);
        for(int r=0;r<int(a.size());++r)if(r!=row && ((a[r]>>c)&1))a[r]^=a[row];++row;
    }
    a.resize(row);return a;
}
inline std::vector<std::uint64_t> binary_kernel(std::vector<std::uint64_t> a,int width) {
    a=binary_rref(a,width);std::vector<int> pivots;std::uint64_t used=0;
    for(auto row:a){int c=__builtin_ctzll(row);pivots.push_back(c);used|=std::uint64_t(1)<<c;}
    std::vector<std::uint64_t> out;
    for(int c=0;c<width;++c)if(!((used>>c)&1)) {
        std::uint64_t row=std::uint64_t(1)<<c;
        for(int r=0;r<int(a.size());++r)if((a[r]>>c)&1)row|=std::uint64_t(1)<<pivots[r];out.push_back(row);
    }
    return out;
}
inline std::vector<std::uint64_t> masks(int k,int d,std::uint64_t allowed=~std::uint64_t(0)) {
    std::vector<std::uint64_t> out;
    std::function<void(int,int,std::uint64_t)> visit=[&](int a,int left,std::uint64_t mask) {
        if(!left){out.push_back(mask);return;}
        for(int b=a;b<=k-left;++b)if((allowed>>b)&1)visit(b+1,left-1,mask|(std::uint64_t(1)<<b));
    };
    visit(0,d,0);return out;
}
} // namespace toy
