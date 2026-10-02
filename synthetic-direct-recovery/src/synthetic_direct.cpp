// Synthetic test of direct support recovery from one local jet.
//
// A key is a subfield subcode C = GRS_{D+1}(alpha, lambda) ∩ F_q^n with
// K = GF(q^m). The program generates keys of three families (binary Goppa,
// wild Goppa, generic alternant), computes a local jet of the secret curve at
// one public position, with a random reparametrization and on a random
// Frobenius conjugate branch, and then runs the direct recovery using only the
// public generator matrix and this jet:
//
//   1. find forms R, S of degree delta with
//        R(F(Z)) = c (Z - a)^{delta D},
//        S(F(Z)) = (Z - a)^{delta D - 1} (u Z + v),  u a + v != 0,
//      from linear conditions on the jet (binary Goppa codes: delta = 2,
//      orders 2D and 2D - 2, and a square root of the quotient);
//   2. evaluate beta_j = R(y_j) / S(y_j), a Moebius image of the support;
//   3. find a multiplier z with C diag(z) H_RS^T = 0 and set lambda' = 1/z;
//   4. check that GRS_{D+1}(alpha', lambda') ∩ F_q^n equals C.
//
// The secret key is used only to produce the public matrix and the jet, and
// for a final audit of the recovered support.

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

[[noreturn]] void die(const std::string& msg) {
  std::fprintf(stderr, "error: %s\n", msg.c_str());
  std::exit(2);
}

double now_s() {
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

class Rng {
 public:
  explicit Rng(const std::vector<u32>& s) {
    std::seed_seq q(s.begin(), s.end());
    g_.seed(q);
  }
  u64 next() { return g_(); }
  u64 below(u64 b) {
    const u64 threshold = (0 - b) % b;
    for (;;) {
      const u64 x = g_();
      if (x >= threshold) return x % b;
    }
  }
  template <class T>
  void shuffle(std::vector<T>& v) {
    for (std::size_t i = v.size(); i > 1; --i) std::swap(v[i - 1], v[below(i)]);
  }

 private:
  std::mt19937_64 g_;
};

// ------------------------------------------------------------- GF(p^e)
//
// Elements are integers whose base-p digits are the coefficients of a
// polynomial basis. Multiplication uses logarithms; addition in odd
// characteristic uses Zech logarithms.

struct Field {
  u32 p = 0, e = 0, size = 0, Q = 0, ZL = 0, neg1 = 0, gen = 0;
  bool char2 = false;
  std::vector<u32> modulus;  // monic, degree e, coefficients 0..e
  std::vector<u32> ex;       // ex[i] = gen^(i mod Q) for i < 2Q, zero above
  std::vector<u32> lg;       // lg[0] = ZL = 2Q
  std::vector<u32> zech;     // zech[d] = log(1 + gen^d), ZL when 1 + gen^d = 0

  std::vector<u32> digits(u32 a) const {
    std::vector<u32> d(e);
    for (u32 i = 0; i < e; ++i) { d[i] = a % p; a /= p; }
    return d;
  }
  u32 encode(const std::vector<u32>& d) const {
    u32 a = 0;
    for (u32 i = e; i-- > 0;) a = a * p + d[i];
    return a;
  }
  u32 slow_mul(u32 a, u32 b) const {
    const std::vector<u32> A = digits(a), B = digits(b);
    std::vector<u32> P(2 * e - 1, 0);
    for (u32 i = 0; i < e; ++i) {
      for (u32 j = 0; j < e; ++j) P[i + j] = (P[i + j] + A[i] * B[j]) % p;
    }
    for (u32 d = 2 * e - 1; d-- > e;) {
      const u32 c = P[d];
      if (!c) continue;
      for (u32 i = 0; i <= e; ++i) P[d - e + i] = (P[d - e + i] + (p - c) * modulus[i]) % p;
    }
    P.resize(e);
    return encode(P);
  }
  u32 slow_pow(u32 a, u64 k) const {
    u32 r = 1;
    while (k) {
      if (k & 1) r = slow_mul(r, a);
      a = slow_mul(a, a);
      k >>= 1;
    }
    return r;
  }

  Field(u32 p_, u32 e_) : p(p_), e(e_) {
    if (p < 2) die("bad characteristic");
    u64 s = 1;
    for (u32 i = 0; i < e; ++i) s *= p;
    if (s > 65536 || s < 2) die("field size must be at most 2^16");
    size = static_cast<u32>(s);
    Q = size - 1;
    ZL = 2 * Q;
    char2 = p == 2;
    std::vector<u32> primes;
    u32 x = Q;
    for (u32 r = 2; r * r <= x; ++r) {
      if (x % r == 0) {
        primes.push_back(r);
        while (x % r == 0) x /= r;
      }
    }
    if (x > 1) primes.push_back(x);
    auto primitive = [&](u32 c) {
      if (c == 0 || slow_pow(c, Q) != 1) return false;
      for (u32 r : primes) {
        if (slow_pow(c, Q / r) == 1) return false;
      }
      return true;
    };
    // Classic McEliece moduli for m = 12, 13; otherwise the first monic
    // polynomial whose root x is primitive.
    modulus.assign(e + 1, 0);
    modulus[e] = 1;
    bool found = false;
    if (p == 2 && (e == 12 || e == 13)) {
      const u32 bits = e == 12 ? 0x1009 : 0x201B;
      for (u32 i = 0; i <= e; ++i) modulus[i] = bits >> i & 1;
      for (u32 c = 2; c < size && !found; ++c) {
        if (primitive(c)) { gen = c; found = true; }
      }
      if (!found) die("modulus not irreducible");
    } else {
      const u32 xel = e == 1 ? 0 : p;  // the class of x
      for (u32 low = 1; low < size && !found; ++low) {
        u32 v = low;
        for (u32 i = 0; i < e; ++i) { modulus[i] = v % p; v /= p; }
        if (e == 1) {
          for (u32 c = 2; c < size && !found; ++c) {
            if (primitive(c)) { gen = c; found = true; }
          }
          if (size == 2) { gen = 1; found = true; }
        } else if (primitive(xel)) {
          gen = xel;
          found = true;
        }
      }
      if (!found) die("no primitive polynomial found");
    }
    ex.assign(4 * static_cast<std::size_t>(Q) + 1, 0);
    lg.assign(size, ZL);
    u32 v = 1;
    for (u32 i = 0; i < Q; ++i) {
      if (lg[v] != ZL) die("generator has small order");
      ex[i] = v;
      ex[i + Q] = v;
      lg[v] = i;
      v = slow_mul(v, gen);
    }
    zech.assign(Q, ZL);
    for (u32 d = 0; d < Q; ++d) {
      const u32 a = ex[d];
      const u32 c0 = a % p;
      zech[d] = lg[a - c0 + (c0 + 1) % p];
    }
    neg1 = char2 ? 0 : lg[p - 1];
  }

  u32 mul(u32 a, u32 b) const { return ex[lg[a] + lg[b]]; }
  u32 add(u32 a, u32 b) const {
    if (char2) return a ^ b;
    if (!a) return b;
    if (!b) return a;
    const u32 la = lg[a], lb = lg[b];
    const u32 d = lb >= la ? lb - la : lb + Q - la;
    const u32 z = zech[d];
    return z == ZL ? 0 : ex[la + z];
  }
  u32 neg(u32 a) const { return (char2 || !a) ? a : ex[lg[a] + neg1]; }
  u32 sub(u32 a, u32 b) const { return add(a, neg(b)); }
  u32 inv(u32 a) const {
    if (!a) die("inverse of zero");
    return ex[Q - lg[a]];
  }
  u32 div(u32 a, u32 b) const { return mul(a, inv(b)); }
  u32 pow(u32 a, u64 k) const {
    if (k == 0) return 1;
    if (!a) return 0;
    return ex[(static_cast<u64>(lg[a]) * (k % Q)) % Q];
  }
  // v[x] += c * w[x] for x in [from, len), with w given by logarithms.
  void axpy_log(u32* v, const u32* wlog, u32 c, int from, int len) const {
    if (!c) return;
    const u32 lc = lg[c];
    const u32* E = ex.data();
    if (char2) {
      for (int x = from; x < len; ++x) v[x] ^= E[wlog[x] + lc];
    } else {
      for (int x = from; x < len; ++x) {
        const u32 t = E[wlog[x] + lc];
        if (t) v[x] = add(v[x], t);
      }
    }
  }
};

// The subfield F_q of K, q = p^s, and F_q-coordinates of K.
struct Subfield {
  u32 q = 0, s = 0, m = 0;
  std::vector<u32> elems;  // elems[0] = 0, elems[1] = 1, then the other elements
  std::vector<u32> coord;  // coord[x] = sum_i d_i q^i with x = sum_i elems[d_i] theta^i
  std::vector<u32> index;  // index[x] = i with elems[i] = x, for x in F_q

  Subfield(const Field& K, u32 s_) : s(s_) {
    if (K.e % s) die("subfield degree must divide the field degree");
    q = 1;
    for (u32 i = 0; i < s; ++i) q *= K.p;
    m = K.e / s;
    elems.push_back(0);
    const u32 step = K.Q / (q - 1);
    for (u32 i = 0; i < q - 1; ++i) elems.push_back(K.ex[i * step]);
    index.assign(K.size, 0xFFFFFFFFu);
    for (u32 i = 0; i < q; ++i) index[elems[i]] = i;
    // theta = gen has degree m over F_q.
    std::vector<u32> th(m);
    th[0] = 1;
    for (u32 i = 1; i < m; ++i) th[i] = K.mul(th[i - 1], K.gen);
    coord.assign(K.size, 0xFFFFFFFFu);
    std::vector<u32> d(m, 0);
    for (u32 idx = 0; idx < K.size; ++idx) {
      u32 r = idx, x = 0;
      for (u32 i = 0; i < m; ++i) {
        x = K.add(x, K.mul(elems[r % q], th[i]));
        r /= q;
      }
      if (coord[x] != 0xFFFFFFFFu) die("subfield basis is not a basis");
      coord[x] = idx;
    }
  }
  bool contains(u32 x) const { return index[x] != 0xFFFFFFFFu; }
  u32 component(u32 x, u32 c) const {
    u32 r = coord[x];
    for (u32 i = 0; i < c; ++i) r /= q;
    return elems[r % q];
  }
  // phi^sigma(x) = x^(q^sigma), the Frobenius automorphism over F_q.
  u32 frob(const Field& K, u32 x, u32 sigma) const {
    u64 k = 1;
    for (u32 i = 0; i < sigma; ++i) k *= q;
    return K.pow(x, k);
  }
};

// ----------------------------------------------------------- polynomials

using Poly = std::vector<u32>;

void trim(Poly& a) {
  while (!a.empty() && a.back() == 0) a.pop_back();
}
int deg(const Poly& a) { return static_cast<int>(a.size()) - 1; }

void pmod(const Field& K, Poly& a, const Poly& g) {  // g monic
  const int t = deg(g);
  for (int d = deg(a); d >= t; --d) {
    const u32 c = a[d];
    if (!c) continue;
    for (int j = 0; j <= t; ++j) a[d - t + j] = K.sub(a[d - t + j], K.mul(c, g[j]));
  }
  if (deg(a) >= t) a.resize(t);
  trim(a);
}

Poly monic(const Field& K, Poly a) {
  trim(a);
  if (a.empty()) return a;
  const u32 il = K.inv(a.back());
  for (u32& c : a) c = K.mul(c, il);
  return a;
}

Poly pgcd(const Field& K, Poly a, Poly b) {
  trim(a);
  trim(b);
  while (!b.empty()) {
    b = monic(K, b);
    pmod(K, a, b);
    std::swap(a, b);
  }
  return monic(K, a);
}

u32 peval(const Field& K, const Poly& a, u32 x) {
  u32 r = 0;
  for (int i = deg(a); i >= 0; --i) r = K.add(K.mul(r, x), a[i]);
  return r;
}

// a^p mod g in characteristic p: sum c_i^p x^(i p).
Poly frob_mod(const Field& K, const Poly& a, const Poly& g) {
  Poly r(a.empty() ? 0 : (a.size() - 1) * K.p + 1, 0);
  for (std::size_t i = 0; i < a.size(); ++i) r[i * K.p] = K.pow(a[i], K.p);
  trim(r);
  pmod(K, r, g);
  return r;
}

// Ben-Or irreducibility test over K.
bool irreducible(const Field& K, const Poly& g) {
  const int t = deg(g);
  if (t <= 0) return false;
  if (t == 1) return true;
  Poly h = {0, 1};
  for (int i = 1; i <= t / 2; ++i) {
    for (u32 s = 0; s < K.e; ++s) h = frob_mod(K, h, g);  // h = x^(|K|^i)
    Poly d = h;
    if (d.size() < 2) d.resize(2, 0);
    d[1] = K.sub(d[1], 1);
    trim(d);
    if (deg(pgcd(K, g, d)) > 0) return false;
  }
  return true;
}

// Coefficients of a(beta + X).
Poly taylor_shift(const Field& K, Poly a, u32 beta) {
  const int n = static_cast<int>(a.size());
  for (int i = 0; i < n; ++i) {
    for (int j = n - 2; j >= i; --j) a[j] = K.add(a[j], K.mul(beta, a[j + 1]));
  }
  return a;
}

// Newton interpolation at fixed nodes, for many value vectors.
class Interpolator {
 public:
  Interpolator(const Field& K, std::vector<u32> xs) : K_(&K), xs_(std::move(xs)) {
    n_ = static_cast<int>(xs_.size());
    inv_.assign(static_cast<std::size_t>(n_) * n_, 0);
    for (int j = 1; j < n_; ++j) {
      for (int i = j; i < n_; ++i) inv_[static_cast<std::size_t>(j) * n_ + i] = K.inv(K.sub(xs_[i], xs_[i - j]));
    }
  }
  Poly operator()(std::vector<u32> c) const {
    const Field& K = *K_;
    for (int j = 1; j < n_; ++j) {
      const u32* iv = &inv_[static_cast<std::size_t>(j) * n_];
      for (int i = n_ - 1; i >= j; --i) c[i] = K.mul(K.sub(c[i], c[i - 1]), iv[i]);
    }
    Poly p(n_, 0);
    int len = 1;
    p[0] = c[n_ - 1];
    for (int i = n_ - 2; i >= 0; --i) {
      const u32 mx = K.neg(xs_[i]);
      for (int d = len; d >= 1; --d) p[d] = K.add(p[d - 1], K.mul(p[d], mx));
      p[0] = K.add(K.mul(p[0], mx), c[i]);
      ++len;
    }
    trim(p);
    return p;
  }

 private:
  const Field* K_;
  std::vector<u32> xs_;
  int n_ = 0;
  std::vector<u32> inv_;
};

// ------------------------------------------------------- linear algebra

// Semi-echelon basis of rows of length `len`; row r is stored by
// logarithms, vanishes before its pivot col[r], equals 1 there, and vanishes
// at the pivots of earlier rows. Columns >= `vars` hold right-hand sides.
struct Echelon {
  const Field* K;
  int len = 0, vars = 0;
  std::vector<u32> rows;  // logarithms
  std::vector<int> col;
  int rank = 0;
  bool inconsistent = false;

  Echelon(const Field& f, int length, int nvars) : K(&f), len(length), vars(nvars) {}

  void reduce(u32* v) const {
    for (int r = 0; r < rank; ++r) {
      const int c = col[r];
      const u32 x = v[c];
      if (!x) continue;
      K->axpy_log(v, &rows[static_cast<std::size_t>(r) * len], K->neg(x), c, len);
    }
  }
  // Returns true if the row increased the rank.
  bool insert(std::vector<u32> v) {
    reduce(v.data());
    int c = 0;
    while (c < vars && !v[c]) ++c;
    if (c == vars) {
      for (int x = vars; x < len; ++x) {
        if (v[x]) inconsistent = true;
      }
      return false;
    }
    const u32 li = K->Q - K->lg[v[c]];
    rows.resize(static_cast<std::size_t>(rank + 1) * len);
    u32* p = &rows[static_cast<std::size_t>(rank) * len];
    for (int x = 0; x < len; ++x) p[x] = v[x] ? (K->lg[v[x]] + li) % K->Q : K->ZL;
    col.push_back(c);
    ++rank;
    return true;
  }
  // One solution of the system for right-hand side column `rhs` (free
  // variables zero). Assumes the system is consistent.
  std::vector<u32> solve(int rhs) const {
    std::vector<u32> h(vars, 0);
    for (int r = rank - 1; r >= 0; --r) {
      const u32* p = &rows[static_cast<std::size_t>(r) * len];
      u32 acc = K->ex[p[rhs]];
      for (int c = col[r] + 1; c < vars; ++c) {
        if (h[c] && p[c] != K->ZL) acc = K->sub(acc, K->mul(K->ex[p[c]], h[c]));
      }
      h[col[r]] = acc;
    }
    return h;
  }
  // A nonzero kernel vector when the kernel is one-dimensional.
  std::vector<u32> kernel_vector() const {
    std::vector<char> piv(vars, 0);
    for (int r = 0; r < rank; ++r) piv[col[r]] = 1;
    int f = -1;
    for (int c = 0; c < vars; ++c) {
      if (!piv[c]) { f = c; break; }
    }
    std::vector<u32> h(vars, 0);
    if (f < 0) return h;
    h[f] = 1;
    for (int r = rank - 1; r >= 0; --r) {
      const u32* p = &rows[static_cast<std::size_t>(r) * len];
      u32 acc = 0;
      for (int c = col[r] + 1; c < vars; ++c) {
        if (h[c] && p[c] != K->ZL) acc = K->sub(acc, K->mul(K->ex[p[c]], h[c]));
      }
      h[col[r]] = acc;
    }
    return h;
  }
};

// Kernel over F_q of a matrix with entries in F_q (stored as K elements),
// in reduced echelon form; the basis is systematic on the free columns.
std::vector<std::vector<u32>> fq_kernel(const Field& K, std::vector<std::vector<u32>> A, int n) {
  int r = 0;
  std::vector<int> pivcol;
  const int rows = static_cast<int>(A.size());
  for (int c = 0; c < n && r < rows; ++c) {
    int p = -1;
    for (int i = r; i < rows; ++i) {
      if (A[i][c]) { p = i; break; }
    }
    if (p < 0) continue;
    std::swap(A[p], A[r]);
    const u32 il = K.inv(A[r][c]);
    for (int x = 0; x < n; ++x) A[r][x] = K.mul(A[r][x], il);
    for (int i = 0; i < rows; ++i) {
      if (i == r || !A[i][c]) continue;
      const u32 f = K.neg(A[i][c]);
      for (int x = 0; x < n; ++x) {
        if (A[r][x]) A[i][x] = K.add(A[i][x], K.mul(f, A[r][x]));
      }
    }
    pivcol.push_back(c);
    ++r;
  }
  std::vector<char> isp(n, 0);
  for (int c : pivcol) isp[c] = 1;
  std::vector<std::vector<u32>> basis;
  for (int f = 0; f < n; ++f) {
    if (isp[f]) continue;
    std::vector<u32> v(n, 0);
    v[f] = 1;
    for (int i = 0; i < r; ++i) v[pivcol[i]] = K.neg(A[i][f]);
    basis.push_back(v);
  }
  return basis;
}

// Rank over F_q of a matrix with entries in F_q.
int fq_rank(const Field& K, std::vector<std::vector<u32>> A, int n) {
  return n - static_cast<int>(fq_kernel(K, std::move(A), n).size());
}

// Expands a K-matrix row-wise over F_q: each row gives m rows of components.
std::vector<std::vector<u32>> expand(const Field& K, const Subfield& Fq,
                                     const std::vector<std::vector<u32>>& H, int n) {
  (void)K;
  std::vector<std::vector<u32>> out;
  for (const auto& row : H) {
    for (u32 c = 0; c < Fq.m; ++c) {
      std::vector<u32> v(n);
      for (int j = 0; j < n; ++j) v[j] = Fq.component(row[j], c);
      out.push_back(v);
    }
  }
  return out;
}

// Pi'(a_j) = prod_{l != j} (a_j - a_l).
std::vector<u32> support_derivative(const Field& K, const std::vector<u32>& a) {
  const int n = static_cast<int>(a.size());
  std::vector<u32> d(n);
  for (int j = 0; j < n; ++j) {
    u64 s = 0;
    for (int l = 0; l < n; ++l) {
      if (l != j) s += K.lg[K.sub(a[j], a[l])];
    }
    d[j] = K.ex[s % K.Q];
  }
  return d;
}

// Parity-check matrix of GRS_{D+1}(a, lam): rows a^l / (lam_j Pi'(a_j)).
std::vector<std::vector<u32>> grs_parity(const Field& K, const std::vector<u32>& a,
                                         const std::vector<u32>& lam, int D) {
  const int n = static_cast<int>(a.size());
  const std::vector<u32> d = support_derivative(K, a);
  std::vector<std::vector<u32>> H(n - D - 1, std::vector<u32>(n));
  for (int j = 0; j < n; ++j) {
    u32 v = K.inv(K.mul(lam[j], d[j]));
    for (int l = 0; l < n - D - 1; ++l) {
      H[l][j] = v;
      v = K.mul(v, a[j]);
    }
  }
  return H;
}

// ----------------------------------------------------------------- keys

enum class Family { BinaryGoppa, WildGoppa, Alternant };

const char* family_name(Family f) {
  switch (f) {
    case Family::BinaryGoppa: return "binary-goppa";
    case Family::WildGoppa: return "wild-goppa";
    default: return "alternant";
  }
}

struct Params {
  Family family = Family::Alternant;
  u32 q = 2, m = 8;
  int n = 0, t = 0, r = 0;
  int delta = 0;  // 0 = smallest admissible
};

struct Key {
  int n = 0, k = 0, D = 0;
  Poly g;
  std::vector<u32> alpha, lambda;    // C ⊆ GRS_{D+1}(alpha, lambda)
  std::vector<std::vector<u32>> Y;   // k x n public generator, systematic
};

bool make_key(const Field& K, const Subfield& Fq, const Params& P, u64 seed, Key& key,
              std::string& why) {
  key = Key();
  const int n = P.n;
  key.n = n;
  Rng rng({static_cast<u32>(P.family), P.q, P.m, static_cast<u32>(n), static_cast<u32>(P.t),
           static_cast<u32>(P.r), static_cast<u32>(seed), static_cast<u32>(seed >> 32), 0x5D1Eu});
  const bool goppa = P.family != Family::Alternant;
  if (goppa) {
    if (P.t < 1) { why = "Goppa degree t must be positive"; return false; }
    for (int tries = 0;; ++tries) {
      if (tries > 100000) { why = "no irreducible Goppa polynomial found"; return false; }
      Poly g(P.t + 1);
      for (int i = 0; i < P.t; ++i) g[i] = static_cast<u32>(rng.below(K.size));
      g[P.t] = 1;
      if (irreducible(K, g)) { key.g = g; break; }
    }
  }
  std::vector<u32> cand;
  for (u32 x = 0; x < K.size; ++x) {
    if (!goppa || peval(K, key.g, x) != 0) cand.push_back(x);
  }
  if (static_cast<int>(cand.size()) < n) { why = "support larger than the field"; return false; }
  rng.shuffle(cand);
  key.alpha.assign(cand.begin(), cand.begin() + n);
  const std::vector<u32> dpi = support_derivative(K, key.alpha);

  std::vector<std::vector<u32>> H;
  if (P.family == Family::BinaryGoppa) {
    if (P.q != 2) { why = "binary Goppa codes need q = 2"; return false; }
    // Gamma(alpha, G) = Gamma(alpha, G^2) ⊆ GRS_{n-2t}(alpha, G(a)^2 / Pi'(a)).
    H.assign(P.t, std::vector<u32>(n));
    key.D = n - 2 * P.t - 1;
    key.lambda.resize(n);
    for (int j = 0; j < n; ++j) {
      const u32 gv = peval(K, key.g, key.alpha[j]);
      u32 v = K.inv(gv);
      for (int l = 0; l < P.t; ++l) { H[l][j] = v; v = K.mul(v, key.alpha[j]); }
      key.lambda[j] = K.div(K.mul(gv, gv), dpi[j]);
    }
  } else if (P.family == Family::WildGoppa) {
    if (P.q < 3) { why = "wild Goppa codes need q > 2"; return false; }
    // Gamma(alpha, g^(q-1)) = Gamma(alpha, g^q) ⊆ GRS_{n-qt}(alpha, g(a)^q / Pi'(a)).
    const int rows = static_cast<int>(P.q - 1) * P.t;
    H.assign(rows, std::vector<u32>(n));
    key.D = n - static_cast<int>(P.q) * P.t - 1;
    key.lambda.resize(n);
    for (int j = 0; j < n; ++j) {
      const u32 gv = peval(K, key.g, key.alpha[j]);
      u32 v = K.inv(K.pow(gv, P.q - 1));
      for (int l = 0; l < rows; ++l) { H[l][j] = v; v = K.mul(v, key.alpha[j]); }
      key.lambda[j] = K.div(K.pow(gv, P.q), dpi[j]);
    }
  } else {
    if (P.r < 1) { why = "alternant codes need r >= 1"; return false; }
    key.D = n - P.r - 1;
    key.lambda.resize(n);
    for (int j = 0; j < n; ++j) key.lambda[j] = 1 + static_cast<u32>(rng.below(K.Q));
    H = grs_parity(K, key.alpha, key.lambda, key.D);
  }
  if (key.D < 2) { why = "degree bound D too small"; return false; }
  key.Y = fq_kernel(K, expand(K, Fq, H, n), n);
  key.k = static_cast<int>(key.Y.size());
  if (key.k < 2) { why = "code dimension below 2"; return false; }
  return true;
}

// Curve polynomials f_a of degree at most D with y_j = lambda_j F(alpha_j).
bool curve_polys(const Field& K, const Key& key, std::vector<Poly>& f) {
  const int L = key.D + 1;
  Interpolator interp(K, std::vector<u32>(key.alpha.begin(), key.alpha.begin() + L));
  f.assign(key.k, Poly());
  for (int a = 0; a < key.k; ++a) {
    std::vector<u32> ys(L);
    for (int j = 0; j < L; ++j) ys[j] = K.div(key.Y[a][j], key.lambda[j]);
    f[a] = interp(ys);
    if (deg(f[a]) > key.D) return false;
    for (int j = L; j < key.n; ++j) {
      if (peval(K, f[a], key.alpha[j]) != K.div(key.Y[a][j], key.lambda[j])) return false;
    }
  }
  return true;
}

// --------------------------------------------------------------- jets

struct Jet {
  int pos = 0;
  u32 sigma = 0;
  int len = 0;                     // coefficients 0..len-1
  std::vector<std::vector<u32>> V;  // k series
};

// First `len` coefficients of u(T) phi^sigma(lambda_i) F^sigma(phi^sigma(alpha_i) + psi(T)).
Jet make_jet(const Field& K, const Subfield& Fq, const Key& key, const std::vector<Poly>& f,
             int len, Rng& rng) {
  Jet J;
  J.len = len;
  J.pos = static_cast<int>(rng.below(key.n));
  J.sigma = static_cast<u32>(rng.below(Fq.m));
  const u32 beta = Fq.frob(K, key.alpha[J.pos], J.sigma);
  const u32 lam = Fq.frob(K, key.lambda[J.pos], J.sigma);
  const int deg_aux = 4;
  std::vector<u32> u(deg_aux + 1), psi(deg_aux + 1);
  u[0] = 1;
  psi[0] = 0;
  psi[1] = 1 + static_cast<u32>(rng.below(K.Q));
  for (int i = 1; i <= deg_aux; ++i) u[i] = static_cast<u32>(rng.below(K.size));
  for (int i = 2; i <= deg_aux; ++i) psi[i] = static_cast<u32>(rng.below(K.size));
  // Powers psi^e mod T^len for e <= min(D, len - 1).
  const int emax = std::min(key.D, len - 1);
  std::vector<std::vector<u32>> pw(emax + 1, std::vector<u32>(len, 0));
  pw[0][0] = 1;
  for (int e = 1; e <= emax; ++e) {
    for (int i = 0; i < len; ++i) {
      const u32 c = pw[e - 1][i];
      if (!c) continue;
      for (int d = 1; d <= deg_aux && i + d < len; ++d) {
        if (psi[d]) pw[e][i + d] = K.add(pw[e][i + d], K.mul(c, psi[d]));
      }
    }
  }
  J.V.assign(key.k, std::vector<u32>(len, 0));
  std::vector<u32> s(len);
  for (int a = 0; a < key.k; ++a) {
    Poly fa = f[a];
    for (u32& c : fa) c = Fq.frob(K, c, J.sigma);
    const Poly g = taylor_shift(K, fa, beta);
    std::fill(s.begin(), s.end(), 0);
    for (int e = 0; e <= emax && e < static_cast<int>(g.size()); ++e) {
      if (!g[e]) continue;
      for (int i = e; i < len; ++i) {
        if (pw[e][i]) s[i] = K.add(s[i], K.mul(g[e], pw[e][i]));
      }
    }
    std::vector<u32>& v = J.V[a];
    for (int i = 0; i < len; ++i) {
      if (!s[i]) continue;
      const u32 c = K.mul(lam, s[i]);
      for (int d = 0; d <= deg_aux && i + d < len; ++d) {
        if (u[d]) v[i + d] = K.add(v[i + d], K.mul(c, u[d]));
      }
    }
  }
  return J;
}

// ----------------------------------------------------------- recovery

struct Recovery {
  int delta = 0, monomials_total = 0, monomials_used = 0;
  bool R_found = false, S_found = false;
  int infinities = 0;
  bool distinct = false;
  int multiplier_nullity = -1;
  bool multiplier_nonzero = false;
  bool contains = false;  // C ⊆ GRS_{D+1}(alpha', lambda')
  bool equal = false;     // GRS_{D+1}(alpha', lambda') ∩ F_q^n = C
  std::string failure;
  std::vector<u32> alpha, lambda;
  double t_forms = 0, t_support = 0, t_multiplier = 0, t_verify = 0;
};

u64 binom(int a, int b) {
  if (b < 0 || b > a) return 0;
  u64 r = 1;
  for (int i = 1; i <= b; ++i) r = r * static_cast<u64>(a - b + i) / static_cast<u64>(i);
  return r;
}

// Degree-delta monomials as nondecreasing index tuples; pure powers first.
std::vector<std::vector<int>> monomials(int k, int delta, std::size_t need, Rng& rng,
                                        std::size_t& total) {
  total = static_cast<std::size_t>(binom(k + delta - 1, delta));
  std::vector<std::vector<int>> pure, rest;
  std::vector<int> cur(delta, 0);
  std::function<void(int, int)> rec = [&](int pos, int from) {
    if (pos == delta) {
      bool p = true;
      for (int i = 1; i < delta; ++i) p = p && cur[i] == cur[0];
      (p ? pure : rest).push_back(cur);
      return;
    }
    for (int a = from; a < k; ++a) { cur[pos] = a; rec(pos + 1, a); }
  };
  if (total > 20000000ULL) die("too many monomials");
  rec(0, 0);
  rng.shuffle(pure);
  rng.shuffle(rest);
  std::vector<std::vector<int>> out = pure;
  for (std::size_t i = 0; i < rest.size() && out.size() < need; ++i) out.push_back(rest[i]);
  return out;
}

// Direct recovery from the public generator Y and the jet V at position pos.
Recovery direct_recover(const Field& K, const Subfield& Fq, const std::vector<std::vector<u32>>& Y,
                        int n, int D, const Jet& J, bool binary, int delta, Rng& rng) {
  Recovery out;
  out.delta = delta;
  const int k = static_cast<int>(Y.size());
  const int pos = J.pos;
  // Orders: R needs [T^e] R(V) = 0 for e < LR; S needs zero below LS and 1 at LS.
  const int LR = delta * D;
  const int LS = binary ? 2 * D - 2 : delta * D - 1;
  if (J.len < LR) die("jet too short");
  const double t0 = now_s();
  // Normalization column for R: a public position other than the jet position.
  int j0 = static_cast<int>(rng.below(n - 1));
  if (j0 >= pos) ++j0;

  std::size_t total = 0;
  std::size_t need = static_cast<std::size_t>(LR) + 1 + 64;
  std::vector<u32> hR, hS;
  std::vector<std::vector<int>> mons;
  for (;;) {
    mons = monomials(k, delta, need, rng, total);
    const int N = static_cast<int>(mons.size());
    // Series of each monomial in the jet, truncated at LR.
    std::vector<std::vector<u32>> cols(N, std::vector<u32>(LR, 0));
    std::vector<std::vector<u32>> Vlog(k, std::vector<u32>(LR));
    for (int a = 0; a < k; ++a) {
      for (int i = 0; i < LR; ++i) Vlog[a][i] = K.lg[J.V[a][i]];
    }
    std::vector<u32> acc(LR), nxt(LR);
    for (int c = 0; c < N; ++c) {
      const std::vector<int>& mu = mons[c];
      acc = J.V[mu[0]];
      acc.resize(LR);
      for (int d = 1; d < delta; ++d) {
        std::fill(nxt.begin(), nxt.end(), 0);
        for (int i = 0; i < LR; ++i) {
          if (acc[i]) K.axpy_log(nxt.data() + i, Vlog[mu[d]].data(), acc[i], 0, LR - i);
        }
        acc.swap(nxt);
      }
      cols[c] = acc;
    }
    // Rows [T^e] h(V) for e < LR, with two right-hand-side columns (R, S).
    const int len = N + 2;
    Echelon shared(K, len, N);
    std::vector<u32> row(len);
    for (int e = 0; e < LS; ++e) {
      for (int c = 0; c < N; ++c) row[c] = cols[c][e];
      row[N] = row[N + 1] = 0;
      shared.insert(row);
    }
    Echelon ES = shared, ER = shared;
    // S: [T^LS] S(V) = 1.
    for (int c = 0; c < N; ++c) row[c] = cols[c][LS];
    row[N] = 0;
    row[N + 1] = 1;
    ES.insert(row);
    // R: zero up to LR - 1, and R(y_j0) = 1.
    for (int e = LS; e < LR; ++e) {
      for (int c = 0; c < N; ++c) row[c] = cols[c][e];
      row[N] = row[N + 1] = 0;
      ER.insert(row);
    }
    for (int c = 0; c < N; ++c) {
      u32 v = 1;
      for (int a : mons[c]) v = K.mul(v, Y[a][j0]);
      row[c] = v;
    }
    row[N] = 1;
    row[N + 1] = 0;
    ER.insert(row);
    out.R_found = !ER.inconsistent;
    out.S_found = !ES.inconsistent;
    if ((out.R_found && out.S_found) || mons.size() >= total) {
      if (out.R_found) hR = ER.solve(N);
      if (out.S_found) hS = ES.solve(N + 1);
      break;
    }
    need = std::min<std::size_t>(total, need * 2);
  }
  out.monomials_total = static_cast<int>(total);
  out.monomials_used = static_cast<int>(mons.size());
  out.t_forms = now_s() - t0;
  if (!out.R_found || !out.S_found) {
    out.failure = !out.R_found ? "no R" : "no S";
    return out;
  }

  // Quotients at the public columns.
  const double t1 = now_s();
  const int N = static_cast<int>(mons.size());
  std::vector<u32> beta(n, 0);
  std::vector<char> inf(n, 0);
  const u64 sqrt_exp = K.size / 2;  // x^(2^(e-1)) is the square root in GF(2^e)
  for (int j = 0; j < n; ++j) {
    if (j == pos) continue;
    u32 r = 0, s = 0;
    for (int c = 0; c < N; ++c) {
      u32 v = 1;
      for (int a : mons[c]) v = K.mul(v, Y[a][j]);
      if (!v) continue;
      if (hR[c]) r = K.add(r, K.mul(hR[c], v));
      if (hS[c]) s = K.add(s, K.mul(hS[c], v));
    }
    if (!s) {
      inf[j] = 1;
      ++out.infinities;
      continue;
    }
    u32 b = K.div(r, s);
    if (binary) b = K.pow(b, sqrt_exp);
    beta[j] = b;
  }
  // Moebius move to a finite support when a point went to infinity.
  if (out.infinities) {
    std::vector<char> used(K.size, 0);
    for (int j = 0; j < n; ++j) {
      if (!inf[j]) used[beta[j]] = 1;
    }
    u32 p0 = 0;
    while (p0 < K.size && used[p0]) ++p0;
    if (p0 == K.size) { out.failure = "no free point for the Moebius move"; return out; }
    for (int j = 0; j < n; ++j) beta[j] = inf[j] ? 0 : K.inv(K.sub(beta[j], p0));
  }
  {
    std::vector<u32> sorted = beta;
    std::sort(sorted.begin(), sorted.end());
    out.distinct = std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();
  }
  out.alpha = beta;
  out.t_support = now_s() - t1;
  if (!out.distinct) { out.failure = "support values not distinct"; return out; }

  // Multiplier: c_a ⊙ z ∈ RS_{D+1}(alpha') for the basis codewords.
  const double t2 = now_s();
  const std::vector<u32> ones(n, 1);
  const std::vector<std::vector<u32>> Hrs = grs_parity(K, beta, ones, D);
  const int red = n - D - 1;
  int use = std::min(k, (n + 16 + red - 1) / red + 2);
  std::vector<u32> z;
  for (;;) {
    Echelon E(K, n, n);
    std::vector<u32> row(n);
    for (int a = 0; a < use; ++a) {
      for (int l = 0; l < red; ++l) {
        for (int j = 0; j < n; ++j) row[j] = K.mul(Hrs[l][j], Y[a][j]);
        E.insert(row);
      }
    }
    out.multiplier_nullity = n - E.rank;
    if (out.multiplier_nullity == 1 || use == k) {
      if (out.multiplier_nullity == 1) z = E.kernel_vector();
      break;
    }
    use = std::min(k, 2 * use);
  }
  out.t_multiplier = now_s() - t2;
  if (out.multiplier_nullity != 1) { out.failure = "multiplier not unique"; return out; }
  out.multiplier_nonzero = std::all_of(z.begin(), z.end(), [](u32 x) { return x != 0; });
  if (!out.multiplier_nonzero) { out.failure = "multiplier has zero entries"; return out; }
  out.lambda.resize(n);
  for (int j = 0; j < n; ++j) out.lambda[j] = K.inv(z[j]);

  // Verification against the public code only.
  const double t3 = now_s();
  const std::vector<std::vector<u32>> Hk = grs_parity(K, beta, out.lambda, D);
  out.contains = true;
  for (int a = 0; a < k && out.contains; ++a) {
    for (int l = 0; l < red && out.contains; ++l) {
      u32 sacc = 0;
      for (int j = 0; j < n; ++j) {
        if (Y[a][j]) sacc = K.add(sacc, K.mul(Hk[l][j], Y[a][j]));
      }
      if (sacc) out.contains = false;
    }
  }
  const int sub_dim = n - fq_rank(K, expand(K, Fq, Hk, n), n);
  out.equal = out.contains && sub_dim == k;
  out.t_verify = now_s() - t3;
  if (!out.equal) out.failure = out.contains ? "subfield subcode larger than C" : "C not contained";
  return out;
}

// alpha' = mu(phi^sigma(alpha)) for a Moebius map mu: all cross-ratios agree.
bool moebius_equivalent(const Field& K, const std::vector<u32>& a, const std::vector<u32>& b) {
  const int n = static_cast<int>(a.size());
  if (n < 4) return true;
  auto cr = [&](const std::vector<u32>& z, int j) {
    // (z0, z1; z2, zj) = (z2 - z0)(zj - z1) / ((z2 - z1)(zj - z0))
    const u32 num = K.mul(K.sub(z[2], z[0]), K.sub(z[j], z[1]));
    const u32 den = K.mul(K.sub(z[2], z[1]), K.sub(z[j], z[0]));
    return K.div(num, den);
  };
  for (int j = 3; j < n; ++j) {
    if (cr(a, j) != cr(b, j)) return false;
  }
  return true;
}

// --------------------------------------------------------------- driver

struct Options {
  Params P;
  u64 seed_lo = 1, seed_hi = 1;
  bool adjacent = false;  // binary Goppa control: adjacent orders instead of the square variant
  std::string label, out, dump;
};

u32 prime_of(u32 q, u32& s) {
  u32 p = 0;
  for (u32 d = 2; d <= q; ++d) {
    if (q % d == 0) { p = d; break; }
  }
  s = 0;
  u32 x = q;
  while (x % p == 0) { x /= p; ++s; }
  if (x != 1) die("q must be a prime power");
  return p;
}

int choose_delta(const Params& P, int k, int D) {
  if (P.family == Family::BinaryGoppa) return 2;
  if (P.delta) return P.delta;
  for (int d = 2; d <= 4; ++d) {
    if (binom(k + d - 1, d) >= static_cast<u64>(d * D + 2)) return d;
  }
  return 4;
}

int run_one(const Options& o, const Field& K, const Subfield& Fq, u64 seed, std::FILE* out) {
  const double t0 = now_s();
  Key key;
  std::string why;
  if (!make_key(K, Fq, o.P, seed, key, why)) die(why);
  std::vector<Poly> f;
  const bool curve_ok = curve_polys(K, key, f);
  if (!curve_ok) die("the public code does not lie on a curve of degree D");
  // The square-root variant is used for binary Goppa codes unless the
  // control with adjacent orders is requested.
  const bool binary = o.P.family == Family::BinaryGoppa && !o.adjacent;
  const int delta = choose_delta(o.P, key.k, key.D);
  const int len = delta * key.D;  // jet order len - 1
  Rng rng({static_cast<u32>(seed), static_cast<u32>(seed >> 32), static_cast<u32>(key.n),
           o.P.q, o.P.m, static_cast<u32>(o.P.family), 0x1E7u});
  const double t1 = now_s();
  const Jet J = make_jet(K, Fq, key, f, len, rng);
  bool jet_ok = true;
  for (int a = 0; a < key.k; ++a) jet_ok = jet_ok && J.V[a][0] == key.Y[a][J.pos];
  const double t2 = now_s();
  const Recovery R = direct_recover(K, Fq, key.Y, key.n, key.D, J, binary, delta, rng);
  const double t3 = now_s();
  // Audit with the secret key: the recovered support is a Moebius image of
  // the support on the branch of the jet.
  bool audit = false;
  if (R.equal) {
    std::vector<u32> conj(key.n);
    for (int j = 0; j < key.n; ++j) conj[j] = Fq.frob(K, key.alpha[j], J.sigma);
    audit = moebius_equivalent(K, conj, R.alpha);
  }
  const bool success = R.equal && audit;

  std::ostringstream js;
  js << "{\"label\":\"" << o.label << "\",\"family\":\"" << family_name(o.P.family)
     << "\",\"q\":" << o.P.q << ",\"m\":" << o.P.m << ",\"n\":" << key.n << ",\"t\":" << o.P.t
     << ",\"r\":" << o.P.r << ",\"k\":" << key.k << ",\"D\":" << key.D << ",\"seed\":" << seed
     << ",\"variant\":\"" << (binary ? "square-root" : "adjacent") << "\""
     << ",\"delta\":" << delta << ",\"jet_order\":" << len - 1 << ",\"position\":" << J.pos
     << ",\"branch\":" << J.sigma << ",\"monomials_total\":" << R.monomials_total
     << ",\"monomials_used\":" << R.monomials_used << ",\"R_found\":" << (R.R_found ? "true" : "false")
     << ",\"S_found\":" << (R.S_found ? "true" : "false") << ",\"infinities\":" << R.infinities
     << ",\"distinct\":" << (R.distinct ? "true" : "false")
     << ",\"multiplier_nullity\":" << R.multiplier_nullity
     << ",\"contains\":" << (R.contains ? "true" : "false")
     << ",\"equal\":" << (R.equal ? "true" : "false") << ",\"audit\":" << (audit ? "true" : "false")
     << ",\"success\":" << (success ? "true" : "false") << ",\"failure\":\"" << R.failure << "\""
     << ",\"checks\":{\"curve\":" << (curve_ok ? "true" : "false")
     << ",\"jet_starts_at_public_column\":" << (jet_ok ? "true" : "false") << "}"
     << ",\"seconds\":{\"key\":" << (t1 - t0) << ",\"jet\":" << (t2 - t1)
     << ",\"forms\":" << R.t_forms << ",\"support\":" << R.t_support
     << ",\"multiplier\":" << R.t_multiplier << ",\"verify\":" << R.t_verify
     << ",\"total\":" << (t3 - t0) << "}}";
  std::fprintf(out, "%s\n", js.str().c_str());
  std::fflush(out);
  if (!o.dump.empty() && R.equal) {
    // Public code, recovered key, and the secret data needed for the audit.
    const std::string path = o.dump + "/" + o.label + "-s" + std::to_string(seed) + ".txt";
    std::FILE* d = std::fopen(path.c_str(), "w");
    if (!d) die("cannot write " + path);
    std::fprintf(d, "p %u\ne %u\nq %u\nm %u\nn %d\nk %d\nD %d\nbranch %u\nmodulus", K.p, K.e,
                 o.P.q, o.P.m, key.n, key.k, key.D, J.sigma);
    for (u32 c : K.modulus) std::fprintf(d, " %u", c);
    auto line = [&](const char* name, const std::vector<u32>& v) {
      std::fprintf(d, "\n%s", name);
      for (u32 x : v) std::fprintf(d, " %u", x);
    };
    for (int a = 0; a < key.k; ++a) line("Y", key.Y[a]);
    line("alpha_recovered", R.alpha);
    line("lambda_recovered", R.lambda);
    line("alpha_secret", key.alpha);
    std::fprintf(d, "\n");
    std::fclose(d);
  }
  std::fprintf(stderr, "%s %s q=%u m=%u n=%d k=%d D=%d seed=%" PRIu64
               ": delta %d, jet order %d, branch %u, %s%s (%.2fs)\n",
               o.label.c_str(), family_name(o.P.family), o.P.q, o.P.m, key.n, key.k, key.D, seed,
               delta, len - 1, J.sigma, success ? "recovered" : "FAILED: ",
               success ? "" : (R.failure.empty() ? "audit" : R.failure.c_str()), t3 - t0);
  return (jet_ok && curve_ok) ? 0 : 1;
}

// ------------------------------------------------------------ self-test

int self_test() {
  int failures = 0;
  auto expect = [&](bool c, const std::string& what) {
    if (!c) {
      std::fprintf(stderr, "self-test FAILED: %s\n", what.c_str());
      ++failures;
    }
  };
  const std::pair<u32, u32> fields[] = {{2, 1}, {2, 4}, {2, 8}, {2, 12}, {2, 13}, {3, 1}, {3, 2},
                                        {3, 5}, {3, 7}, {5, 1}, {5, 3}, {5, 4}, {7, 3}, {11, 2},
                                        {13, 2}, {2, 16}};
  for (auto [p, e] : fields) {
    Field K(p, e);
    Rng rng({p, e, 1});
    for (int i = 0; i < 2000; ++i) {
      const u32 a = static_cast<u32>(rng.below(K.size)), b = static_cast<u32>(rng.below(K.size)),
                c = static_cast<u32>(rng.below(K.size));
      const std::string tag = " in GF(" + std::to_string(p) + "^" + std::to_string(e) + ")";
      expect(K.mul(a, b) == K.slow_mul(a, b), "multiplication" + tag);
      expect(K.add(a, b) == K.add(b, a), "commutative addition" + tag);
      expect(K.add(K.add(a, b), c) == K.add(a, K.add(b, c)), "associative addition" + tag);
      expect(K.mul(a, K.add(b, c)) == K.add(K.mul(a, b), K.mul(a, c)), "distributivity" + tag);
      expect(K.add(a, K.neg(a)) == 0, "negation" + tag);
      expect(K.pow(K.add(a, b), p) == K.add(K.pow(a, p), K.pow(b, p)), "Frobenius" + tag);
      if (a) expect(K.mul(a, K.inv(a)) == 1, "inverse" + tag);
    }
  }
  // Subfields, Taylor shift and interpolation.
  {
    Field K(3, 6);
    for (u32 s : {1u, 2u, 3u}) {
      Subfield Fq(K, s);
      for (u32 x : Fq.elems) expect(K.pow(x, Fq.q) == x, "subfield elements fixed by Frobenius");
    }
    Rng rng({5, 5});
    Poly a(9);
    for (u32& c : a) c = static_cast<u32>(rng.below(K.size));
    const u32 beta = static_cast<u32>(rng.below(K.size));
    const Poly sh = taylor_shift(K, a, beta);
    for (int i = 0; i < 20; ++i) {
      const u32 x = static_cast<u32>(rng.below(K.size));
      expect(peval(K, sh, x) == peval(K, a, K.add(beta, x)), "Taylor shift");
    }
    std::vector<u32> xs(9), ys(9);
    for (int i = 0; i < 9; ++i) { xs[i] = static_cast<u32>(3 * i + 1); ys[i] = peval(K, a, xs[i]); }
    Poly b = Interpolator(K, xs)(ys);
    b.resize(9, 0);
    expect(b == a, "interpolation");
  }
  // Number of monic irreducible polynomials of degree t over GF(4) and GF(9).
  for (auto [p, e] : std::vector<std::pair<u32, u32>>{{2, 2}, {3, 2}}) {
    Field K(p, e);
    for (int t = 2; t <= 3; ++t) {
      const u64 q = K.size;
      const u64 expected = t == 2 ? (q * q - q) / 2 : (q * q * q - q) / 3;
      u64 count = 0, total = 1;
      for (int i = 0; i < t; ++i) total *= q;
      for (u64 idx = 0; idx < total; ++idx) {
        Poly g(t + 1);
        u64 r = idx;
        for (int i = 0; i < t; ++i) { g[i] = static_cast<u32>(r % q); r /= q; }
        g[t] = 1;
        count += irreducible(K, g);
      }
      expect(count == expected, "irreducible count over GF(" + std::to_string(q) + ")");
    }
  }
  // End-to-end recovery on small keys of each family.
  struct Case { Family f; u32 q, m; int n, t, r; };
  const Case cases[] = {{Family::BinaryGoppa, 2, 8, 100, 8, 0}, {Family::WildGoppa, 3, 5, 120, 4, 0},
                        {Family::WildGoppa, 4, 4, 120, 4, 0},   {Family::Alternant, 2, 8, 120, 0, 10},
                        {Family::Alternant, 3, 5, 120, 0, 8},   {Family::Alternant, 5, 3, 100, 0, 6}};
  int recovered = 0, total = 0;
  for (const Case& c : cases) {
    u32 s;
    const u32 p = prime_of(c.q, s);
    Field K(p, s * c.m);
    Subfield Fq(K, s);
    for (u64 seed = 1; seed <= 4; ++seed) {
      Params P;
      P.family = c.f;
      P.q = c.q;
      P.m = c.m;
      P.n = c.n;
      P.t = c.t;
      P.r = c.r;
      Key key;
      std::string why;
      if (!make_key(K, Fq, P, seed, key, why)) die(why);
      std::vector<Poly> f;
      expect(curve_polys(K, key, f), "curve of degree D");
      Rng rng({static_cast<u32>(seed), 3});
      const int delta = choose_delta(P, key.k, key.D);
      const Jet J = make_jet(K, Fq, key, f, delta * key.D, rng);
      for (int a = 0; a < key.k; ++a) expect(J.V[a][0] == key.Y[a][J.pos], "jet starts at y_i");
      const Recovery R = direct_recover(K, Fq, key.Y, key.n, key.D, J, c.f == Family::BinaryGoppa,
                                        delta, rng);
      std::vector<u32> conj(key.n);
      for (int j = 0; j < key.n; ++j) conj[j] = Fq.frob(K, key.alpha[j], J.sigma);
      ++total;
      if (R.equal && moebius_equivalent(K, conj, R.alpha)) ++recovered;
      else std::fprintf(stderr, "self-test: %s seed %" PRIu64 " not recovered (%s)\n",
                        family_name(c.f), seed, R.failure.c_str());
      // A support that is not a Moebius image must fail the audit.
      if (R.equal && key.n > 4) {
        std::vector<u32> bad = R.alpha;
        std::swap(bad[1], bad[3]);
        expect(!moebius_equivalent(K, conj, bad), "audit rejects a permuted support");
      }
    }
  }
  expect(recovered == total, "all small keys recovered");
  std::fprintf(stderr, "self-test: %d of %d small keys recovered, %s\n", recovered, total,
               failures ? "FAILED" : "ok");
  return failures ? 1 : 0;
}

void parse_range(const std::string& s, long& lo, long& hi) {
  const auto dash = s.find('-');
  if (dash == std::string::npos) {
    lo = hi = std::stol(s);
  } else {
    lo = std::stol(s.substr(0, dash));
    hi = std::stol(s.substr(dash + 1));
  }
  if (lo > hi) die("bad range " + s);
}

void usage() {
  std::fprintf(stderr,
               "usage:\n"
               "  synthetic_direct self-test\n"
               "  synthetic_direct run --family F --q Q --m M --n N (--t T | --r R) [options]\n"
               "families: binary-goppa (q = 2), wild-goppa (q > 2), alternant\n"
               "options:\n"
               "  --seeds A[-B]   seeds (default 1)\n"
               "  --delta D       degree of the forms (default: smallest admissible)\n"
               "  --adjacent      binary Goppa control: adjacent orders, no square root\n"
               "  --dump DIR      write recovered keys for crosscheck.py\n"
               "  --label L  --out FILE\n");
  std::exit(2);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) usage();
  const std::string cmd = argv[1];
  if (cmd == "self-test") return self_test();
  if (cmd != "run") usage();
  Options o;
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    auto val = [&]() -> std::string {
      if (i + 1 >= argc) die("missing value for " + a);
      return argv[++i];
    };
    long lo, hi;
    if (a == "--family") {
      const std::string v = val();
      if (v == "binary-goppa") o.P.family = Family::BinaryGoppa;
      else if (v == "wild-goppa") o.P.family = Family::WildGoppa;
      else if (v == "alternant") o.P.family = Family::Alternant;
      else die("unknown family " + v);
    } else if (a == "--q") o.P.q = static_cast<u32>(std::stoul(val()));
    else if (a == "--m") o.P.m = static_cast<u32>(std::stoul(val()));
    else if (a == "--n") o.P.n = std::stoi(val());
    else if (a == "--t") o.P.t = std::stoi(val());
    else if (a == "--r") o.P.r = std::stoi(val());
    else if (a == "--delta") o.P.delta = std::stoi(val());
    else if (a == "--seeds") {
      parse_range(val(), lo, hi);
      o.seed_lo = static_cast<u64>(lo);
      o.seed_hi = static_cast<u64>(hi);
    } else if (a == "--label") o.label = val();
    else if (a == "--out") o.out = val();
    else if (a == "--dump") o.dump = val();
    else if (a == "--adjacent") o.adjacent = true;
    else die("unknown option " + a);
  }
  if (o.P.n == 0) usage();
  u32 s;
  const u32 p = prime_of(o.P.q, s);
  if (o.label.empty()) o.label = family_name(o.P.family);
  Field K(p, s * o.P.m);
  Subfield Fq(K, s);
  std::FILE* out = stdout;
  if (!o.out.empty()) {
    out = std::fopen(o.out.c_str(), "a");
    if (!out) die("cannot open " + o.out);
  }
  int bad = 0;
  for (u64 seed = o.seed_lo; seed <= o.seed_hi; ++seed) bad += run_one(o, K, Fq, seed, out);
  if (out != stdout) std::fclose(out);
  return bad ? 1 : 0;
}
