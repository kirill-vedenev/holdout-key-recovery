// Ground-truth check of the existence condition for binary direct recovery.
//
// For a binary square-free Goppa code Gamma(alpha, G) of length n with
// deg G = t, the public code lies in GRS_{D+1}(alpha, lambda) with
// D = n - 2t - 1. Write f_1, ..., f_k for the polynomials of degree at most D
// that correspond to an F_2-basis of the code, and let
//
//     W = span_K { f_a f_b : 1 <= a <= b <= k },   K = GF(2^m).
//
// Every product satisfies G^2 | (f_a f_b)', so dim W <= 2D - t + 1. If
// equality holds, then (Z - alpha_i)^{2D} and (Z - alpha_i)^{2D-2} lie in W
// for every position i, and the quadratic forms R and S of the binary direct
// recovery exist. This program generates keys, computes dim W from the
// secret curve, and stops as soon as the bound 2D - t + 1 is reached.
//
// Products are represented by their values at N = 2D + 1 distinct points of
// K: the n support points and N - n further points. Evaluation is injective
// on polynomials of degree at most 2D, so rank is preserved.

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using u16 = std::uint16_t;
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

// ------------------------------------------------------------------ RNG

class Rng {
 public:
  explicit Rng(const std::vector<u32>& s) {
    std::seed_seq q(s.begin(), s.end());
    g_.seed(q);
  }
  u64 next() { return g_(); }
  u64 below(u64 b) {  // unbiased integer in [0, b)
    const u64 threshold = (0 - b) % b;
    for (;;) {
      const u64 x = g_();
      if (x >= threshold) return x % b;
    }
  }
  template <class T>
  void shuffle(std::vector<T>& v) {  // explicit Fisher-Yates
    for (std::size_t i = v.size(); i > 1; --i) {
      std::swap(v[i - 1], v[below(i)]);
    }
  }

 private:
  std::mt19937_64 g_;
};

// ------------------------------------------------------------ GF(2^m)

u32 default_modulus(int m) {
  switch (m) {
    case 2: return 0x7;
    case 3: return 0xB;
    case 4: return 0x13;
    case 5: return 0x25;
    case 6: return 0x43;      // x^6 + x + 1
    case 7: return 0x83;      // x^7 + x + 1
    case 8: return 0x11D;     // x^8 + x^4 + x^3 + x^2 + 1
    case 9: return 0x211;     // x^9 + x^4 + 1
    case 10: return 0x409;    // x^10 + x^3 + 1
    case 11: return 0x805;    // x^11 + x^2 + 1
    case 12: return 0x1009;   // x^12 + x^3 + 1 (Classic McEliece)
    case 13: return 0x201B;   // x^13 + x^4 + x^3 + x + 1 (Classic McEliece)
    case 14: return 0x4443;   // x^14 + x^10 + x^6 + x + 1
    case 15: return 0x8003;   // x^15 + x + 1
    default: die("m must be in [2, 15]");
  }
}

struct GF {
  int m = 0;
  u32 size = 0, Q = 0, modulus = 0, gen = 0, ZL = 0;
  std::vector<u16> ex;  // ex[i] = gen^(i mod Q) for i < 2Q, zero for i >= 2Q
  std::vector<u16> lg;  // lg[x] = discrete log, lg[0] = ZL = 2Q

  static u32 slow_mul(u32 a, u32 b, int m, u32 mod) {
    u32 r = 0;
    while (b) {
      if (b & 1) r ^= a;
      b >>= 1;
      a <<= 1;
      if (a >> m & 1) a ^= mod;
    }
    return r;
  }

  GF(int m_, u32 mod_) : m(m_), modulus(mod_) {
    if (m < 2 || m > 15) die("m must be in [2, 15]");
    if ((modulus >> m) != 1) die("modulus must have degree m");
    size = 1u << m;
    Q = size - 1;
    ZL = 2 * Q;
    std::vector<u32> primes;
    u32 x = Q;
    for (u32 p = 2; p * p <= x; ++p) {
      if (x % p == 0) {
        primes.push_back(p);
        while (x % p == 0) x /= p;
      }
    }
    if (x > 1) primes.push_back(x);
    auto spow = [&](u32 a, u64 e) {
      u32 r = 1;
      while (e) {
        if (e & 1) r = slow_mul(r, a, m, modulus);
        a = slow_mul(a, a, m, modulus);
        e >>= 1;
      }
      return r;
    };
    for (u32 c = 2; c < size && gen == 0; ++c) {
      if (spow(c, Q) != 1) continue;
      bool primitive = true;
      for (u32 p : primes) {
        if (spow(c, Q / p) == 1) primitive = false;
      }
      if (primitive) gen = c;
    }
    // An element of order 2^m - 1 exists only if the modulus is irreducible.
    if (gen == 0) die("modulus is not irreducible");
    ex.assign(4 * Q + 1, 0);
    lg.assign(size, static_cast<u16>(ZL));
    u32 v = 1;
    for (u32 i = 0; i < Q; ++i) {
      if (lg[v] != ZL) die("generator has small order");
      ex[i] = static_cast<u16>(v);
      ex[i + Q] = static_cast<u16>(v);
      lg[v] = static_cast<u16>(i);
      v = slow_mul(v, gen, m, modulus);
    }
    if (v != 1) die("generator order mismatch");
  }

  u16 mul(u16 a, u16 b) const { return ex[lg[a] + lg[b]]; }
  u16 inv(u16 a) const {
    if (a == 0) die("inverse of zero");
    return ex[Q - lg[a]];
  }
  u16 div(u16 a, u16 b) const { return mul(a, inv(b)); }
  u16 pow(u16 a, u64 e) const {
    if (e == 0) return 1;
    if (a == 0) return 0;
    return ex[(static_cast<u64>(lg[a]) * (e % Q)) % Q];
  }
};

// --------------------------------------------------- polynomials over K

using Poly = std::vector<u16>;  // coefficients from degree 0 upward

void trim(Poly& p) {
  while (!p.empty() && p.back() == 0) p.pop_back();
}
int deg(const Poly& p) { return static_cast<int>(p.size()) - 1; }

Poly pmul(const GF& F, const Poly& a, const Poly& b) {
  if (a.empty() || b.empty()) return {};
  Poly r(a.size() + b.size() - 1, 0);
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (!a[i]) continue;
    for (std::size_t j = 0; j < b.size(); ++j) r[i + j] ^= F.mul(a[i], b[j]);
  }
  trim(r);
  return r;
}

// a <- a mod g, for a monic g of degree >= 1.
void pmod(const GF& F, Poly& a, const Poly& g) {
  const int t = deg(g);
  for (int d = deg(a); d >= t; --d) {
    const u16 c = a[d];
    if (!c) continue;
    const u32 lc = F.lg[c];
    for (int j = 0; j <= t; ++j) a[d - t + j] ^= F.ex[F.lg[g[j]] + lc];
  }
  if (deg(a) >= t) a.resize(t);
  trim(a);
}

Poly monic(const GF& F, Poly p) {
  trim(p);
  if (p.empty()) return p;
  const u16 il = F.inv(p.back());
  for (auto& c : p) c = F.mul(c, il);
  return p;
}

Poly pgcd(const GF& F, Poly a, Poly b) {
  trim(a);
  trim(b);
  while (!b.empty()) {
    b = monic(F, b);
    pmod(F, a, b);
    std::swap(a, b);
  }
  return monic(F, a);
}

u16 peval(const GF& F, const Poly& p, u16 x) {
  u16 r = 0;
  for (int i = deg(p); i >= 0; --i) r = static_cast<u16>(F.mul(r, x) ^ p[i]);
  return r;
}

Poly pderiv(const Poly& p) {  // characteristic two
  Poly r(p.size() > 1 ? p.size() - 1 : 0, 0);
  for (std::size_t i = 1; i < p.size(); i += 2) r[i - 1] = p[i];
  trim(r);
  return r;
}

// Ben-Or test for a monic polynomial over K = GF(q), q = 2^m.
bool irreducible(const GF& F, const Poly& g) {
  const int t = deg(g);
  if (t <= 0) return false;
  if (t == 1) return true;
  Poly h = {0, 1};
  for (int i = 1; i <= t / 2; ++i) {
    for (int s = 0; s < F.m; ++s) {  // h <- h^2 mod g, m times: h = x^{q^i}
      Poly sq(2 * h.size(), 0);
      for (std::size_t j = 0; j < h.size(); ++j) sq[2 * j] = F.mul(h[j], h[j]);
      trim(sq);
      pmod(F, sq, g);
      h = sq;
    }
    Poly d = h;
    if (d.size() < 2) d.resize(2, 0);
    d[1] ^= 1;
    trim(d);
    if (deg(pgcd(F, g, d)) > 0) return false;
  }
  return true;
}

bool squarefree(const GF& F, const Poly& g) {
  return deg(pgcd(F, g, pderiv(g))) == 0;
}

// Newton interpolation through (xs[i], ys[i]); returns monomial coefficients.
Poly interpolate(const GF& F, const std::vector<u16>& xs, std::vector<u16> c) {
  const int n = static_cast<int>(xs.size());
  for (int j = 1; j < n; ++j) {
    for (int i = n - 1; i >= j; --i) {
      c[i] = F.div(static_cast<u16>(c[i] ^ c[i - 1]),
                   static_cast<u16>(xs[i] ^ xs[i - j]));
    }
  }
  Poly p = {c[n - 1]};
  for (int i = n - 2; i >= 0; --i) {
    Poly q(p.size() + 1, 0);  // q = p * (x - xs[i]) + c[i]
    for (std::size_t j = 0; j < p.size(); ++j) {
      q[j + 1] ^= p[j];
      q[j] ^= F.mul(p[j], xs[i]);
    }
    q[0] ^= c[i];
    p = q;
  }
  trim(p);
  return p;
}

// ------------------------------------------------------------- instance

enum class GoppaKind { Irreducible, SquareFree };

const char* kind_name(GoppaKind k) {
  return k == GoppaKind::Irreducible ? "irreducible" : "squarefree";
}

struct Instance {
  int m = 0, t = 0, n = 0, k = 0, D = 0, N = 0;
  u64 seed = 0;
  GoppaKind kind = GoppaKind::Irreducible;
  Poly G;
  bool coeff = false;     // products in coefficient form instead of values
  std::vector<u16> pts;   // pts[0..n) = support, pts[n..N) = further points
                          // (only the support in coefficient form)
  int words = 0;          // 64-bit words per codeword
  std::vector<u64> code;  // k codewords: an F_2-basis of Gamma(alpha, G)

  bool bit(int a, int j) const {
    return code[static_cast<std::size_t>(a) * words + j / 64] >> (j % 64) & 1;
  }
};

// Returns false (with a reason) when the parameters admit no instance.
// A nonzero `stream` (the full Classic McEliece length) separates parameter
// sets whose shortened codes have the same (m, t, n).
//
// Products are represented by their values at 2D + 1 points when the field
// has enough points, and by their coefficients otherwise (or when forced).
bool build_instance(const GF& F, int t, int n, u64 seed, GoppaKind kind, u32 stream,
                    bool force_coeff, Instance& I, std::string& why) {
  I = Instance();
  I.m = F.m;
  I.t = t;
  I.n = n;
  I.seed = seed;
  I.kind = kind;
  I.D = n - 2 * t - 1;
  I.N = 2 * I.D + 1;
  if (t < 1) { why = "t must be positive"; return false; }
  if (I.D < t) { why = "need D = n - 2t - 1 >= t"; return false; }
  I.coeff = force_coeff || I.N < n || static_cast<u64>(I.N) > F.size;
  const int npts = I.coeff ? n : I.N;

  std::vector<u32> seed_words = {static_cast<u32>(F.m), static_cast<u32>(t), static_cast<u32>(n),
                                 static_cast<u32>(seed), static_cast<u32>(seed >> 32),
                                 static_cast<u32>(kind), 0x6077A};
  if (stream) seed_words.push_back(stream);
  Rng rng(seed_words);
  for (int tries = 0;; ++tries) {
    if (tries > 1000000) { why = "no Goppa polynomial found"; return false; }
    Poly g(t + 1);
    for (int i = 0; i < t; ++i) g[i] = static_cast<u16>(rng.next() >> (64 - F.m));
    g[t] = 1;
    const bool ok = kind == GoppaKind::Irreducible ? irreducible(F, g) : squarefree(F, g);
    if (ok) { I.G = g; break; }
  }
  std::vector<u16> cand;
  for (u32 x = 0; x < F.size; ++x) {
    if (peval(F, I.G, static_cast<u16>(x)) != 0) cand.push_back(static_cast<u16>(x));
  }
  if (static_cast<int>(cand.size()) < npts) { why = "too few non-roots of G"; return false; }
  rng.shuffle(cand);
  I.pts.assign(cand.begin(), cand.begin() + npts);

  // Binary parity-check matrix: bit b of alpha_j^i / G(alpha_j), 0 <= i < t.
  const int rows = F.m * t;
  const int W = (n + 63) / 64;
  std::vector<u64> H(static_cast<std::size_t>(rows) * W, 0);
  for (int j = 0; j < n; ++j) {
    u16 e = F.inv(peval(F, I.G, I.pts[j]));
    for (int i = 0; i < t; ++i) {
      for (int b = 0; b < F.m; ++b) {
        if (e >> b & 1) H[static_cast<std::size_t>(i * F.m + b) * W + j / 64] |= u64{1} << (j % 64);
      }
      e = F.mul(e, I.pts[j]);
    }
  }
  std::vector<int> pivcol;
  int r = 0;
  for (int c = 0; c < n && r < rows; ++c) {
    int p = -1;
    for (int i = r; i < rows; ++i) {
      if (H[static_cast<std::size_t>(i) * W + c / 64] >> (c % 64) & 1) { p = i; break; }
    }
    if (p < 0) continue;
    if (p != r) {
      for (int w = 0; w < W; ++w) std::swap(H[static_cast<std::size_t>(p) * W + w], H[static_cast<std::size_t>(r) * W + w]);
    }
    for (int i = 0; i < rows; ++i) {
      if (i != r && (H[static_cast<std::size_t>(i) * W + c / 64] >> (c % 64) & 1)) {
        for (int w = 0; w < W; ++w) H[static_cast<std::size_t>(i) * W + w] ^= H[static_cast<std::size_t>(r) * W + w];
      }
    }
    pivcol.push_back(c);
    ++r;
  }
  std::vector<char> is_piv(n, 0);
  for (int c : pivcol) is_piv[c] = 1;
  I.words = W;
  for (int f = 0; f < n; ++f) {
    if (is_piv[f]) continue;
    std::vector<u64> v(W, 0);
    v[f / 64] |= u64{1} << (f % 64);
    for (int i = 0; i < r; ++i) {
      if (H[static_cast<std::size_t>(i) * W + f / 64] >> (f % 64) & 1) {
        v[pivcol[i] / 64] |= u64{1} << (pivcol[i] % 64);
      }
    }
    I.code.insert(I.code.end(), v.begin(), v.end());
  }
  I.k = n - r;
  if (I.k < 1) { why = "the Goppa code is zero"; return false; }
  return true;
}

// Checks Sum_j c_j / (Z - alpha_j) == 0 mod G^2 for every basis codeword,
// i.e. Gamma(alpha, G) = Gamma(alpha, G^2). This is what makes f_c / G^2 a
// polynomial of degree at most D, used for the values at the further points.
bool check_goppa_square(const GF& F, const Instance& I) {
  const Poly G2 = pmul(F, I.G, I.G);
  const int L = 2 * I.t;
  std::vector<u16> inv_lin(static_cast<std::size_t>(I.n) * L);
  for (int j = 0; j < I.n; ++j) {
    // (x - a)^{-1} = ((G2(x) - G2(a)) / (x - a)) / G2(a)   mod G2
    const u16 a = I.pts[j];
    std::vector<u16> q(L, 0);
    u16 acc = G2[L];
    for (int d = L - 1; d >= 0; --d) {
      q[d] = acc;
      acc = static_cast<u16>(F.mul(acc, a) ^ G2[d]);
    }
    const u16 s = F.inv(acc);  // acc = G2(a)
    for (int d = 0; d < L; ++d) inv_lin[static_cast<std::size_t>(j) * L + d] = F.mul(q[d], s);
  }
  std::vector<u16> sum(L);
  for (int a = 0; a < I.k; ++a) {
    std::fill(sum.begin(), sum.end(), 0);
    for (int j = 0; j < I.n; ++j) {
      if (!I.bit(a, j)) continue;
      const u16* p = &inv_lin[static_cast<std::size_t>(j) * L];
      for (int d = 0; d < L; ++d) sum[d] ^= p[d];
    }
    for (int d = 0; d < L; ++d) {
      if (sum[d]) return false;
    }
  }
  return true;
}

// Log-domain values flog[a * N + x] = log f_a(pts[x]) of the secret curve.
// The stride is the number of points, pts.size().
std::vector<u16> curve_values(const GF& F, const Instance& I, int threads) {
  const int n = I.n, N = static_cast<int>(I.pts.size()), E = N - n, k = I.k;
  const u32 Q = F.Q;
  // sigma(x) = Pi'(alpha_x) / G(alpha_x)^2 on the support,
  // sigma(b) = Pi(b) / G(b)^2 at the further points.
  std::vector<u32> sig(N);
  for (int x = 0; x < N; ++x) {
    u64 s = 0;
    for (int l = 0; l < n; ++l) {
      if (l != x) s += F.lg[I.pts[x] ^ I.pts[l]];
    }
    s += 2 * static_cast<u64>(Q - F.lg[peval(F, I.G, I.pts[x])]);
    sig[x] = static_cast<u32>(s % Q);
  }
  // invT[j * E + e] = 1 / (b_e - alpha_j)
  std::vector<u16> invT(static_cast<std::size_t>(n) * E);
  for (int j = 0; j < n; ++j) {
    for (int e = 0; e < E; ++e) invT[static_cast<std::size_t>(j) * E + e] = F.inv(I.pts[n + e] ^ I.pts[j]);
  }
  std::vector<u16> flog(static_cast<std::size_t>(k) * N);
  auto work = [&](int tid) {
    std::vector<u16> acc(E);
    for (int a = tid; a < k; a += threads) {
      u16* row = &flog[static_cast<std::size_t>(a) * N];
      std::fill(acc.begin(), acc.end(), 0);
      for (int j = 0; j < n; ++j) {
        if (I.bit(a, j)) {
          row[j] = static_cast<u16>(sig[j]);
          const u16* p = &invT[static_cast<std::size_t>(j) * E];
          for (int e = 0; e < E; ++e) acc[e] ^= p[e];
        } else {
          row[j] = static_cast<u16>(F.ZL);
        }
      }
      for (int e = 0; e < E; ++e) {
        row[n + e] = acc[e] ? static_cast<u16>((F.lg[acc[e]] + sig[n + e]) % Q) : static_cast<u16>(F.ZL);
      }
    }
  };
  std::vector<std::thread> pool;
  for (int i = 0; i < threads; ++i) pool.emplace_back(work, i);
  for (auto& th : pool) th.join();
  return flog;
}


// Newton interpolation of many value vectors at the same nodes.
class Interpolator {
 public:
  Interpolator(const GF& F, std::vector<u16> xs) : F_(&F), xs_(std::move(xs)), n_(static_cast<int>(xs_.size())) {
    inv_.assign(static_cast<std::size_t>(n_) * n_, 0);
    for (int j = 1; j < n_; ++j) {
      for (int i = j; i < n_; ++i) inv_[static_cast<std::size_t>(j) * n_ + i] = F.lg[F.inv(xs_[i] ^ xs_[i - j])];
    }
  }
  Poly operator()(std::vector<u16> c) const {
    const GF& F = *F_;
    for (int j = 1; j < n_; ++j) {
      const u16* inv = &inv_[static_cast<std::size_t>(j) * n_];
      for (int i = n_ - 1; i >= j; --i) c[i] = F.ex[F.lg[c[i] ^ c[i - 1]] + inv[i]];
    }
    Poly p(n_, 0);  // Horner from the top Newton coefficient
    int len = 1;
    p[0] = c[n_ - 1];
    for (int i = n_ - 2; i >= 0; --i) {
      for (int d = len; d >= 1; --d) p[d] = static_cast<u16>(p[d - 1] ^ F.mul(p[d], xs_[i]));
      p[0] = static_cast<u16>(F.mul(p[0], xs_[i]) ^ c[i]);
      ++len;
    }
    trim(p);
    return p;
  }

 private:
  const GF* F_;
  std::vector<u16> xs_;
  int n_;
  std::vector<u16> inv_;
};

// The secret curve in the representation used for the products.
struct Curve {
  bool coeff = false;
  int N = 0, D = 0;
  int stride = 0;          // number of points in `vals`
  std::vector<u16> vals;   // log f_a(pts[x]), x < stride
  std::vector<u16> coef;   // coefficient form: log of the coefficients of f_a, degrees 0..D
  bool degree_ok = true;   // coefficient form: the interpolant matches every support value
};

Curve build_curve(const GF& F, const Instance& I, int threads) {
  Curve C;
  C.coeff = I.coeff;
  C.N = I.N;
  C.D = I.D;
  C.stride = static_cast<int>(I.pts.size());
  C.vals = curve_values(F, I, threads);
  if (!C.coeff) return C;
  // f_a has degree at most D: interpolate at the first D + 1 support points and
  // check the values at the remaining 2t support points.
  const int L = I.D + 1;
  Interpolator interp(F, std::vector<u16>(I.pts.begin(), I.pts.begin() + L));
  C.coef.assign(static_cast<std::size_t>(I.k) * L, static_cast<u16>(F.ZL));
  std::vector<char> ok(I.k, 1);
  auto work = [&](int tid) {
    for (int a = tid; a < I.k; a += threads) {
      const u16* va = &C.vals[static_cast<std::size_t>(a) * C.stride];
      std::vector<u16> ys(L);
      for (int x = 0; x < L; ++x) ys[x] = F.ex[va[x]];
      const Poly f = interp(ys);
      if (deg(f) > I.D) ok[a] = 0;
      for (int j = L; j < I.n; ++j) {
        if (peval(F, f, I.pts[j]) != F.ex[va[j]]) ok[a] = 0;
      }
      for (int d = 0; d < L && d < static_cast<int>(f.size()); ++d) {
        C.coef[static_cast<std::size_t>(a) * L + d] = F.lg[f[d]];
      }
    }
  };
  std::vector<std::thread> pool;
  for (int i = 0; i < threads; ++i) pool.emplace_back(work, i);
  for (auto& th : pool) th.join();
  C.degree_ok = std::all_of(ok.begin(), ok.end(), [](char c) { return c != 0; });
  return C;
}

// Writes the product f_a f_b (values or coefficients) to v[0..N).
void make_row(const GF& F, const Curve& C, int a, int b, u16* v) {
  const u16* ex = F.ex.data();
  if (!C.coeff) {
    const u16* fa = &C.vals[static_cast<std::size_t>(a) * C.stride];
    const u16* fb = &C.vals[static_cast<std::size_t>(b) * C.stride];
    for (int x = 0; x < C.N; ++x) v[x] = ex[fa[x] + fb[x]];
    return;
  }
  const int L = C.D + 1;
  std::fill(v, v + C.N, 0);
  const u16* ca = &C.coef[static_cast<std::size_t>(a) * L];
  const u16* cb = &C.coef[static_cast<std::size_t>(b) * L];
  for (int i = 0; i < L; ++i) {
    const u32 la = ca[i];
    if (la == F.ZL) continue;
    u16* w = v + i;
    for (int j = 0; j < L; ++j) w[j] ^= ex[la + cb[j]];
  }
}

// ----------------------------------------------------------- rank engine

// Semi-echelon basis. Row r is stored in log form, is zero before its pivot
// column col[r], has value 1 at col[r], and is zero at col[s] for s < r.
struct Basis {
  const GF* F;
  int N;
  int cap;
  std::vector<u16> rows;
  std::vector<int> col;
  int rank = 0;

  Basis(const GF& f, int n, int c)
      : F(&f), N(n), cap(c), rows(static_cast<std::size_t>(c) * n) {}

  void reduce(u16* v, int r0, int r1) const {
    const u16* ex = F->ex.data();
    const u16* lg = F->lg.data();
    for (int r = r0; r < r1; ++r) {
      const int c = col[r];
      const u16 x = v[c];
      if (!x) continue;
      const u32 lc = lg[x];
      const u16* p = rows.data() + static_cast<std::size_t>(r) * N;
      for (int i = c; i < N; ++i) v[i] ^= ex[p[i] + lc];
    }
  }

  bool insert(const u16* v) {
    int c = 0;
    while (c < N && !v[c]) ++c;
    if (c == N) return false;
    if (rank == cap) die("basis capacity exceeded");
    const u32 li = F->Q - F->lg[v[c]];
    u16* p = rows.data() + static_cast<std::size_t>(rank) * N;
    for (int i = 0; i < N; ++i) {
      p[i] = v[i] ? static_cast<u16>((F->lg[v[i]] + li) % F->Q) : static_cast<u16>(F->ZL);
    }
    col.push_back(c);
    ++rank;
    return true;
  }

  bool contains(std::vector<u16> v) const {
    reduce(v.data(), 0, rank);
    for (u16 x : v) {
      if (x) return false;
    }
    return true;
  }
};

using Pair = std::pair<u16, u16>;

struct RankRun {
  long rows_used = 0;
  bool stopped_at_bound = false;
};

// Adds product rows in the given order; stops at `bound` unless `full`.
RankRun product_rank(const GF& F, const Curve& C, const std::vector<Pair>& pairs,
                     int bound, bool full, int threads, Basis& B) {
  RankRun out;
  const int N = C.N;
  const long P = static_cast<long>(pairs.size());
  const int bsz = std::max(16, 8 * threads);
  std::vector<u16> buf(static_cast<std::size_t>(bsz) * N);
  long next = 0;
  while (next < P && (full || B.rank < bound)) {
    const int cnt = static_cast<int>(std::min<long>(bsz, P - next));
    const int r0 = B.rank;
    auto work = [&](int tid) {
      for (int i = tid; i < cnt; i += threads) {
        u16* v = &buf[static_cast<std::size_t>(i) * N];
        const Pair pr = pairs[next + i];
        make_row(F, C, pr.first, pr.second, v);
        B.reduce(v, 0, r0);
      }
    };
    if (threads == 1) {
      work(0);
    } else {
      std::vector<std::thread> pool;
      for (int i = 0; i < threads; ++i) pool.emplace_back(work, i);
      for (auto& th : pool) th.join();
    }
    for (int i = 0; i < cnt; ++i) {
      u16* v = &buf[static_cast<std::size_t>(i) * N];
      B.reduce(v, r0, B.rank);
      B.insert(v);
      out.rows_used = next + i + 1;
      if (!full && B.rank >= bound) {
        out.stopped_at_bound = true;
        break;
      }
    }
    next += cnt;
  }
  if (out.rows_used == 0) out.rows_used = next;
  return out;
}

// Rank over F_2 of the binary Schur square {c_a * c_b}, stopping at n.
int square_code_rank(const Instance& I, const std::vector<Pair>& pairs) {
  const int W = I.words;
  std::vector<u64> piv;
  std::vector<int> pc;
  std::vector<u64> v(W);
  int rank = 0;
  for (const Pair& pr : pairs) {
    const u64* ca = &I.code[static_cast<std::size_t>(pr.first) * W];
    const u64* cb = &I.code[static_cast<std::size_t>(pr.second) * W];
    for (int w = 0; w < W; ++w) v[w] = ca[w] & cb[w];
    for (int r = 0; r < rank; ++r) {
      if (v[pc[r] / 64] >> (pc[r] % 64) & 1) {
        const u64* p = &piv[static_cast<std::size_t>(r) * W];
        for (int w = 0; w < W; ++w) v[w] ^= p[w];
      }
    }
    int c = -1;
    for (int w = 0; w < W; ++w) {
      if (v[w]) { c = 64 * w + __builtin_ctzll(v[w]); break; }
    }
    if (c < 0) continue;
    piv.insert(piv.end(), v.begin(), v.end());
    pc.push_back(c);
    if (++rank == I.n) break;
  }
  return rank;
}

// Columns of the generator matrix whose rows are the basis codewords, as
// points of the curve: column j is F(alpha_j) up to scaling. The leading
// coefficient of f_a is the weight of c_a mod 2, so F(infinity) is the parity
// column, and infinity behaves like one more support point.
struct ColumnStats {
  int zero_columns = 0;          // finite base points: c_j = 0
  int duplicate_excess = 0;      // finite double points: sum of (class size - 1)
  bool parity_zero = false;      // base point at infinity: all weights even
  int duplicate_excess_ext = 0;  // double points among the n + 1 columns
};

ColumnStats column_stats(const Instance& I) {
  ColumnStats s;
  const int kw = (I.k + 63) / 64;
  std::map<std::vector<u64>, int> classes;
  std::vector<u64> parity(kw, 0);
  for (int j = 0; j < I.n; ++j) {
    std::vector<u64> key(kw, 0);
    bool zero = true;
    for (int a = 0; a < I.k; ++a) {
      if (I.bit(a, j)) {
        key[a / 64] |= u64{1} << (a % 64);
        parity[a / 64] ^= u64{1} << (a % 64);
        zero = false;
      }
    }
    if (zero) ++s.zero_columns;
    else ++classes[key];
  }
  for (const auto& kv : classes) s.duplicate_excess += kv.second - 1;
  s.parity_zero = std::all_of(parity.begin(), parity.end(), [](u64 w) { return w == 0; });
  s.duplicate_excess_ext = s.duplicate_excess + (!s.parity_zero && classes.count(parity) ? 1 : 0);
  return s;
}

struct TangentStats {
  int finite = 0;            // support positions with F'(alpha_j) parallel to F(alpha_j) != 0
  bool at_infinity = false;  // the same at infinity, with F(infinity) != 0
};

// A degenerate tangent at a point forces U' = 0 there for every U in W.
TangentStats degenerate_tangents(const GF& F, const Instance& I, const Curve& C) {
  const int n = I.n, D = I.D, N = C.stride;
  const std::vector<u16>& flog = C.vals;
  std::vector<u16> xs(I.pts.begin(), I.pts.begin() + D + 1);
  std::vector<std::vector<u16>> dval(I.k, std::vector<u16>(n));
  std::vector<u16> lead(I.k), sub(I.k);
  for (int a = 0; a < I.k; ++a) {
    std::vector<u16> ys(D + 1);
    for (int x = 0; x <= D; ++x) ys[x] = F.ex[flog[static_cast<std::size_t>(a) * N + x]];
    Poly f = interpolate(F, xs, ys);
    f.resize(D + 1, 0);
    lead[a] = f[D];
    sub[a] = f[D - 1];
    const Poly d = pderiv(f);
    for (int j = 0; j < n; ++j) dval[a][j] = peval(F, d, I.pts[j]);
  }
  TangentStats ts;
  for (int j = 0; j < n; ++j) {
    bool any = false, prop = true;
    u16 ratio = 0;
    bool have = false;
    for (int a = 0; a < I.k; ++a) {
      if (I.bit(a, j)) {
        any = true;
        if (!have) { ratio = dval[a][j]; have = true; }
        else if (dval[a][j] != ratio) prop = false;
      } else if (dval[a][j] != 0) {
        prop = false;
      }
    }
    if (any && prop) ++ts.finite;
  }
  // At infinity, in the parameter w = 1/Z, the value is the vector of degree-D
  // coefficients and the tangent the vector of degree-(D-1) coefficients.
  int piv = -1;
  for (int a = 0; a < I.k; ++a) {
    if (lead[a]) { piv = a; break; }
  }
  if (piv >= 0) {
    const u16 mu = F.div(sub[piv], lead[piv]);
    ts.at_infinity = true;
    for (int a = 0; a < I.k; ++a) {
      if (sub[a] != F.mul(mu, lead[a])) ts.at_infinity = false;
    }
  }
  return ts;
}

// (Z - alpha)^e in the representation of the curve: values at the points,
// or coefficients binom(e, i) alpha^(e - i), with binom(e, i) odd iff the
// binary digits of i are among those of e.
std::vector<u16> power_vector(const GF& F, const Instance& I, const Curve& C, u16 alpha, u64 e) {
  std::vector<u16> v(C.N, 0);
  if (!C.coeff) {
    for (int x = 0; x < C.N; ++x) v[x] = F.pow(static_cast<u16>(I.pts[x] ^ alpha), e);
  } else {
    for (u64 i = 0; i <= e && i < static_cast<u64>(C.N); ++i) {
      if ((i & ~e) == 0) v[i] = F.pow(alpha, e - i);
    }
  }
  return v;
}

// Interpolates selected curve polynomials and products from their values and
// checks deg f_a <= D and G^2 | (f_a f_b)'.
void interpolation_checks(const GF& F, const Instance& I, const Curve& C,
                          int samples, Rng& rng, bool& f_ok, bool& s_ok) {
  f_ok = s_ok = true;
  const int N = C.N;
  const Poly G2 = pmul(F, I.G, I.G);
  if (C.coeff) {
    f_ok = C.degree_ok;
    std::vector<u16> v(N);
    for (int s = 0; s < samples; ++s) {
      make_row(F, C, static_cast<int>(rng.below(I.k)), static_cast<int>(rng.below(I.k)), v.data());
      Poly du = pderiv(Poly(v.begin(), v.end()));
      pmod(F, du, G2);
      if (!du.empty()) s_ok = false;
    }
    return;
  }
  const std::vector<u16>& flog = C.vals;
  std::vector<u16> xs(I.pts.begin(), I.pts.end());
  for (int s = 0; s < samples; ++s) {
    const int a = static_cast<int>(rng.below(I.k));
    const int b = static_cast<int>(rng.below(I.k));
    std::vector<u16> ya(N), yu(N);
    for (int x = 0; x < N; ++x) {
      ya[x] = F.ex[flog[static_cast<std::size_t>(a) * N + x]];
      yu[x] = F.ex[flog[static_cast<std::size_t>(a) * N + x] + flog[static_cast<std::size_t>(b) * N + x]];
    }
    if (deg(interpolate(F, xs, ya)) > I.D) f_ok = false;
    Poly du = pderiv(interpolate(F, xs, yu));
    pmod(F, du, G2);
    if (!du.empty()) s_ok = false;
  }
}

// ------------------------------------------------------------ reporting

struct Options {
  int m = 0, t = 0;
  std::vector<int> ks;  // shortened dimensions to sweep (n = m t + k)
  int n_override = 0;
  u64 seed_lo = 1, seed_hi = 1;
  GoppaKind kind = GoppaKind::Irreducible;
  bool full = false;
  int positions = 4;  // -1 = all
  int threads = 0;
  bool diagnose = false;
  bool cross_rank = false;  // also compute dim span{f_a f_b : a < b}
  bool coeff = false;       // force the coefficient representation
  std::string label;
  std::string out;
  std::string dump;
  std::string scheme;     // e.g. Weis (d, s) = (7, 6)
  int full_n = 0;         // full length of the Classic McEliece code
};

std::string json_str(const std::string& s) { return "\"" + s + "\""; }

template <class T>
std::string json_list(const std::vector<T>& v) {
  std::ostringstream o;
  o << "[";
  for (std::size_t i = 0; i < v.size(); ++i) o << (i ? "," : "") << v[i];
  o << "]";
  return o.str();
}

void dump_instance(const Options& o, const Instance& I, const GF& F, int rank, int bound,
                   const std::vector<int>& positions, const std::vector<int>& r_ok,
                   const std::vector<int>& s_ok, const std::vector<int>& neg_in) {
  const std::string path = o.dump + "/" + o.label + "-n" + std::to_string(I.n) + "-s" +
                           std::to_string(I.seed) + ".txt";
  std::ofstream f(path);
  if (!f) die("cannot write " + path);
  f << "m " << I.m << "\nmodulus " << F.modulus << "\nt " << I.t << "\nn " << I.n
    << "\nk " << I.k << "\nD " << I.D << "\nrank " << rank << "\nbound " << bound << "\nG";
  for (u16 c : I.G) f << " " << c;
  f << "\nalpha";
  for (int j = 0; j < I.n; ++j) f << " " << I.pts[j];
  f << "\npositions";
  for (int p : positions) f << " " << p;
  f << "\nR_in_W";
  for (int x : r_ok) f << " " << x;
  f << "\nS_in_W";
  for (int x : s_ok) f << " " << x;
  f << "\nneg_in_W";
  for (int x : neg_in) f << " " << x;
  f << "\n";
}

int run_one(const Options& o, const GF& F, int n, u64 seed, std::FILE* out) {
  const double t0 = now_s();
  Instance I;
  std::string why;
  if (!build_instance(F, o.t, n, seed, o.kind, static_cast<u32>(o.full_n), o.coeff, I, why)) die(why);
  const double t_gen = now_s();
  const bool goppa_sq = check_goppa_square(F, I);
  const Curve curve = build_curve(F, I, std::max(1, o.threads));
  const double t_eval = now_s();

  const int bound = 2 * I.D - I.t + 1;
  const long P = static_cast<long>(I.k) * (I.k + 1) / 2;
  if (P > 60000000L) die("too many products");
  // The k squares f_a^2 come first, then the cross products f_a f_b, a < b,
  // each group in random order. In our runs the cross products span only a
  // subspace of codimension k in W, so a fully random order reaches the bound
  // only after the last square has appeared.
  Rng rng({static_cast<u32>(seed), static_cast<u32>(seed >> 32), static_cast<u32>(n),
           static_cast<u32>(o.t), static_cast<u32>(F.m), 0x5eed});
  std::vector<Pair> pairs, cross;
  pairs.reserve(P);
  cross.reserve(P - I.k);
  for (int a = 0; a < I.k; ++a) {
    pairs.emplace_back(static_cast<u16>(a), static_cast<u16>(a));
    for (int b = a + 1; b < I.k; ++b) cross.emplace_back(static_cast<u16>(a), static_cast<u16>(b));
  }
  rng.shuffle(pairs);
  rng.shuffle(cross);
  pairs.insert(pairs.end(), cross.begin(), cross.end());

  const int threads = I.N < 256 ? 1 : std::max(1, o.threads);
  const int cap = o.full ? static_cast<int>(std::min<long>(P, I.N)) : static_cast<int>(std::min<long>(P, bound));
  Basis B(F, I.N, cap);
  const RankRun rr = product_rank(F, curve, pairs, bound, o.full, threads, B);
  const double t_rank = now_s();
  const int rank = B.rank;
  const bool reached = rank == bound;
  const bool generic = rank == static_cast<int>(std::min<long>(P, bound));

  // Optional: rank of the cross products alone (every cross product is used).
  int cross_rank = -1;
  if (o.cross_rank) {
    Basis C(F, I.N, static_cast<int>(std::min<long>(static_cast<long>(cross.size()), I.N)));
    product_rank(F, curve, cross, bound, true, threads, C);
    cross_rank = C.rank;
  }
  const double t_cross = now_s();

  // Membership of (Z - a)^{2D}, (Z - a)^{2D-2} and the control (Z - a)^{2D-1}.
  // W is spanned by the basis unless the run stopped early, and then W = S.
  std::vector<int> positions;
  const bool all_positions = o.positions < 0 || (!reached && I.n <= 512);
  if (all_positions) {
    for (int j = 0; j < I.n; ++j) positions.push_back(j);
  } else {
    std::vector<int> idx(I.n);
    for (int j = 0; j < I.n; ++j) idx[j] = j;
    rng.shuffle(idx);
    positions.assign(idx.begin(), idx.begin() + std::min(o.positions, I.n));
    std::sort(positions.begin(), positions.end());
  }
  std::vector<int> r_ok, s_ok, neg_in;
  int r_count = 0, s_count = 0, both = 0, neg_count = 0;
  for (int j : positions) {
    const u16 a = I.pts[j];
    const int r = B.contains(power_vector(F, I, curve, a, 2 * static_cast<u64>(I.D)));
    const int s = B.contains(power_vector(F, I, curve, a, 2 * static_cast<u64>(I.D) - 2));
    const int ng = B.contains(power_vector(F, I, curve, a, 2 * static_cast<u64>(I.D) - 1));
    r_ok.push_back(r);
    s_ok.push_back(s);
    neg_in.push_back(ng);
    r_count += r;
    s_count += s;
    both += r && s;
    neg_count += ng;
  }
  const double t_member = now_s();

  bool f_ok = true, s_in = true;
  interpolation_checks(F, I, curve, 2, rng, f_ok, s_in);

  int sq_rank = -1;
  ColumnStats cs;
  TangentStats ts;
  bool have_tangents = false;
  if (!reached) {
    sq_rank = square_code_rank(I, pairs);
    cs = column_stats(I);
    if (o.diagnose || I.n <= 400) {
      ts = degenerate_tangents(F, I, curve);
      have_tangents = true;
    }
  }
  const double t_end = now_s();

  // The proposition: equality forces both memberships and excludes the control.
  const bool consistent = rank <= bound && neg_count == 0 && goppa_sq && f_ok && s_in &&
                          (!reached || both == static_cast<int>(positions.size()));

  std::ostringstream js;
  js << "{\"label\":" << json_str(o.label) << ",\"scheme\":" << json_str(o.scheme)
     << ",\"full_n\":" << o.full_n << ",\"m\":" << I.m << ",\"modulus\":" << F.modulus
     << ",\"t\":" << I.t << ",\"n\":" << I.n << ",\"k\":" << I.k << ",\"D\":" << I.D
     << ",\"N\":" << I.N << ",\"bound\":" << bound << ",\"products\":" << P
     << ",\"seed\":" << seed << ",\"goppa\":" << json_str(kind_name(o.kind))
     << ",\"representation\":" << json_str(I.coeff ? "coefficients" : "values")
     << ",\"mode\":" << json_str(o.full ? "full" : "stop-at-bound")
     << ",\"rank\":" << rank << ",\"reached_bound\":" << (reached ? "true" : "false")
     << ",\"generic\":" << (generic ? "true" : "false")
     << ",\"deficit\":" << (std::min<long>(P, bound) - rank)
     << ",\"rows_used\":" << rr.rows_used
     << ",\"cross_rank\":" << cross_rank
     << ",\"positions_tested\":" << positions.size() << ",\"R_in_W\":" << r_count
     << ",\"S_in_W\":" << s_count << ",\"R_and_S\":" << both
     << ",\"control_in_W\":" << neg_count;
  if (!reached) {
    // Predicted deficit: one per extra copy of a point (double points), two per
    // base point, one per degenerate tangent; infinity is included.
    const int predicted = cs.duplicate_excess_ext + 2 * (cs.zero_columns + cs.parity_zero) +
                          ts.finite + ts.at_infinity;
    js << ",\"square_code_rank\":" << sq_rank << ",\"zero_columns\":" << cs.zero_columns
       << ",\"duplicate_excess\":" << cs.duplicate_excess
       << ",\"parity_zero\":" << (cs.parity_zero ? "true" : "false")
       << ",\"duplicate_excess_ext\":" << cs.duplicate_excess_ext
       << ",\"degenerate_tangents\":" << (have_tangents ? ts.finite : -1)
       << ",\"cusp_at_infinity\":" << (have_tangents ? (ts.at_infinity ? 1 : 0) : -1)
       << ",\"predicted_deficit\":" << (have_tangents ? predicted : -1);
  }
  js << ",\"checks\":{\"goppa_square\":" << (goppa_sq ? "true" : "false")
     << ",\"f_degree\":" << (f_ok ? "true" : "false")
     << ",\"product_in_S\":" << (s_in ? "true" : "false")
     << ",\"consistent\":" << (consistent ? "true" : "false") << "}"
     << ",\"threads\":" << threads
     << ",\"seconds\":{\"generate\":" << (t_gen - t0) << ",\"evaluate\":" << (t_eval - t_gen)
     << ",\"rank\":" << (t_rank - t_eval) << ",\"cross_rank\":" << (t_cross - t_rank)
     << ",\"membership\":" << (t_member - t_cross)
     << ",\"total\":" << (t_end - t0) << "}}";
  std::fprintf(out, "%s\n", js.str().c_str());
  std::fflush(out);
  std::fprintf(stderr,
               "%s m=%d t=%d n=%d k=%d seed=%" PRIu64 ": rank %d / bound %d (%s), rows %ld, "
               "R&S %d/%zu, control %d, %.2fs%s\n",
               o.label.c_str(), I.m, I.t, I.n, I.k, seed, rank, bound,
               reached ? "reached" : "DEFICIENT", rr.rows_used, both, positions.size(), neg_count,
               t_end - t0, consistent ? "" : "  INCONSISTENT");
  if (!o.dump.empty()) dump_instance(o, I, F, rank, bound, positions, r_ok, s_ok, neg_in);
  return consistent ? 0 : 1;
}

// ----------------------------------------------------------- self-test

int self_test() {
  int failures = 0;
  auto expect = [&](bool c, const std::string& what) {
    if (!c) {
      std::fprintf(stderr, "self-test FAILED: %s\n", what.c_str());
      ++failures;
    }
  };
  // Field arithmetic against carry-less multiplication.
  for (int m = 2; m <= 15; ++m) {
    GF F(m, default_modulus(m));
    Rng rng({static_cast<u32>(m), 1});
    for (int i = 0; i < 2000; ++i) {
      const u16 a = static_cast<u16>(rng.below(F.size)), b = static_cast<u16>(rng.below(F.size)),
                c = static_cast<u16>(rng.below(F.size));
      expect(F.mul(a, b) == GF::slow_mul(a, b, m, F.modulus), "mul m=" + std::to_string(m));
      expect(F.mul(a, b ^ c) == (F.mul(a, b) ^ F.mul(a, c)), "distributivity");
      if (a) expect(F.mul(a, F.inv(a)) == 1, "inverse");
    }
  }
  // Number of monic irreducible polynomials of degree t over GF(q).
  {
    auto moebius = [](int d) {
      int r = 1;
      for (int p = 2; p * p <= d; ++p) {
        if (d % p == 0) {
          d /= p;
          if (d % p == 0) return 0;
          r = -r;
        }
      }
      if (d > 1) r = -r;
      return r;
    };
    for (int m : {2, 3}) {
      GF F(m, default_modulus(m));
      const long q = F.size;
      for (int t = 2; t <= (m == 2 ? 4 : 3); ++t) {
        long expected = 0, qp = 1;
        for (int d = 1; d <= t; ++d) {
          if (t % d) continue;
          qp = 1;
          for (int e = 0; e < t / d; ++e) qp *= q;
          expected += moebius(d) * qp;
        }
        expected /= t;
        long count = 0, total = 1;
        for (int e = 0; e < t; ++e) total *= q;
        for (long idx = 0; idx < total; ++idx) {
          Poly g(t + 1);
          long r = idx;
          for (int i = 0; i < t; ++i) { g[i] = static_cast<u16>(r % q); r /= q; }
          g[t] = 1;
          count += irreducible(F, g);
        }
        expect(count == expected, "irreducible count q=" + std::to_string(q) + " t=" + std::to_string(t));
      }
    }
  }
  // Rank engine against plain dense elimination of all products, with the
  // batch path forced to use several threads.
  struct Case { int m, t, n; GoppaKind kind; };
  const Case cases[] = {{6, 3, 26, GoppaKind::SquareFree}, {6, 3, 30, GoppaKind::Irreducible},
                        {7, 4, 38, GoppaKind::Irreducible}, {8, 5, 52, GoppaKind::SquareFree}};
  int instances = 0, deficient = 0;
  for (const Case& cs : cases) {
    GF F(cs.m, default_modulus(cs.m));
    for (u64 seed = 1; seed <= 12; ++seed) {
      Instance I;
      std::string why;
      if (!build_instance(F, cs.t, cs.n, seed, cs.kind, 0, false, I, why)) die(why);
      ++instances;
      expect(!I.coeff, "value representation for small codes");
      expect(check_goppa_square(F, I), "Gamma(G) = Gamma(G^2)");
      const Curve curve = build_curve(F, I, 2);
      const std::vector<u16>& flog = curve.vals;
      Rng rng({static_cast<u32>(seed), 7});
      bool f_ok, s_in;
      interpolation_checks(F, I, curve, 6, rng, f_ok, s_in);
      expect(f_ok, "deg f_a <= D");
      expect(s_in, "G^2 | (f_a f_b)'");
      const int N = I.N, bound = 2 * I.D - I.t + 1;
      std::vector<Pair> pairs;
      for (int a = 0; a < I.k; ++a) {
        for (int b = a; b < I.k; ++b) pairs.emplace_back(static_cast<u16>(a), static_cast<u16>(b));
      }
      // dense reference
      std::vector<std::vector<u16>> M;
      for (const Pair& pr : pairs) {
        std::vector<u16> row(N);
        for (int x = 0; x < N; ++x) {
          row[x] = F.mul(F.ex[flog[static_cast<std::size_t>(pr.first) * N + x]],
                         F.ex[flog[static_cast<std::size_t>(pr.second) * N + x]]);
        }
        M.push_back(row);
      }
      int dense = 0;
      for (int c = 0; c < N && dense < static_cast<int>(M.size()); ++c) {
        int p = -1;
        for (int i = dense; i < static_cast<int>(M.size()); ++i) {
          if (M[i][c]) { p = i; break; }
        }
        if (p < 0) continue;
        std::swap(M[p], M[dense]);
        const u16 il = F.inv(M[dense][c]);
        for (int x = 0; x < N; ++x) M[dense][x] = F.mul(M[dense][x], il);
        for (int i = 0; i < static_cast<int>(M.size()); ++i) {
          if (i == dense || !M[i][c]) continue;
          const u16 f = M[i][c];
          for (int x = 0; x < N; ++x) M[i][x] ^= F.mul(f, M[dense][x]);
        }
        ++dense;
      }
      expect(dense <= bound, "dim W <= 2D - t + 1");
      Basis Bf(F, N, static_cast<int>(std::min<std::size_t>(pairs.size(), N)));
      product_rank(F, curve, pairs, bound, true, 3, Bf);
      expect(Bf.rank == dense, "engine full rank = dense rank");
      Basis Bs(F, N, std::min(static_cast<int>(pairs.size()), bound));
      product_rank(F, curve, pairs, bound, false, 3, Bs);
      expect(Bs.rank == std::min(dense, bound), "early stop rank");
      if (dense < bound) ++deficient;
      // The same key in coefficient form: products are multiplied as
      // polynomials instead of evaluated, an independent computation.
      Instance Ic;
      if (!build_instance(F, cs.t, cs.n, seed, cs.kind, 0, true, Ic, why)) die(why);
      expect(Ic.coeff && Ic.G == I.G && std::equal(Ic.pts.begin(), Ic.pts.end(), I.pts.begin()),
             "same key in coefficient form");
      const Curve cc = build_curve(F, Ic, 2);
      expect(cc.degree_ok, "coefficient form: deg f_a <= D");
      bool cf_ok, cs_in;
      interpolation_checks(F, Ic, cc, 6, rng, cf_ok, cs_in);
      expect(cs_in, "coefficient form: G^2 | (f_a f_b)'");
      Basis Bc(F, N, static_cast<int>(std::min<std::size_t>(pairs.size(), N)));
      product_rank(F, cc, pairs, bound, true, 2, Bc);
      expect(Bc.rank == dense, "coefficient-form rank = dense rank");
      for (int j = 0; j < I.n; ++j) {
        for (u64 e : {2 * static_cast<u64>(I.D), 2 * static_cast<u64>(I.D) - 2, 2 * static_cast<u64>(I.D) - 1}) {
          const bool in_values = Bf.contains(power_vector(F, I, curve, I.pts[j], e));
          const bool in_coeffs = Bc.contains(power_vector(F, Ic, cc, I.pts[j], e));
          expect(in_values == in_coeffs, "membership agrees between representations");
        }
        const auto r = power_vector(F, I, curve, I.pts[j], 2 * static_cast<u64>(I.D));
        const auto ng = power_vector(F, I, curve, I.pts[j], 2 * static_cast<u64>(I.D) - 1);
        expect(Bf.contains(r) == Bs.contains(r) || dense < bound, "membership agreement");
        expect(!Bf.contains(ng), "control (Z-a)^{2D-1} not in W");
        if (dense == bound) expect(Bf.contains(r), "(Z-a)^{2D} in W at full rank");
      }
    }
  }
  std::fprintf(stderr, "self-test: %d instances (%d deficient), %s\n", instances, deficient,
               failures ? "FAILED" : "ok");
  return failures ? 1 : 0;
}

// ------------------------------------------------------------------ CLI

struct Preset {
  const char* name;
  int m, t, k, full_n;
  const char* scheme;
};

// Shortened dimensions k_l of the distinguisher, from GIAJS (d, s) = (8, 7)
// and Weis (d, s) = (7, 6), as in the cost table of the paper.
const Preset kPresets[] = {
    {"mceliece348864-weis", 12, 64, 154, 3488, "Weis (d,s)=(7,6)"},
    {"mceliece460896-weis", 13, 96, 195, 4608, "Weis (d,s)=(7,6)"},
    {"mceliece6688128-weis", 13, 128, 224, 6688, "Weis (d,s)=(7,6)"},
    {"mceliece6960119-weis", 13, 119, 216, 6960, "Weis (d,s)=(7,6)"},
    {"mceliece8192128-weis", 13, 128, 224, 8192, "Weis (d,s)=(7,6)"},
    {"mceliece348864-giajs", 12, 64, 217, 3488, "GIAJS (d,s)=(8,7)"},
    {"mceliece460896-giajs", 13, 96, 274, 4608, "GIAJS (d,s)=(8,7)"},
    {"mceliece6688128-giajs", 13, 128, 315, 6688, "GIAJS (d,s)=(8,7)"},
    {"mceliece6960119-giajs", 13, 119, 304, 6960, "GIAJS (d,s)=(8,7)"},
    {"mceliece8192128-giajs", 13, 128, 315, 8192, "GIAJS (d,s)=(8,7)"},
};

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
               "  product_rank self-test\n"
               "  product_rank run --m M --t T --k K[-K2] [options]\n"
               "  product_rank preset NAME [options]\n"
               "  product_rank presets\n"
               "options:\n"
               "  --seeds A[-B]      seeds (default 1)\n"
               "  --n N              code length instead of n = m t + k\n"
               "  --goppa KIND       irreducible (default) or squarefree\n"
               "  --full             process every product (no stop at the bound)\n"
               "  --positions P|all  positions tested for R and S (default 4)\n"
               "  --threads T        worker threads (default: hardware)\n"
               "  --diagnose         degenerate-tangent count for deficient instances\n"
               "  --cross-rank       also compute the rank of the products f_a f_b, a < b\n"
               "  --coeff            products in coefficient form (automatic if 2D+1 > 2^m)\n"
               "  --label L  --out FILE  --dump DIR\n");
  std::exit(2);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) usage();
  const std::string cmd = argv[1];
  if (cmd == "self-test") return self_test();
  if (cmd == "presets") {
    for (const Preset& p : kPresets) {
      std::printf("%-24s m=%d t=%d k=%d n=%d (full n=%d), %s\n", p.name, p.m, p.t, p.k,
                  p.m * p.t + p.k, p.full_n, p.scheme);
    }
    return 0;
  }
  Options o;
  int argi = 2;
  if (cmd == "preset") {
    if (argc < 3) usage();
    const std::string name = argv[2];
    const Preset* found = nullptr;
    for (const Preset& p : kPresets) {
      if (name == p.name) found = &p;
    }
    if (!found) die("unknown preset " + name);
    o.m = found->m;
    o.t = found->t;
    o.ks = {found->k};
    o.label = found->name;
    o.scheme = found->scheme;
    o.full_n = found->full_n;
    argi = 3;
  } else if (cmd != "run") {
    usage();
  }
  for (; argi < argc; ++argi) {
    const std::string a = argv[argi];
    auto val = [&]() -> std::string {
      if (argi + 1 >= argc) die("missing value for " + a);
      return argv[++argi];
    };
    long lo, hi;
    if (a == "--m") o.m = std::stoi(val());
    else if (a == "--t") o.t = std::stoi(val());
    else if (a == "--k") {
      parse_range(val(), lo, hi);
      o.ks.clear();
      for (long k = lo; k <= hi; ++k) o.ks.push_back(static_cast<int>(k));
    } else if (a == "--n") o.n_override = std::stoi(val());
    else if (a == "--seeds") {
      parse_range(val(), lo, hi);
      o.seed_lo = static_cast<u64>(lo);
      o.seed_hi = static_cast<u64>(hi);
    } else if (a == "--goppa") {
      const std::string v = val();
      if (v == "irreducible") o.kind = GoppaKind::Irreducible;
      else if (v == "squarefree") o.kind = GoppaKind::SquareFree;
      else die("unknown Goppa kind " + v);
    } else if (a == "--full") o.full = true;
    else if (a == "--positions") {
      const std::string v = val();
      o.positions = v == "all" ? -1 : std::stoi(v);
    } else if (a == "--threads") o.threads = std::stoi(val());
    else if (a == "--diagnose") o.diagnose = true;
    else if (a == "--cross-rank") o.cross_rank = true;
    else if (a == "--coeff") o.coeff = true;
    else if (a == "--label") o.label = val();
    else if (a == "--out") o.out = val();
    else if (a == "--dump") o.dump = val();
    else die("unknown option " + a);
  }
  if (o.m == 0 || o.t == 0 || (o.ks.empty() && o.n_override == 0)) usage();
  if (o.threads <= 0) o.threads = std::max(1u, std::thread::hardware_concurrency());
  if (o.label.empty()) o.label = "m" + std::to_string(o.m) + "-t" + std::to_string(o.t);
  GF F(o.m, default_modulus(o.m));
  std::FILE* out = stdout;
  if (!o.out.empty()) {
    out = std::fopen(o.out.c_str(), "a");
    if (!out) die("cannot open " + o.out);
  }
  std::vector<int> lengths;
  if (o.n_override) lengths.push_back(o.n_override);
  else for (int k : o.ks) lengths.push_back(o.m * o.t + k);
  int bad = 0;
  for (int n : lengths) {
    for (u64 s = o.seed_lo; s <= o.seed_hi; ++s) bad += run_one(o, F, n, s, out);
  }
  if (out != stdout) std::fclose(out);
  return bad ? 1 : 0;
}
