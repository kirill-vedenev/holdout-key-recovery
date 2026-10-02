// Adapted from Markku-Juhani O. Saarinen's holdout_cuda_worker.cu.
// Copyright (c) 2026 Markku-Juhani O. Saarinen. MIT license: LICENSE.txt.
// Changes: standalone namespace, bounds checks, and a shared host/device table API.
#pragma once
#include <array>
#include <cstdint>
#include <stdexcept>

namespace saarinen {
using u64 = std::uint64_t;
constexpr int MAX_K = 63;
constexpr int MAX_D = 8;
struct HostChoose {
    std::array<u64, (MAX_K + 1) * (MAX_D + 1)> values{};
    HostChoose() {
        for (int n = 0; n <= MAX_K; ++n) {
            at(n, 0) = 1;
            for (int r = 1; r <= n && r <= MAX_D; ++r)
                at(n, r) = r == n ? 1 : get(n - 1, r - 1) + get(n - 1, r);
        }
    }
    u64 get(int n, int r) const {
        if (n < 0 || n > MAX_K || r < 0 || r > MAX_D || r > n) return 0;
        return values[std::size_t(n) * (MAX_D + 1) + r];
    }
    u64& at(int n, int r) { return values[std::size_t(n) * (MAX_D + 1) + r]; }
};
inline const HostChoose HOST_CHOOSE;
inline u64 colex_rank_host(u64 mask) {
    u64 rank = 0;
    int index = 1;
    while (mask) {
        const int bit = __builtin_ctzll(mask);
        mask &= mask - 1;
        rank += HOST_CHOOSE.get(bit, index++);
    }
    return rank;
}
inline u64 colex_unrank_host(u64 rank, int size, int k) {
    if (k < 0 || k > MAX_K || size < 0 || size > MAX_D || size > k
        || rank >= HOST_CHOOSE.get(k, size)) throw std::runtime_error("invalid colex index");
    u64 mask = 0;
    int upper = k;
    for (int index = size; index >= 1; --index) {
        int bit = upper - 1;
        while (bit >= index && HOST_CHOOSE.get(bit, index) > rank) --bit;
        mask |= u64{1} << bit;
        rank -= HOST_CHOOSE.get(bit, index);
        upper = bit;
    }
    return mask;
}
} // namespace saarinen
