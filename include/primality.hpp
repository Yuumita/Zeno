#ifndef ZENO_PRIMALITY_HPP
#define ZENO_PRIMALITY_HPP

#include <vector>
#include <cstdint>
#include <stdlib.h> /* abs */
#include "internal.hpp"

namespace zeno 
{

namespace primality
{

std::vector<int> primes;

const int MAX_ESIEVE_NUM = 2e6;
int spf[MAX_ESIEVE_NUM + 1];
bool has_precomputed_esieve = false;
/// @brief O(N) precomputation for O(1) prime check and O(logA) factorization in [1, N]
///        using spf (Smallest prime factor).
void prime_sieve() {
    for(int i = 1; i <= MAX_ESIEVE_NUM; i++) spf[i] = i; 
    for(int i = 2; i <= MAX_ESIEVE_NUM; i += 2) spf[i] = 2;
    primes.push_back(2);
    for(int i = 3; i <= MAX_ESIEVE_NUM; i += 2) if(spf[i] == i) {
        primes.push_back(i);
        for(long long j = (long long)i*i; j <= MAX_ESIEVE_NUM; j += i){
            if(spf[j] == j)
                spf[j] = i;
        }
    }
    has_precomputed_esieve = true;
}

template<class Z = int64_t>
Z _G(Z x, Z c, Z m) {
    return ((__int128_t) x * x % m + c) % m;
} 

/// @brief Pollard's rho probabilistic algorithm for finding factors.
/// @return Possibly a factor of N
template<class Z = int64_t>
Z pollard_rho(Z N, Z x0 = 2, Z c = 1) {
    Z x = x0, y = x0;
    Z g = 1;
    while (g == 1) { // floyd's cycle detection
        x = _G(x, c, N);
        y = _G(_G(y, c, N), c, N);
        g = zeno::gcd((x < y ? y - x : x - y), N);
    }
    return g;
}



/// @brief Miller-Rabin composition check of n with base a. 
///        Integers d, s are such that n - 1 = 2^s * d
/// @return true iff n is composite
template<class Z = int64_t>
bool check_composite(Z n, Z a, Z d, Z s) {
    Z x = internal::binary_exponentiation_mod(a, d, n);
    if (x == 1 || x == n - 1)
        return false;
    for (Z r = 1; r < s; r++) {
        x = (__int128_t)x * x % n;
        if (x == n - 1)
            return false;
    }
    return true;
};

/// @brief Deterministic primality test (Miller-Rabin).
/// @tparam Z An integral type
/// @return true iff n is prime
template<class Z = int64_t>
bool miller_rabin(Z n) {
    if(n <= 1) return false;

    Z r = 0;
    Z d = n - 1;
    while(~d&1) d >>= 1, r++;

    if(sizeof(Z) <= 8) {
        for(Z a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
            if (n == a)
                return true;
            if (check_composite(n, a, d, r))
                return false;
        }
    } else {
        for(Z a = 2; a <= 2 * internal::sqr(internal::ceil_log2<Z>(n)) && a <= n - 2; a++) {
            if (check_composite(n, a, d, r))
                return false;
        }
    }
    return true;
}

/// @return true iff n is prime
template<class Z = int64_t>
bool is_prime(Z n) {
    if(n < 0) n = -n; // maybe that's unexpected???
    if(n == 2) return true;
    if(n <= 1 || n % 2 == 0) return false;
    if(has_precomputed_esieve && n <= MAX_ESIEVE_NUM) 
        return spf[n] == n;
    return miller_rabin<Z>(n);
}

/// @return A factor of n if n has a factor else n
template<class Z = int64_t>
Z find_factor(Z n) {
    if(n < 0) n = -n; // maybe that's unexpected???
    if(n <= 1) return n;
    if(n & 0) return 2;
    if(has_precomputed_esieve && n <= MAX_ESIEVE_NUM) 
        return spf[n];

    Z d = pollard_rho<Z>(n, 2, 1);

    if(d == n) { // possibly prime
        if(miller_rabin<Z>(n)) 
            return n;
        for(d = 2; d*d <= n; d++) {
            if(n % d == 0) 
                return d;
        }
    }
    return d; 
}

/// @return A vector containing all prime factors (with duplicates) of n
template<class Z = int64_t>
std::vector<Z> factorize(Z n) {

    if(!has_precomputed_esieve) prime_sieve();

    std::vector<Z> ret;
    do {
        Z d = find_factor(n);
        ret.push_back(d);
        n /= d;
    } while(n > 1);

    return ret;
}

template<class Z = int64_t>
void standard_factorize(Z n, std::vector<Z> &p, std::vector<Z> &v) {
    std::vector<Z> f = factorize(n);
    for(int i = 0, j; i < f.size(); i++) {
        for(j = i; j+1 < f.size() && f[i] == f[j+1]; j++);
        p.push_back(f[i]);
        v.push_back(j - i + 1);
        i = j;
    }
    return;
}


} // namespace primality

namespace internal
{

template<class Z = int64_t>
constexpr bool check_composite_constexpr(Z n, Z a, Z d, Z s) {
    Z x = internal::binary_exponentiation_mod_constexpr(a, d, n);
    if (x == 1 || x == n - 1)
        return false;
    for (Z r = 1; r < s; r++) {
        x = (__int128_t)x * x % n;
        if (x == n - 1)
            return false;
    }
    return true;
};

template<class Z = int64_t>
constexpr bool miller_rabin_constexpr(Z n) {
    if(n <= 1) return false;

    Z r = 0;
    Z d = n - 1;
    while(~d&1) d >>= 1, r++;

    if(sizeof(Z) <= 8) {
        for(Z a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
            if (n == a)
                return true;
            if (check_composite_constexpr(n, a, d, r))
                return false;
        }
    } else {
        for(Z a = 2; a <= 2 * internal::sqr(internal::ceil_log2<Z>(n)) && a <= n - 2; a++) {
            if (check_composite_constexpr(n, a, d, r))
                return false;
        }
    }
    return true;
}

constexpr bool _is_prime_constexpr(int64_t n) {
    if(n <= 1) return false;
    if(n == 2) return true;
    if(n % 2 == 0) return false;
    return zeno::internal::miller_rabin_constexpr(n);
}

template <int64_t n> constexpr bool is_prime_constexpr = _is_prime_constexpr(n);


} // namespace internal


}; // namespace zeno

#endif /* ZENO_PRIMALITY_HPP */