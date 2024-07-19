#ifndef ZENO_INTERNAL_HPP
#define ZENO_INTERNAL_HPP

#include <cstdint>

namespace zeno
{

namespace internal
{


/// @brief Computes n! (where n! = 0 for n < 0).
template<typename Z = int64_t>
Z fact(int n) {
    if(n < 0) return Z(0);
    static std::vector<Z> _factorial({1});
    while(_factorial.size() <= n)
        _factorial.push_back(_factorial.back() * _factorial.size());
    return _factorial[n];
}

/// @brief Computes 1 / n! (where 1/n! = 0 for n < 0).
template<typename Z = int64_t>
Z inv_fact(int n) {
    if(n < 0) return Z(0);
    static std::vector<Z> _inv_factorial({1});
    while(_inv_factorial.size() <= n)
        _inv_factorial.push_back(_inv_factorial.back() / _inv_factorial.size());
    return _inv_factorial[n];
}

/// @brief Compute the binomial coefficient C(n, m) (where C(n, m) = 0 for n or m < 0).
template<typename Z = int64_t>
Z binom(int n, int m) {
    if(n < 0 || m < 0 || n - m < 0) return Z(0);
    if(std::is_integral_v<Z>) {
        static std::vector<std::vector<Z>> _binom(std::vector({}));
        while(_binom.size() <= n)
            _binom.push_back({});
        while(_binom[n].size() <= m)
            _binom[n].push_back(-1);
        if(_binom[n][m] == -1)
            return binom[n][m] = binom(n, m-1) + binom(n-1, m-1);
        return binom[n][m];
    }
    return fact(n) * inv_fact(m) * inv_fact(n - m);
}

template<class Z = int64_t>
Z binary_exponentiation_mod(Z base, Z e, Z mod) {
    Z result = 1;
    base %= mod;
    while (e > 0) {
        if (e & 1) result = (__int128_t)result * base % mod;
        base = (__int128_t)base * base % mod;
        e >>= 1;
    }
    return result;
}

template<class Z = int64_t>
constexpr Z binary_exponentiation_mod_constexpr(Z base, Z e, Z mod) {
    Z result = 1;
    base %= mod;
    while (e > 0) {
        if (e & 1) result = (__int128_t)result * base % mod;
        base = (__int128_t)base * base % mod;
        e >>= 1;
    }
    return result;
}

/* G should form a group, Z should be an integer struct */
template<class G = int64_t, class Z = int64_t>
G binary_exponentiation(G g, Z n) {
    if(n < 0) g = 1/g, n = -n;
    G x = g, y = G(1);
    while(n > 0) {
        if(n % 2 == 1) y *= x;
        n /= 2, x *= x;
    }
    return y;
}

template<class G = int64_t, class Z = int64_t>
G pow(G g, Z n) {
    return binary_exponentiation(g, n);
}

template<typename Z = int64_t>
Z sqr(Z n) { return n * n; }





template<typename Z = int64_t>
Z ceil_log2(Z n) {
    using uZ = std::make_unsigned_t<Z>;
    uZ x = 0;
    while((uZ(1) << x) < (Z)(n)) x++;    
    return x;
}

template<typename Z = int64_t>
Z popcount(Z x) {
    Z ret = 0;
    while(x) ret += (x&1), x >>= 1;
    return ret;
}


// int64_t ceil_log2(int64_t n) {
//     int64_t x = 0;
//     while((1u << x) < (uint64_t)(n)) x++;    
//     return x;
// }



}; // namespace internal

}; // namespace zeno


#endif /* ZENO_INTERNAL_HPP */