#ifndef ZENO_INTERNAL_HPP
#define ZENO_INTERNAL_HPP


#include <cstdint>

namespace zeno
{

namespace internal
{


// Computes factorial
template<typename Z = int64_t>
Z fact(int n) {
    static std::vector<Z> _factorial({1});
    while(_factorial.size() <= n)
        _factorial.push_back(_factorial.back() * _factorial.size());
    return _factorial[n];
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
    int64_t x = 0;
    while((uZ(1) << x) < (Z)(n)) x++;    
    return x;
}

// int64_t ceil_log2(int64_t n) {
//     int64_t x = 0;
//     while((1u << x) < (uint64_t)(n)) x++;    
//     return x;
// }



}; // namespace internal

}; // namespace zeno


#endif /* ZENO_INTERNAL_HPP */