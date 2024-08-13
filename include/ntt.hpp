#pragma once

#include <vector>
#include <cstdint>
#include <stdlib.h> /* abs */
#include <type_traits>
#include <unordered_map>
#include "internal.hpp"
#include "internal_fft.hpp"
#include "modular.hpp"
#include "groups.hpp"
#include <algorithm> // std::swap

namespace zeno
{

namespace fft
{
    


/// @brief The NTT class. Asserts that for the given prime number there are primitive roots of unity.
/// @tparam M Modular ring class.
/// @tparam Z Integer ring holding numbers at least as large as the cardinality of M.
template <class M, typename Z = int64_t, std::enable_if_t<zeno::is_modular_v<M>>* = nullptr>
class NumberTheoreticTransform {
    using NTT = NumberTheoreticTransform;

    int ordlog; // log(p-1), needed for finding the biggest k such that there is a (2^k)-th unity root
    Z p;        // odd prime p = c2^k + 1
    M g;        // primitive unity root of p
    Z c;

    mutable std::vector<M> w, iw; // w[i] is the 2^i root of unity and iw[i] is its inverse

    NumberTheoreticTransform(Z pp) : p(pp), ordlog(__builtin_ctzll(pp-1)) {
        g = M(zeno::primitive_root(p));
        assert(g != M(0));
        w.resize(ordlog + 1);
        iw.resize(ordlog + 1);

        
    
        c = (p-1)/(Z(1) << ordlog);
        w[ordlog] = g.pow(c); 
        iw[ordlog] = w[ordlog].inv();
        for(int i = ordlog - 1; i >= 0; i--) {
            w[i] = w[i + 1] * w[i + 1];
            iw[i] = w[i].inv();
        }

        assert((p - 1) / (Z(1) << ordlog) * (Z(1) << ordlog) + 1 == p);

        for(int j = 1; j <= ordlog; j++) {
            assert(w[j].pow((1ll << j)) == M(1));
            for(int i = 1; i < j; i++) {
                assert(w[j].pow((1ll << i)) != M(1));
            }
        }

        assert(w[0] == M(1));

    }

    NumberTheoreticTransform() : NumberTheoreticTransform(M::mod()) {}
    

public:

    int get_ordlog() { return ordlog; }
    Z get_p() { return p; }
    M get_g() { return g; }
    Z get_c() { return c; }


    template <typename U = M>
    static const NTT& get_info(Z m = 0) {
    // static std::enable_if_t<!std::is_same<U, static_modular<>::value, const NTT&> get_info(Z m = 0) {
        static NTT info(M::mod()); 
        return info;
    }

    // template <typename U = M>
    // static std::enable_if_t<std::is_same<U, dynamic_modular<Z>>::value, NTT&> get_info(Z m) {
    //     static std::unordered_map<Z, NTT> instances;
    //     auto it = instances.find(m);
    //     if (it != instances.end()) {
    //         return it->second;
    //     } else {
    //         NTT info(m);
    //         return instances[m] = info;
    //     }
    // }

    const std::vector<M> &roots()  const { return w; }
    const std::vector<M> &iroots() const { return iw; }


    void transform(std::vector<M> &a, bool inverse = false) {
        size_t n = a.size(), s = 0;
        if(n <= 1) return;
        while((1 << s) < n) s++;
        assert((n & (n - 1)) == 0); // power of 2

        // bit-reverse
        for(int i = 1, j = 0; i < n; i++) { // j = bit_reverse(i, log n)
            // j := bit_reverse(i+1) == this_loop(bit_reverse(i)): from msb to lsb flip the trailing 1s and the next 0
            int bit = n >> 1;
            for(; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;

            if(i < j) 
                std::swap(a[i], a[j]);
        }

        for (size_t l = 1; l <= s; l++) {
            size_t len = size_t(1) << l;
            M wlen = inverse ? iw[l] : w[l];
            for (size_t i = 0; i < n; i += len) {
                M w(1);
                for(size_t j = 0; j < len / 2; j++) {
                    M u = a[i + j];
                    M v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }

        if(inverse) {
            for(M &e: a)
                e /= M(n);
        }
    }
    inline void inv_transform(std::vector<M> &a) {
        static_assert(zeno::is_modular_v<M>, "Template parameter M must be modular (modint).");
        transform(a, true);
    }
};


/// @brief Applies the (non-cyclic) convolution of a and b. Assertion rises if 
///        the modular class M does not have enough roots of unity for the result.
/// @tparam M The modular class.
/// @return The convolution of a and b.
template <class M, std::enable_if_t<zeno::is_modular_v<M>>* = nullptr>
std::vector<M> convolution_ntt(std::vector<M> const &a, std::vector<M> const &b) {

    if(a.empty() || b.empty()) return {};

    size_t n = size_t(a.size()), m = size_t(b.size());
    size_t N = zeno::fft::compute_convolution_size(n, m);

    NumberTheoreticTransform<M> ntt = NumberTheoreticTransform<M>::get_info(a[0].mod());
    assert(N <= (size_t(1) << ntt.roots().size()));

    std::vector<M> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    fa.resize(N, M(0));
    fb.resize(N, M(0));

    ntt.transform(fa);
    ntt.transform(fb);
    for(int i = 0; i < N; i++) 
        fa[i] *= fb[i];
    ntt.inv_transform(fa);

    fa.resize(n + m - 1);
    return fa;
}

} // namespace fft
    
} // namespace zeno