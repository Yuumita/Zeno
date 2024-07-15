#ifndef ZENO_NTT_HPP
#define ZENO_NTT_HPP

#include <vector>
#include <cstdint>
#include <stdlib.h> /* abs */
#include <type_traits>
#include <unordered_map>
#include "internal.hpp"
#include "internal_fft.hpp"
#include "modular.hpp"
#include "convolution.hpp"
#include "cyclic_groups.hpp"
#include <algorithm> // std::swap

namespace zeno
{

namespace fft
{
    


template <class M, typename Z = int64_t>
class NumberTheoreticTransform {
    using NTT = NumberTheoreticTransform;
    static_assert(zeno::is_modular_v<M>, "Template parameter M must be modular (modint).");

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

        // assert((p - 1) / (Z(1) << ordlog) * (Z(1) << ordlog) + 1 == p);
        
        c = (p-1)/(Z(1) << ordlog);
        w[ordlog] = g.pow(c); 
        iw[ordlog] = w[ordlog].inv();
        for(int i = ordlog - 1; i >= 0; i--) {
            w[i] = w[i + 1] * w[i + 1];
            iw[i] = w[i].inv();
        }

        assert(w[0] == M(1));

    }

    NumberTheoreticTransform() {
        NumberTheoreticTransform(M::mod());
    }

public:

    int get_ordlog() { return ordlog; }
    Z get_p() { return p; }
    M get_g() { return g; }
    Z get_c() { return c; }


    template <typename U = M>
    static std::enable_if_t<!std::is_same<U, dynamic_modular<Z>>::value, const NTT&> get_info(Z m = 0) {
        static NTT info(M::mod()); 
        return info;
    }

    template <typename U = M>
    static std::enable_if_t<std::is_same<U, dynamic_modular<Z>>::value, NTT&> get_info(Z m) {
        static std::unordered_map<Z, NTT> instances;
        auto it = instances.find(m);
        if (it != instances.end()) {
            return it->second;
        } else {
            NTT info(m);
            return instances[m] = info;
        }
    }

    const std::vector<M> &roots() const {
        return w;
    }

    const std::vector<M> &iroots() const {
        return iw;
    }


};

template <class M>
void ntt(std::vector<M> &a, bool inverse = false) {
    static_assert(zeno::is_modular_v<M>, "Template parameter M must be modular (modint).");


    size_t n = a.size(), s = 0;
    if(n <= 1) return;
    while((1 << s) < n) s++;
    assert(n == (1 << s));



    NumberTheoreticTransform<M> ntt_info = NumberTheoreticTransform<M>::get_info();
    const std::vector<M> &root = ntt_info.roots(), &iroot = ntt_info.iroots();

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

    for (size_t l = 1; (size_t(1) << l) <= n; l++) {
        size_t len = size_t(1) << l;
        assert(l < root.size());
        M wlen = inverse ? iroot[l] : root[l];

        M wlenT = inverse ? ntt_info.get_g().pow(ntt_info.get_c()).inv() : ntt_info.get_g().pow(ntt_info.get_c());
        for(size_t i = len; i < (size_t(1) << ntt_info.get_ordlog()); i <<= 1) {
            wlenT = wlenT * wlenT;
        }
        assert(wlenT == wlen);

        assert(wlen.pow(len) == M(1));
        for(int i = 1; i < len; i++) 
            assert(wlen.pow(i) != M(1));


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


template <class M>
void inv_ntt(std::vector<M> &a) {
    static_assert(zeno::is_modular_v<M>, "Template parameter M must be modular (modint).");
    ntt(a, true);
}

template <class M>
std::vector<M> convolution_ntt(std::vector<M> const &a, std::vector<M> const &b) {
    static_assert(zeno::is_modular_v<M>, "Template parameter M must be modular (modint).");

    if(a.empty() || b.empty()) return {};

    size_t n = size_t(a.size()), m = size_t(b.size());
    size_t N = zeno::fft::compute_convolution_size(n, m);

    std::vector<M> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    fa.resize(N, M(0));
    fb.resize(N, M(0));

    ntt(fa);

    ntt(fb);

    for(int i = 0; i < N; i++) 
        fa[i] *= fb[i];

    inv_ntt(fa);

    fa.resize(n + m - 1);

    return fa;
}

} // namespace fft
    
} // namespace zeno


#endif /* ZENO_NTT_HPP */