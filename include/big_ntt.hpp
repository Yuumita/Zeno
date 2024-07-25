#pragma once

#include <vector>
#include "modular.hpp"
#include "ntt.hpp"

namespace zeno
{
    

/// @ref nyaannyaan, https://cp-algorithms.com/algebra/garners-algorithm.html
class BigNTT {
    using int128_t = __int128_t;
private:
    /// @brief Garner's algorithm
    static const int m0 = 167772161;
    static const int m1 = 469762049;
    static const int m2 = 754974721;
    using mint0 = modint<m0>;
    using mint1 = modint<m1>;
    using mint2 = modint<m2>;
    static int r01; // = mint1(m0).inv().val();
    static int r02; // = mint2(m0).inv().val();
    static int r12; // = mint2(m1).inv().val();
public:

    template <typename T, class M>
    static std::vector<M> _convolution_mint(std::vector<T> const &a, std::vector<T> const &b) {
        std::vector<M> A(a.size()), B(b.size());
        for(int i = 0; i < a.size(); i++) A[i] = M(a[i]);
        for(int i = 0; i < b.size(); i++) B[i] = M(b[i]);
        return fft::convolution_ntt(A, B);
    }

    template <typename T>
    static std::vector<int128_t> convolution_int128(std::vector<T> const &a, std::vector<T> const &b) {
        if(a.empty() || b.empty()) return {};
        if (std::min(a.size(), b.size()) <= fft::magic_number) 
            return convolution_naive<int128>(a, b);
        
        std::vector<mint0> x0 = _convolution_mint<T, mint0>(a, b);
        std::vector<mint1> x1 = _convolution_mint<T, mint1>(a, b);
        std::vector<mint2> x2 = _convolution_mint<T, mint2>(a, b);
        // xi[j] = (a * b)[j] (mod mi)
        // retrieve x = (a * b)[j] (mod m0 * m1 * m2) through chinese remainder theorem (garner's algorithm)

        std::vector<int128> ret(x0.size());
        for(int i = 0; i < x0.size(); i++) {
            /// int is probably enough for a0, a1, a2
            int64_t a0 = x0[i].val();
            int64_t a1 = (int64_t)(x1[i] + m1 - a0) * r01 % m1;
            int64_t a2 = ((int64_t)(x2[i] + m2 - a0) * r02 * r12 % m2 + (int64_t)(m2 - a1) * r12) % m2;

            /// TODO: check that the following doesn't contain danger of overflow
            ret[i] = a0 + (int64_t)(a1) * m0 + (int128_t)(a2) * m0 * m1; 
        }

        return ret;
    }
};

int BigNTT::r01 = modint<m1>(m0).inv().val();
int BigNTT::r02 = modint<m2>(m0).inv().val();
int BigNTT::r12 = modint<m2>(m1).inv().val();


} // namespace zeno