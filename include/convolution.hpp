#pragma once

#include "internal.hpp"
#include "internal_fft.hpp"
#include "modular.hpp"
#include "ntt.hpp"
#include <vector>
#include <math.h> // acos
#include <algorithm> // std::swap

namespace zeno {

namespace fft {

    using ftype = long double;
    const ftype PI = acos(-1); 

    /// @brief Transform the sequence a to its Fourier transformation. 
    ///        If inverse is true then apply the inverse Fourier transform.
    /// @pre n is a power of 2
    template<typename T>
    void fft_complex(std::vector<complex<T>> &a, bool inverse = false) {
        size_t n = a.size();
        if(n <= 1) return;
        std::vector<complex<T>> a0(n/2), a1(n/2);
        for(int i = 0; i < n / 2; i++) {
            a0[i] = a[2*i];
            a1[i] = a[2*i + 1];
        }
        fft_complex(a0, inverse);
        fft_complex(a1, inverse);

        T angle = 2 * PI / n * (inverse ? -1 : +1);
        complex<T> wn(cos(angle), sin(angle)), w(1);
        for(int i = 0; i < n/2; i++) {
            a[i] = a0[i] + w * a1[i];
            a[i + n/2] = a0[i] - w * a1[i];
            if(inverse) {
                a[i] /= T(2);
                a[i + n/2] /= T(2);
            }
            w *= wn;
        }
    }

    template<typename T>
    void inv_fft_complex(std::vector<T> &a) { fft_complex(a, true); }

    template<typename T = long double>
    std::vector<T> convolution_fft_complex(std::vector<T> const &a, std::vector<T> const &b) {
        if(a.empty() || b.empty()) return {};

        size_t n = size_t(a.size()), m = size_t(b.size());
        size_t N = compute_convolution_size(n, m);

        std::vector<complex<ftype>> fa(N, complex<ftype>(0)), fb(N, complex<ftype>(0)); // maybe replace T with ftype ???
        for(int i = 0; i < n; i++) fa[i] = complex<ftype>(a[i]);
        for(int i = 0; i < m; i++) fb[i] = complex<ftype>(b[i]);

        fft_complex(fa);
        fft_complex(fb);

        for(int i = 0; i < N; i++)
            fa[i] *= fb[i];

        inv_fft_complex(fa);
        
        std::vector<T> ans(n + m - 1);
        for(size_t i = 0; i < n + m - 1; i++) 
            ans[i] = T(fa[i].Re());

        return ans;
    }

    /// @tparam F A field
    template<typename F>
    void fft(std::vector<F> &a, std::vector<F> &root, std::vector<F> &iroot, bool inverse = false) {
        size_t n = a.size(), s = 0;
        if(n <= 1) return;
        while((1 << s) < n) s++;
        assert(n == (1 << s));

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
            F wlen = inverse ? iroot[l] : root[l];
            for (size_t i = 0; i < n; i += len) {
                F w(1);
                for(size_t j = 0; j < len / 2; j++) {
                    F u = a[i + j];
                    F v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }

        if(inverse) {
            for(F &e: a)
                e /= F(n);
        }

    }

    template<typename F>
    inline void inv_fft(std::vector<F> &a, std::vector<F> &root, std::vector<F> &iroot) { 
        fft(a, true); 
    }


    template<typename T>
    std::vector<T> convolution_wcut(std::vector<T> const &a, std::vector<T> const &b) {
        if(a.empty() || b.empty()) return {};

        size_t na = size_t(a.size()), nb = size_t(b.size());
        size_t N = compute_convolution_size(na, nb);


        T cut = zeno::internal::pow(T(2), sizeof(T) * 4u);

        // Az(x) = A1(x) + iA2(x) = (A%cut)(x) + i(A/cut)(x)   [same for B(x)]
        std::vector<complex<ftype>> Az(N), Bz(N);
        for(int i = 0; i < na; i++) {
            T x0 = std::floor(a[i] / cut), x1 = a[i] - x0;
            Az[i] = complex<ftype>(x0, x1);

        }

        for(int i = 0; i < nb; i++) {
            T x0 = std::floor(b[i] / cut), x1 = b[i] - x0;
            Bz[i] = complex<ftype>(x0, x1);
        }


        fft_complex(Az);
        fft_complex(Bz);

        // P1(x) = A1(x) * Bz(x), P2(x) = A2(x) * Bz(x)
        std::vector<complex<ftype>> P1(N), P2(N);

        // We have to compute P1, P2 using element-wise operations of DFT(A) and DFT(B)
        // Re{DFT(Az)[i]} = (Az[i] + conj(Az)[i])/2  = (Az[i] + conj(Az[n-i]))/2
        // Im{DFT(Az)[i]} = (Az[i] - conj(Az)[i])/2i = (Az[i] - conj(Az[n-i]))/2i
        for(int i = 0; i < N; i++) {
            int j = (i ? N - i : 0);
            P1[i] = (Az[i] + Az[j].conj()) * complex<ftype>(0.5, 0) * Bz[i];
            P2[i] = (Az[i] - Az[j].conj()) * complex<ftype>(0, -0.5) * Bz[i];
        }

        inv_fft_complex(P1); 
        inv_fft_complex(P2);

        std::vector<T> ans(N);
        for(int i = 0; i < N; i++) {
            T v0, v1, v2;

            // (A1 * B1)(x)
            v0 = T(P1[i].Re());
            // (A1 * B2 + A2 * B1)(x)
            v1 = T(P1[i].Im() + P2[i].Re());
            // (A2 * B2)(x)
            v2 = T(P2[i].Im());

            ans[i] = v0 + v1*cut + v2*cut*cut;
        }

        return ans;
    }



    template<typename T>
    std::vector<T> convolution_fft(std::vector<T> const &a, std::vector<T> const &b) {
        if(a.empty() || b.empty()) return {};

        size_t n = size_t(a.size()), m = size_t(b.size());
        size_t N = compute_convolution_size(n, m);

        std::vector<T> fa(a.begin(), a.end()), fb(b.begin(), b.end());

        fft(fa);
        fft(fb);
        for(int i = 0; i < N; i++) 
            fa[i] *= fb[i];
        inv_fft(fa);

        fa.resize(n + m - 1);
        return fa;
    }





    /// @ref https://core.ac.uk/download/pdf/82485770.pdf, https://www.luogu.com/article/wje8kchr, https://rushcheyo.blog.uoj.ac/blog/6547
    /// @brief Computes the multivariate convolution G(x1, ..., xN) H(x1, ..., xN) (mod <x1^n1, ..., xN^nN>).
    /// @pre [x^I]A (as multivariate) -> [x^i]B (as univariate) where I = (i1, ..., iN) is in mixed-radix wrt n and i = i1 + n1(i2 + n2(...)).
    /// @return The resulting FPS (in mixed-radix wrt n)
    template <typename FPS, size_t N>
    FPS multivariate_convolution(FPS const &A, FPS const &B, std::array<int, N> n) {
        if(N == 0) return FPS(A[0] * B[0]);

        size_t l = std::max(A.size(), B.size());
        size_t L = compute_convolution_size(l, l);

        std::vector<int> chi(l, 0);
        for(size_t i = 0; i < l; ++i) {
            for(size_t j = 0, tmp = i / n[j]; j < N; ++j, tmp /= n[j])
                chi[i] += tmp;
            // chi[i] %= N
        }

        /// TODO: ...

    }




} // namespace fft


template <typename T>
    std::vector<T> convolution_naive(std::vector<T> const &a, std::vector<T> const &b) {
    if (a.empty() || b.empty()) return {};
    size_t n = a.size(), m = b.size();
    std::vector<T> ret(n + m - 1);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j) 
            ret[i + j] += static_cast<T>(a[i]) * static_cast<T>(b[j]);
    return ret;
}


template <typename T>
std::vector<T> convolution(std::vector<T> const &a, std::vector<T> const &b) {
    if (std::min(a.size(), b.size()) <= fft::magic_number) 
        return convolution_naive(a, b);
    if constexpr (zeno::is_modular_v<T>) {
        return fft::convolution_ntt<T>(a, b);
    } 
    /* 
    else if(std::is_integral_v<T>) {
        // strassen ??? 
        std::vector<bigint> A(a), B(b), R;
        R = convolution_
    } 
    */
    return fft::convolution_fft_complex(a, b);
}

} // namespace zeno