/**
 * @file fft.hpp
 * @author Yuumita
 * @brief Fast Fourier Tranformation Interface
 * @version 0.1
 * @date 2023-03-15
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef ZENO_FFT_HPP
#define ZENO_FFT_HPP

#include "internal.hpp"
#include "internal_fft.hpp"
#include "modular.hpp"
#include <vector>
#include <math.h> // acos
#include <algorithm> // std::swap

namespace zeno {

namespace fft {

    size_t compute_convolution_size(size_t n, size_t m) {
        return (1u << internal::ceil_log2(n + m - 1));
    }


    using ftype = double;
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
        fft(a0, inverse);
        fft(a1, inverse);

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

    template<typename T>
    std::vector<T> convolution_fft_complex(std::vector<T> const &a, std::vector<T> const &b) {
        if(a.empty() || b.empty()) return {};

        size_t n = size_t(a.size()), m = size_t(b.size());
        size_t N = compute_convolution_size(n, m);

        std::vector<complex<T>> fa(N, complex<T>(0)), fb(N, complex<T>(0)); // maybe replace T with ftype ???
        for(int i = 0; i < n; i++) fa[i] = complex<T>(a[i]);
        for(int i = 0; i < m; i++) fb[i] = complex<T>(b[i]);

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

    template<typename T>
    void fft(std::vector<T> &a, bool inverse = false) {
        
    }

    template<typename T>
    void inv_fft(std::vector<T> &a) { fft(a, true); }

    const int magic_number = 0; // 60;


    const int cut =  (1 << 15);
    template<int m>
    std::vector<modint<m>> convolution_wcut(std::vector<modint<m>> a, std::vector<modint<m>> b) {
        size_t na = size_t(a.size()), nb = size_t(b.size());
        if(na == 0 || nb == 0) return {};
        if(na <= magic_number && nb <= magic_number) {
            if(na < nb) {
                std::swap(na, nb);
                std::swap(a, b);
            }
            std::vector<modint<m>> ans(na + nb - 1, modint<m>(0));
            for(int i = 0; i < na; i++) {
                for(int j = 0; j < nb; j++) {
                    ans[i + j] += a[i] * b[j];
                }
            }
            return ans;
        }

        size_t N = compute_convolution_size(na, nb);
        // Az(x) = A1(x) + iA2(x) = (A%cut)(x) + i(A/cut)(x)   [same for B(x)]
        std::vector<complex<ftype>> Az(N), Bz(N);
        for(int i = 0; i < na; i++) Az[i] = complex<ftype>(int(a[i]) % cut, int(a[i]) / cut);
        for(int i = 0; i < nb; i++) Bz[i] = complex<ftype>(int(b[i]) % cut, int(b[i]) / cut);


        fft(Az);
        fft(Bz);

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

        fft(P1, true); 
        fft(P2, true);

        std::vector<modint<m>> ans(N);
        for(int i = 0; i < N; i++) {
            modint<m> v0, v1, v2;

            // (A1 * B1)(x)
            v0 = llround(P1[i].Re());
            // (A1 * B2 + A2 * B1)(x)
            v1 = llround(P1[i].Im() + P2[i].Re());
            // (A2 * B2)(x)
            v2 = llround(P2[i].Im());


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

    template <typename T>
     std::vector<T> convolution_naive(std::vector<T> const &a, std::vector<T> const &b) {
        if (a.empty() || b.empty()) return {};
        size_t n = a.size(), m = b.size();
        std::vector<T> ret(n + m - 1);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) 
                ret[i + j] += a[i] * b[j];
        return ret;
    }

    template <typename T>
    std::vector<T> convolution(std::vector<T> const &a, std::vector<T> const &b) {
        if (std::min(a.size(), b.size()) <= magic_number) 
            return convolution_naive(a, b);
        if(zeno::is_modular_v<T>) {
            // ntt
        } 
        /* 
        else if(std::is_integral_v<T>) {
            std::vector<bigint> A(a), B(b), R;
            R = convolution_
        } 
        */
        return convolution_fft_complex<T>(a, b);
    }

} // namespace fft

} // namespace zeno

#endif