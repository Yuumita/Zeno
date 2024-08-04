#pragma once

#include <vector>
#include <iostream>
#include <assert.h>

#include "fps.hpp"
#include "internal.hpp"


namespace zeno {

/// @brief Implementation of common Set Power Series operations.
/// @ref https://www.algorithmas.org/articles/fast_subset_convolution.html
/// @tparam R A commutative ring.
template <class R>
class SetPowerSeries {
public:

    /// @brief Applies the zeta transform to f:2^[n]->R in O(n2^n).
    static void zeta(std::vector<R> &f) {
        size_t N = f.size();
        assert((N & (N - 1)) == 0); // power of two
        for(int j = 1; j < N; j <<= 1) {
            for(int S = 0; S < N; S++)  {
                if(S & j)  
                    f[S] += f[S ^ j];
            }
        }
    }

    /// @brief Applies the mobius transform to f:2^[n]->R in O(n2^n).
    static void mobius(std::vector<R> &f) {
        size_t N = f.size();
        assert((N & (N - 1)) == 0); // power of two
        for(int j = N >> 1; j >= 1; j >>= 1) {
            for(int S = N - 1; S >= 0; S--)  {
                if(S & j) 
                    f[S] -= f[S ^ j];
            }
        }
    }


    /// @return The OR convolution (union product) of f,g:2^[n]->R in O(n2^n).
    static std::vector<R> or_convolution(std::vector<R> f, std::vector<R> g) {
        size_t N = f.size();
        assert(N == g.size() && (N & (N-1)) == 0);

        zeta(f);
        zeta(g);
        for(int S = 0; S < N; S++) f[S] *= g[S];
        mobius(f);

        return f;
    }

    /// @return The subset convolution of f,g:2^[n]->R in O(n^2 2^n).
    static std::vector<R> subset_convolution(std::vector<R> f, std::vector<R> g) {
        size_t N = f.size();
        assert(N == g.size() && (N & (N-1)) == 0);
        size_t n = internal::ceil_log2(N);

        std::vector<int> pc(N);
        for(int S = 0; S < N; S++) 
            pc[S] = internal::popcount(S);

        std::vector<FormalPowerSeries<R>> fx, gx;
        for(int S = 0; S < N; S++) {
            fx.coef(pc[S]) = f[S];
            gx.coef(pc[S]) = g[S];
        }

        zeta(fx);
        zeta(gx);
        for(int S = 0; S < N; S++)
            fx[S] *= gx[S];
        mobius(fx);

        for(int S = 0; S < N; S++)
            f[S] = fx.coef(pc[S]);

        return f;
    }


    /// @brief Applies the (superset) zeta transform to f:2^[n]->R in O(n2^n).
    static void super_zeta(std::vector<R> &f) {
        size_t N = f.size();
        assert((N & (N - 1)) == 0); // power of two
        for(int j = 1; j < N; j <<= 1) {
            for(int S = 0; S < N; S++)  {
                if((S & j) == 0)  
                    f[S] += f[S ^ j];
            }
        }
    }

    /// @brief Applies the (superset) mobius transform to f:2^[n]->R in O(n2^n).
    static void super_mobius(std::vector<R> &f) {
        size_t N = f.size();
        assert((N & (N - 1)) == 0); // power of two
        for(int j = N >> 1; j >= 1; j >>= 1) {
            for(int S = N - 1; S >= 0; S--)  {
                if((S & j) == 0) 
                    f[S] -= f[S ^ j];
            }
        }
    }

    /// @return The AND convolution (intersection product) of f,g:2^[n]->R in O(n2^n).
    static std::vector<R> and_convolution(std::vector<R> f, std::vector<R> g) {
        size_t N = f.size();
        assert(N == g.size() && (N &(N-1)) == 0);

        super_zeta(f);
        super_zeta(g);
        for(int S = 0; S < N; S++) f[S] *= g[S];
        super_mobius(f);

        return f;
    }

    /// @brief Computes the Walsh–Hadamard transform (n-dimensional FFT in GF(2)) of f in O(n2^n).
    /// @param inverse If true apply the inverse WHT instead.
    static void walsh_hadamard_transform(std::vector<R> &f, bool inverse = false) {
        size_t N = f.size();
        assert((N &(N-1)) == 0);
        for(int j = 1; j < N; j <<= 1) {
            for(int S = 0; S < N; S++) {
                if((S & j) == 0) {
                    R u = f[S], v = f[S ^ j];
                    f[S]     = u + v;
                    f[S ^ j] = u - v;
                }

            }
        }

        if(inverse) {
            for(auto &e: f) 
                e /= R(N);
        }
    }


    /// @return The XOR convolution of f,g:2^[n]->R in O(n2^n).
    static std::vector<R> xor_convolution(std::vector<R> f, std::vector<R> g) {
        size_t N = f.size();
        assert(N == g.size() && (N &(N-1)) == 0);

        walsh_hadamard_transform(f);
        walsh_hadamard_transform(g);
        for(int S = 0; S < N; S++) f[S] *= g[S];
        walsh_hadamard_transform(f, true);

        return f;
    }

};

} // namespace zeno