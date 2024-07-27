#pragma once

#include <stdlib.h>
#include <array>
#include <algorithm>
#include "fps.hpp"

namespace zeno {

/// TODO: Complete the class
/// @brief Implementation of the ring R[[x_0, ...,x_{N-1}]] / <x_0^n[0], ..., x_{N-1}^n[N-1]>.
/// @tparam R The ring R.
/// @tparam N the number of indetermined variables.
template<class R, size_t N>
class DenseMultivariateFormalPowerSeries {
    using mfps = DenseMultivariateFormalPowerSeries<R, N>;
    using fps = FormalPowerSeries<R>;
private:

    fps F;
    std::array<int, N> n;

    /// [x^i]F = [x^{(i_0, i_1, ..., i_{N-1})}]F where  
    /// (i_0, ..., i_{N-1}) is the mixed radix of i w/r to (n[0], ..., n[N-1])
    /// @ref https://core.ac.uk/download/pdf/82485770.pdf, https://rushcheyo.blog.uoj.ac/blog/6547

    static size_t get_index(std::array<int, N> const &I, std::array<int, N> const &base) {
        size_t i = 0;
        for(size_t j = 0, p = 1; j < N; j++) {
            i += I[j] * p;
            p *= base[j];
        }
        return i;
    }

    size_t get_index(std::array<int, N> const &I) {
        return get_index(I, n);
    }

    static std::array<int, N> get_digits(size_t nu, std::array<int, N> const &base) {
        std::array<int, N> d;
        for(size_t i = 0; i < N; i++) {
            d[i] = nu % base[i];
            nu = (nu - nu % base[i]) / base[i];
        }
        return d;
    }

    static std::array<int, N> merge_radix(std::array<int, N> const &A, std::array<int, N> const &B) {
        std::array<int, N> ret;
        for(size_t i = 0; i < N; i++) ret[i] = std::max(A[i], B[i]);
        return ret;
    }

    /// @brief Replace the quotient ideal with <x_0^{m[0]}, ..., x_{N-1}^{m[N-1]}>
    mfps convert_to_radix(std::array<int, N> m) const {
        size_t nproduct = 1;
        for(size_t i = 0; i < N; i++) {
            m[i] = std::min(m[i], this->n[i]);
            nproduct *= n[i];
        }

        fps &Pn = this->F, Pm(0);
        for(size_t nu = 0; nu < nproduct; nu++) {
            std::array<int, N> nud = get_digits(nu, n);
            bool ok = true;
            for(size_t i = 0; i < N && ok; i++)
                if(nud[i] >= m[i]) ok = false;
            if(!ok) break;
            size_t mu = get_index(nud, m);
            Pm[mu] = Pn[nu];
        }
        return mfps(m, Pm);
    }


    /// TODO: rename the following
    /// @brief Make G the same radix as H
    static void make_same_radix(mfps &G, mfps &H) {
        if(G.n != H.n) {
            std::array<int, N> m = merge_radix(G.n, H.n);
            G = G.convert_to_radix(m);
            H = H.convert_to_radix(m);
        }
    }


public:
    DenseMultivariateFormalPowerSeries();
    DenseMultivariateFormalPowerSeries(std::array<int, N> const &n_)
        : n(n_) {}
    DenseMultivariateFormalPowerSeries(std::array<int, N> const &n_, fps const &F_)
        : n(n_), F(F_){}
    
    template<typename... Args>
    R& operator()(Args... args) {
        static_assert(sizeof...(args) == N, "Number of arguments must be equal to N");
        std::array<int, N> indices = {args...};
        size_t index = get_index(indices);
        return F.coef(index);
    }

    mfps& operator+=(mfps const &rhs) {
        if(this->n != rhs.n) {
            assert(false);
            mfps T = rhs;
            make_same_radix(*this, T);
            this->F += T.F;
        } else {
            this->F += rhs.F;
        }
        return *this;
    }

    mfps& operator-=(mfps const &rhs) {
        if(this->n != rhs.n) {
            assert(false);
            mfps T = rhs;
            make_same_radix(*this, T);
            this->F -= T.F;
        } else {
            this->F -= rhs.F;
        }
        return *this;
    }

    mfps& operator*=(mfps const &rhs) {
        if(this->n != rhs.n) {
            assert(false);
            mfps T = rhs;
            make_same_radix(*this, T);
            this->F = fft::multivariate_convolution<FPS, N>(this->F, T.F, this->n);
        } else {
            this->F = fft::multivariate_convolution<FPS< N>(this->F, rhs.F, this->n);
        }
        return *this;
    }

    mfps& operator/=(mfps const &rhs) {
        return (*this) *= rhs.inv();
    }

    mfps& operator+=(R const &rhs) {
        F.coef(0) += rhs;
        return *this;
    }
    mfps& operator-=(R const &rhs) {
        F.coef(0) -= rhs;
        return *this;
    }
    mfps& operator*=(R const &rhs) {
        for(R &e: F) e *= rhs;
        return *this;
    }
    mfps& operator/=(R const &rhs) {
        for(R &e: F) e /= rhs;
        return *this;
    }


    mfps operator+(mfps const &rhs) const { return mfps(*this) *= rhs; }
    mfps operator*(mfps const &rhs) const { return mfps(*this) *= rhs; }
    mfps operator-(mfps const &rhs) const { return mfps(*this) -= rhs; }
    mfps operator/(mfps const &rhs) const { return mfps(*this) /= rhs; }

    mfps operator+(R const &rhs) const { return mfps(*this) *= rhs; }
    mfps operator*(R const &rhs) const { return mfps(*this) *= rhs; }
    mfps operator-(R const &rhs) const { return mfps(*this) -= rhs; }
    mfps operator/(R const &rhs) const { return mfps(*this) /= rhs; }

    mfps operator+() const { return mfps(*this); }
    mfps operator-() const { return mfps(-F, base); }

    // friend bool operator==(const mfps& lhs, const mfps& rhs) {
    //     return lhs.f == rhs.f && lhs.base == rhs.base;
    // }
    // friend bool operator!=(const mfps& lhs, const mfps& rhs) {
    //     return lhs.f != rhs.f || lhs.base != rhs.base;
    // }



};


template<class R, size_t N>
using dmfps = DenseMultivariateFormalPowerSeries;
    
    
} // namespace zeno