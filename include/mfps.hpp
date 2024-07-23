#ifndef ZENO_MFPS_HPP
#define ZENO_MFPS_HPP

#include <stdlib.h>
#include <array>
#include <algorithm>
#include "fps.hpp"

namespace zeno {

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


    size_t get_index(std::array<int, N> I) {
        size_t i = 0;
        for(size_t j = 0, p = 1; j < N; j++) {
            i += I[j] * p;
            p *= n[j];
        }
        return i;
    }

    static std::array<int, N> merge_radix(std::array<int, N> const &A, std::array<int, N> const &B) {
        std::array<int, N> ret;
        for(size_t i = 0; i < N; i++) ret[i] = std::max(A[i], B[i]);
        return ret;
    }

    /// @brief Replace the quotient ideal with <x_0^{m[0]}, ..., x_{N-1}^{m[N-1]}>
    mfps convert_radix(std::array<int, N> m) const {
        /// TODO:  ...

        fps T;
        for(size_t i = 0; i < N; i++) m[i] = std::min(m[i], this->n[i]);

        std::array<int, N> nn = {0};
        for(size_t i = 0, j = 0; i < F.size(); i++) {

            T[j] = F[i];

            bool carry = true;
            for(size_t k = 0; k < N && carry; k++) {
                ++nn[k], carry = false;
                if(nn[k] >= m[k]) 
                    nn[k] = 0, carry = true;
            }

            if(carry) break;

        }

    }

public:
    DenseMultivariateFormalPowerSeries();
    DenseMultivariateFormalPowerSeries(std::vector<int> const &n_)
        : n(n_) {}
    DenseMultivariateFormalPowerSeries(fps const f_)
        : f(f_), n(n_) {}

    
    template<typename... Args>
    R& operator()(Args... args) {
        static_assert(sizeof...(args) == N, "Number of arguments must be equal to N");
        std::array<int, N> indices = {args...};
        size_t index = get_index(indices);
        return F.coef(index);
    }


    mfps& operator+=(mfps const &rhs) {
        mfps &T = rhs;
        if(this->n != rhs->n)
            mfps T = rhs.convert_radix(merge_radix(this->n, rhs->n));
        for(int i = 0; i < T->F.size(); i++)
            this->F[i] += T->F[i];
        return *this;
    }

    mfps& operator-=(mfps const &rhs) {
        mfps &T = rhs;
        if(this->n != rhs->n)
            mfps T = rhs.convert_radix(merge_radix(this->n, rhs->n));
        for(int i = 0; i < T->F.size(); i++)
            this->F[i] -= T->F[i];
        return *this;
    }

    mfps& operator*=(mfps const &rhs) {
        /// TODO: ...
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

};


template<class R, size_t N>
using dmfps = DenseMultivariateFormalPowerSeries;
    
    
} // namespace zeno


#endif /* ZENO_MFPS_HPP */