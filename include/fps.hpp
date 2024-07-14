/**
 * @file fps.hpp
 * @author Yuumita
 * @brief Formal Power Series Interface
 * @version 0.1
 * @date 2023-03-15
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef ZENO_POLY_HPP
#define ZENO_POLY_HPP

#include <vector>
#include <iostream>
#include <assert.h>
#include <algorithm> // std::min

#include "fft.hpp"
#include "internal.hpp"


namespace zeno {

/// @brief The Polynomial class
/// @tparam R The type of coefficients
template<typename R>
class FormalPowerSeries {
    using FPS = FormalPowerSeries;
private:
    /// @brief Coefficients of the polynomial
    std::vector<R> a;
public:
    /// @brief Reduce the 0 coefficients of the Polynomial (i.e remove all trailing zeroes)
    void normalize(){
        while(!this->a.empty() && this->a.back() == R(0))
            this->a.pop_back();
    }
    /// @brief Construct a new Polynomial object with degree -1 (No coefficients) 
    FPS(): a({}) {}

    /// @brief Construct a new Polynomial object with coefficients from the list coefs 
    /// @param coefs The coefficient list
    /// @param n The size of the list
    FPS(const int n, const R *coefs) {
        this->a.resize(n);
        for(int i = 0; i < n; i++) a[i] = coefs[i];
        this->normalize();
    }

    /// @brief Construct a new Polynomial object with coefficients from the array coefs 
    ///  
    /// @param coefs The coefficient array
    FPS(const std::vector<R> &coefs): a(coefs) {
        this->normalize();
    }

    /// @brief Construct a new (constant, i.e of 0 degree) Polynomial object which is just v
    /// @param v The value of the constant polynomial
    FPS(const R &v) {
        this->a.assign(1, v);
        this->normalize();
    }

    /// @brief Construct a new Polynomial object of degree n-1 with initial values v 
    /// @param n The coefficient count
    /// @param v The value to initialize the coefficients
    FPS(const int n, const R &v) {
        this->a.resize(n);
        for(int i = 0; i < n; i++) a[i] = v;
    }

    /// @brief Return the degree of the polynomial (-1 if there aren't any coefficients)
    /// @return int The degree of the polynomial
    int deg() const {
        return this->a.size() - 1;
    }

    int ord() const {
        for(int i = 0; i < a.size(); i++) 
            if(a[i] != 0) return i;
        return -1;
    }

    int is_zero() const {
        return this->a.empty();
    }

    /// @brief Get the coefficient of x^i
    /// @param i The index of the coefficient
    /// @return R The coefficient
    R get(const int i) const {
        return (0 <= i && i <= deg() ? this->a[i] : R(0));
    }

    /// @brief Set the coefficient of x^i
    /// @param i The index of the coefficient
    /// @param c The value of the new coefficient
    void set(const int i, const R c) {
        assert(0 <= i);
        if(this->a.size() <= i) this->a.resize(i+1);
        this->a[i] = c;
    }

    void resize(const int i) {
        if(deg() < i) this->a.resize(i + 1);
    }

    R leadcoef() const {
        if(this->a.size() == 0) return 0;
        return this->a.back();
    }


    /// @brief Addition of two polynomials (*this and P)
    /// @param P The other polynomial
    /// @return FPS& The resulting polynomial
    FPS plus(const FPS& P) const {
        return *this + P;
    }


    /// @brief Multiplication of two polynomials (*this and P)
    /// @return FPS& The resulting polynomial
    FPS times(const FPS& P) const {
        return *this * P;
    }

    /// @brief Multiplication of this polynomial with the scalar v
    /// @return FPS The resulting polynomial
    FPS times(const R& v) const {
        return *this * v;
    }

    /// @brief Subtraction of two polynomials (*this - P)
    /// @return FPS The resulting polynomial
    FPS minus(const FPS& P) const {
        return *this - P;
    }

    FPS& operator+=(const FPS &rhs) {
        int d = (this->deg() > rhs.deg() ? this->deg() : rhs.deg());
        this->a.resize(d);
        for(int i = 0; i < (int)rhs.a.size(); i++){
            this->a[i] += rhs.a[i];
        }
        this->normalize();
        return *this;
    }

    FPS& operator-=(const FPS &rhs) {
        int d = (this->deg() > rhs.deg() ? this->deg() : rhs.deg());
        this->a.resize(d + 1);
        for(int i = 0; i < (int)rhs.a.size(); i++){
            this->a[i] -= rhs.a[i];
        }
        this->normalize();
        return *this;
    }

    FPS& operator*=(const FPS &rhs) {
        this->a = fft::convolution<R>(this->a, rhs.a);
        this->normalize();
        return *this;
    }

    FPS& operator*=(const R &rhs) {
        for(auto &e: this->a) e *= rhs;
        this->normalize();
        return *this;
    }

    FPS& operator/=(const FPS &rhs) {
        return *this = div(rhs).first;
    }

    FPS& operator/=(const R &rhs) {
        for(auto &e: this->a) e /= rhs;
        return *this;
    }

    FPS& operator%=(const FPS &rhs) {
        return *this = div(rhs).second;
    }

    // TODO: >>=, <<= shifts
    // FPS& operator<<=(const int s) {
    //     return *this = (*this).shift(s);
    // }
    // FPS& operator>>=(const int s) {
    //      return operator<<=(-s);
    // }

    FPS operator-() const {
        auto t = *this;
        for(auto &e: t.a) e = -e;
        return t;
    }

    FPS operator+(const FPS &rhs) const { return FPS(*this) += rhs; }
    FPS operator-(const FPS &rhs) const { return FPS(*this) -= rhs; }
    FPS operator*(const FPS &rhs) const { return FPS(*this) *= rhs; }
    FPS operator*(const R &rhs)    const { return FPS(*this) *= rhs; }
    FPS operator/(const FPS &rhs) const { return FPS(*this) /= rhs; }
    FPS operator/(const R &rhs)    const { return FPS(*this) /= rhs; }
    FPS operator%(const FPS &rhs) const { return FPS(*this) %= rhs; }

    // mutable reference of the coefficient
    R& coef(size_t index) { return this->a[index]; }

    bool operator == (const FPS &rhs) const { return a == rhs.a; }
    bool operator != (const FPS &rhs) const { return a != rhs.a; }

    /// @return FPS A copy of this polynomial mod x^k (getting the first k coefficients)
    FPS mod_xk(int k) const {
        return FPS(std::vector<R>(a.begin(), a.begin() + std::min(k, (int)a.size())));
    }
    inline static FPS mod_xk(FPS &f, int k) const { return f.mod_xk(k); }

    /// @return FPS A copy of this polynomial times x^k (coefficients shifted k steps to the right)
    FPS times_xk(int k) const {
        std::vector<int> coefs = this->a;
        coefs.insert(coefs.begin(), k, 0);
        return FPS(coefs);
    }
    inline static FPS times_xk(FPS &f, int k) const { return f.times_xk(k); }

    /// @return FPS A copy of this polynomial with coefficients shifted k steps to the left
    FPS floordiv_xk(int k) const {
        return FPS(std::vector<R>(a.begin() + k, a.end()));
    }
    inline static FPS floordiv_xk(FPS &f, int k) const { return f.floordiv_xk(k); }

    FPS shift(int k) const {
        if(k > 0)
            return times_xk(k);
        if(k < 0)
            return floordiv_xk(-k);
        return FPS(*this);
    }
    inline static FPS shift(FPS &f, int k) const { return f.shift(k); }


    /// @return FPS P(x^{1/2})
    FPS sqrt_x() const {
        std::vector<R> ret(this->a.size(), R(0));
        for(int i = 0; i < this->a.size(); i++)
            ret[i] = get(2*i);
        return FPS(ret);
    }
    inline static FPS sqrt_x(FPS &f) const { return f.sqrt_x(); }

    /// @return FPS P(x^2)
    FPS sqr_x() const {
        std::vector<R> ret(2 * this->a.size(), R(0));
        for(int i = 0; i < this->a.size(); i++)
            ret[2*i] = this->a[i];
        return FPS(ret);
    }
    inline static FPS sqr_x(FPS &f) const { return f.sqr_x(); }

    /// @return FPS The multiplicative inverse of this polynomial mod x^m
    FPS inv(int m = deg() + 1) const {
        assert(get(0) != R(0));
        if(m == 0) return FPS(0);
        FPS Q(R(1) / get(0));
        for(int i = 1; i < m; i *= 2)
            Q = (Q * R(2) - Q * Q * this->mod_xk(2*i)).mod_xk(2*i);
        return Q.mod_xk(m);
        /*
            if(m == 1) return FPS(R(1) / get(0)); 
            FPS A = this->mod_xk(m);
            FPS Am = A.scale_x(R(-1));
            FPS B = Am * A;
            B = B.sqrt_x().inv((m+1) / 2).sqr_x();
            std::cerr << "inv (m = " << m << ") --> " << (Am * B).mod_xk(m) << std::endl;
            std::cerr << A.mod_xk(m) << " * " << (Am*B).mod_xk(m) << " = " << (A*Am*B).mod_xk(m) << std::endl;
            std::cerr << "inv (m = " << m << ") --> " << (Am * B).mod_xk(m) << std::endl;
            return (Am * B).mod_xk(m);
        */
    }
    inline static FPS inv(FPS &f, int m = deg() + 1) const { return f.inv(m); }

    /// @return FPS  A copy of this polynomial with coefficients reversed
    FPS reverse() const {
        return FPS(std::vector<R>(a.rbegin(), a.rend()));
    }
    inline static FPS reverse(FPS &f) const { return f.reverse(Q); }

    /// @brief Euclidean division of this polynomial with the polynomial P
    /// @return pair<FPS, FPS> the quotiend (first) and remainder (second) polynomials
    std::pair<FPS, FPS> div(const FPS &P) const {
        assert(!P.is_zero());
        if(P.deg() > this->deg())
            return {FPS(0), FPS(*this)};
        FPS A_R = this->reverse(), B_R = P.reverse();
        int d = this->deg() - P.deg();
        FPS D = ((A_R).mod_xk(d + 1) * P.inv(d + 1)).mod_xk(d + 1).reverse();
        return {D, *this - D * P};
    }
    inline static std::pair<FPS, FPS> div(FPS &f, const FPS &Q) const { return f.div(Q); }

    /// @brief component-wise multiplication with v^k
    /// @return FPS A copy of the polynomial with the x value scaled by v
    FPS scale_x(const R &v) const {
        R c = 1;
        FPS P(*this);
        for(int i = 0; i <= P.deg(); i++, c *= v) { 
            P.coef(i) *= c;
        }
        return  P;
    }
    inline static FPS scale_x(FPS &f, const R &v) const { return f.scale_x(v); }



    // Computes the derivative
    FPS deriv(int k = 1) const {
        if(deg() + 1 < k) return FPS(R(0));
        std::vector<R> res(deg() - k +  1, R(0));
        if(k == 1) {
            for(int i = 1; i <= deg(); i++)
                res[i - 1] = get(i) * R(i);
        } else {
            for(int i = k; i <= deg(); i++)
                res[i - k] = get(i) * internal::fact<R>(i) / internal::fact<R>(i - k);
        }
        return FPS(res);
    }

    inline static FPS deriv(FPS &f, int k = 1) const { return f.deriv(k); }


    /// @return The integral of the polynomial (with constant 0)
    FPS integr() const {
        std::vector<R> res(deg() + 2);
        for(int i = 1; i <= deg() + 1; i++) {
            res[i] = get(i-1) / i;
        }
        return FPS(res);
    }
    inline static FPS integr(FPS &f) const { return f.integr(k); }

    /// @return The (natural) logarithm of the polynomial (mod x^m)
    FPS log(int m) const { return (FPS(*this).deriv().mod_xk(m) * FPS(*this).inv(m)).mod_xk(m);  }

    inline static FPS log(FPS &f, int m) const { return f.log(m); }

    FPS exp(int m) const {
        assert(get(0) == 0);
        std::cerr << "T0" << std::endl;
        FPS Q(1), P(*this);
        P.coef(0) += 1;
        for(int i = 1; i < m; i *= 2) { // at the end of each step: F = e^P (mod x^{2^{i+1}})
            std::cerr << "T0.5" << std::endl;
            // Q = (Q * (P - Q.log(2*i))).mod_xk(2*i);
            Q = (Q * (P - Q.mod_xk(2*i))).mod_xk(2*i);
            std::cerr << "T0.8" << std::endl;
        }
        std::cerr << "T1" << std::endl;
        return Q.mod_xk(m);
    }

    inline static FPS exp(FPS &f, int m) const { return f.exp(m); }

    FPS sqrt(int k = 2) const {

    }

    inline static FPS sqrt(FPS &f, int m) const { return f.sqrt(m); }

    /// @brief Evaluate the polynomial at point x0
    /// @return The value of the polynomial at point x0 
    R eval(R x0) const {
        R y(0);
        for(int i = this->deg(); i >= 0; i--) {
            y = (x0 * y) + a[i];
        }
        return y;
    }


    /// @brief Builds the polynomial evaluation tree.
    /// @return The polynomial P(x) = (x - x_L) (x - x_{L+1}) ... (x - x_R)
    static FPS<R> build_poly_tree(std::vector<FPS<R>> &tree, 
        int v, typename std::vector<R>::iterator l, typename std::vector<R>::iterator r) {
            if(r - l == 1) {
                return tree[v] = FPS({-*l, 1});
            } else {
                auto m = l + (r - l) / 2;
                return tree[v] = build_poly_tree(tree, 2*v, l, m) 
                    * build_poly_tree(tree, 2*v+1, m, r);
            }
        }

    std::vector<R> eval(std::vector<FPS<R>> &tree, 
        int v, typename std::vector<R>::iterator l, typename std::vector<R>::iterator r) {
            if(r - l == 1) {
                return {eval(*l)};
            } else {
                auto m = l + (r - l) / 2;
                std::vector<R> A = (*this % tree[2*v]).eval(tree, 2*v, l, m);
                std::vector<R> B = (*this % tree[2*v+1]).eval(tree, 2*v+1, m, r);
                A.insert(end(A), begin(B), end(B));
                return A;
            }
        }

    /// @brief Evaluate the polynomial in points x0, x1, ..., x_{n-1}
    /// @return The value of the polynomial at points x0, x1, ..., x_{n-1}
    std::vector<R> eval(std::vector<R> const &x) const {
        int n = x.size();
        if(is_zero()) return std::vector<R>(n, R(0));
        std::vector<FPS<R>> tree(4*n);
        build_poly_tree(tree, 1, begin(x), end(x));
        return eval(tree, 1, begin(x), end(x));
    }


    static FPS<R> _interpolate_tree(std::vector<FPS<R>> const &tree, 
        int v, typename std::vector<R>::iterator l, typename std::vector<R>::iterator r) {
            if(r - l == 1) {
                return tree[v] = FPS<R>(*l);
            } else {
                auto m = l + (r - l) / 2;
                FPS<R> A = _interpolate_tree(tree, 2*v, l, m);
                FPS<R> B = _interpolate_tree(tree, 2*v+1, m, r);
                return tree[v] = A * tree[2*v+1] + tree[2*v] * B;
            }
    }


    static FPS<R> interpolate(std::vector<R> const &x, std::vector<R> const &y) {
        assert(x.size() == y.size());
        std::vector<FPS<R>> tree(4*x.size());
        std::vector<R> u = FPS<R>::build_poly_tree(tree, 1, begin(x), end(x)).deriv().eval(x);
        for(int i = 0; i < y.size(); i++)
            u[i] = y[i] / u[i];
        return _interpolate_tree(tree, 1, begin(u), end(u));
    }


    /// @return The Borel transformation of the polynomial (a_k x^k -> a_k / k! x^k)
    FPS<R> borel() const {
        FPS<R> ret = *this;
        for(int k = 0; k <= deg(); k++)
            ret->a[k] /= internal::fact(k);
        return ret;
    }

    /// @return The inverse Borel (Laplace) transformation of the polynomial (a_k x^k -> k!a_k x^k)
    FPS<R> invborel() const {
        FPS<R> ret = *this;
        for(int k = 0; k <= deg(); k++)
            ret->a[k] *= internal::fact(k);
        return ret;
    }

};

template<class R>
using FPS = FormalPowerSeries;
template<class R>
using Poly = FormalPowerSeries;

template<typename R>
FPS<R> operator*(const R &lhs, const FPS<R> &rhs) {
    return FPS(rhs) *= lhs;
}

template<typename R>
FPS<R> operator/(const R &lhs, const FPS<R> &rhs) {
    return FPS(lhs) /= rhs;
}


template<typename R>
std::ostream& operator<<(std::ostream& os, const FPS<R> &P) {
    os << "[ ";
    for(int i = 0; i <= P.deg(); i++)
        os << P.get(i) << " ";
    os << "] ";
    return os;
}

namespace polynomial
{

template<class K>
void euclid_division(FPS<K> const &A, FPS<K> const &B, FPS<K> &Q, FPS<K> &R) {
    R = A, Q = FPS<K>(0);
    while(R.deg() >= B.deg()) {
        K s = R.leadcoef() / B.leadcoef();
        int sd = R.deg() - B.deg();

        Q.set(sd, Q.get(sd) + s);
        R = R - s * B.shift(sd);
    }
}

template<class K>
std::pair<FPS<K>, FPS<K>> euclid_division(FPS<K> const &A, FPS<K> const &B) {
    FPS<K> Q, R;
    polynomial::euclid_division(A, B, Q, R);
    return {Q, R};
}

template<class K>
FPS<K> gcd(FPS<K> const &A, FPS<K> const &B) {
    if(B == FPS<K>(0)) return A;
    FPS<K> Q, R;
    polynomial::euclid_division(A, B, Q, R);
    return polynomial::gcd(B, R);
}

template<class K>
FPS<K> extended_gcd(FPS<K> const &A, FPS<K> const &B, FPS<K> &U, FPS<K> &V) {
    FPS<K> D(A), V1 = FPS<K>(0), V3 = FPS<K>(B), T, Q, R;
    U = FPS<K>(1);
    while(!V3.zero()) {
        polynomial::euclid_division(D, V3, Q, R);
        T = U - V1 * Q;
        U = V1, D = V3;
        V1 = T, V3 = R;
    }
    V = (D - A*U) / B;
}

} // namespace polynomial


} // namespace zeno

#endif /* End of ZENO_POLY_HPP*/