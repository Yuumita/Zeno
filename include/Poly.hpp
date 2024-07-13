/**
 * @file Poly.hpp
 * @author Yuumita
 * @brief Polynomial Interface
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
class Poly {
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
    Poly(): a({}) {}

    /// @brief Construct a new Polynomial object with coefficients from the list coefs 
    /// @param coefs The coefficient list
    /// @param n The size of the list
    Poly(const int n, const R *coefs) {
        this->a.resize(n);
        for(int i = 0; i < n; i++) a[i] = coefs[i];
        this->normalize();
    }

    /// @brief Construct a new Polynomial object with coefficients from the array coefs 
    ///  
    /// @param coefs The coefficient array
    Poly(const std::vector<R> &coefs): a(coefs) {
        this->normalize();
    }

    /// @brief Construct a new (constant, i.e of 0 degree) Polynomial object which is just v
    /// @param v The value of the constant polynomial
    Poly(const R &v) {
        this->a.assign(1, v);
        this->normalize();
    }

    /// @brief Construct a new Polynomial object of degree n-1 with initial values v 
    /// @param n The coefficient count
    /// @param v The value to initialize the coefficients
    Poly(const int n, const R &v) {
        this->a.resize(n);
        for(int i = 0; i < n; i++) a[i] = v;
    }

    /// @brief Return the degree of the polynomial (-1 if there aren't any coefficients)
    /// @return int The degree of the polynomial
    int deg() const {
        return this->a.size() - 1;
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
    /// @return Poly& The resulting polynomial
    Poly plus(const Poly& P) const {
        return *this + P;
    }


    /// @brief Multiplication of two polynomials (*this and P)
    /// @return Poly& The resulting polynomial
    Poly times(const Poly& P) const {
        return *this * P;
    }

    /// @brief Multiplication of this polynomial with the scalar v
    /// @return Poly The resulting polynomial
    Poly times(const R& v) const {
        return *this * v;
    }

    /// @brief Subtraction of two polynomials (*this - P)
    /// @return Poly The resulting polynomial
    Poly minus(const Poly& P) const {
        return *this - P;
    }

    /// @return Poly A copy of this polynomial mod x^k (getting the first k coefficients)
    Poly mod_xk(int k) const {
        return Poly(std::vector<R>(a.begin(), a.begin() + std::min(k, (int)a.size())));
    }

    /// @return Poly A copy of this polynomial times x^k (coefficients shifted k steps to the right)
    Poly times_xk(int k) const {
        std::vector<int> coefs = this->a;
        coefs.insert(coefs.begin(), k, 0);
        return Poly(coefs);
    }

    /// @return Poly A copy of this polynomial with coefficients shifted k steps to the left
    Poly floordiv_xk(int k) const {
        return Poly(std::vector<R>(a.begin() + k, a.end()));
    }

    Poly shift(int k) const {
        if(k > 0)
            return times_xk(k);
        if(k < 0)
            return floordiv_xk(-k);
        return Poly(*this);
    }


    /// @return Poly P(x^{1/2})
    Poly sqrt_x() const {
        std::vector<R> ret(this->a.size(), R(0));
        for(int i = 0; i < this->a.size(); i++)
            ret[i] = get(2*i);
        return Poly(ret);
    }

    /// @return Poly P(x^2)
    Poly sqr_x() const {
        std::vector<R> ret(2 * this->a.size(), R(0));
        for(int i = 0; i < this->a.size(); i++)
            ret[2*i] = this->a[i];
        return Poly(ret);
    }

    /// @return Poly The multiplicative inverse of this polynomial mod x^m
    Poly inv(int m = deg() + 1) const {
        assert(get(0) != R(0));
        if(m == 0) return Poly(0);
        Poly Q(R(1) / get(0));
        for(int i = 1; i < m; i *= 2)
            Q = (Q * R(2) - Q * Q * this->mod_xk(2*i)).mod_xk(2*i);
        return Q.mod_xk(m);
        /*
            if(m == 1) return Poly(R(1) / get(0)); 
            Poly A = this->mod_xk(m);
            Poly Am = A.scale_x(R(-1));
            Poly B = Am * A;
            B = B.sqrt_x().inv((m+1) / 2).sqr_x();
            std::cerr << "inv (m = " << m << ") --> " << (Am * B).mod_xk(m) << std::endl;
            std::cerr << A.mod_xk(m) << " * " << (Am*B).mod_xk(m) << " = " << (A*Am*B).mod_xk(m) << std::endl;
            std::cerr << "inv (m = " << m << ") --> " << (Am * B).mod_xk(m) << std::endl;
            return (Am * B).mod_xk(m);
        */
    }

    /// @return Poly  A copy of this polynomial with coefficients reversed
    Poly reverse() const {
        return Poly(std::vector<R>(a.rbegin(), a.rend()));
    }

    /// @brief Euclidean division of this polynomial with the polynomial P
    /// @return pair<Poly, Poly> the quotiend (first) and remainder (second) polynomials
    std::pair<Poly, Poly> div(const Poly& P) const {
        assert(!P.is_zero());
        if(P.deg() > this->deg())
            return {Poly(0), Poly(*this)};
        Poly A_R = this->reverse(), B_R = P.reverse();
        int d = this->deg() - P.deg();
        Poly D = ((A_R).mod_xk(d + 1) * P.inv(d + 1)).mod_xk(d + 1).reverse();
        return {D, *this - D * P};
    }

    /// @brief component-wise multiplication with v^k
    /// @return Poly A copy of the polynomial with the x value scaled by v
    Poly scale_x(const R& v) const {
        R c = 1;
        Poly P(*this);
        for(int i = 0; i <= P.deg(); i++, c *= v) { 
            P.coef(i) *= c;
        }
        return  P;
    }

    Poly& operator+=(const Poly &rhs) {
        int d = (this->deg() > rhs.deg() ? this->deg() : rhs.deg());
        this->a.resize(d);
        for(int i = 0; i < (int)rhs.a.size(); i++){
            this->a[i] += rhs.a[i];
        }
        this->normalize();
        return *this;
    }

    Poly& operator-=(const Poly &rhs) {
        int d = (this->deg() > rhs.deg() ? this->deg() : rhs.deg());
        this->a.resize(d + 1);
        for(int i = 0; i < (int)rhs.a.size(); i++){
            this->a[i] -= rhs.a[i];
        }
        this->normalize();
        return *this;
    }

    Poly& operator*=(const Poly &rhs) {
        this->a = fft::convolution<R>(this->a, rhs.a);
        this->normalize();
        return *this;
    }

    Poly& operator*=(const R &rhs) {
        for(auto &e: this->a) e *= rhs;
        this->normalize();
        return *this;
    }

    Poly& operator/=(const Poly &rhs) {
        return *this = div(rhs).first;
    }

    Poly& operator/=(const R &rhs) {
        for(auto &e: this->a) e /= rhs;
        return *this;
    }

    Poly& operator%=(const Poly &rhs) {
        return *this = div(rhs).second;
    }

    // TODO: >>=, <<= shifts
    // Poly& operator<<=(const int s) {
    //     return *this = (*this).shift(s);
    // }
    // Poly& operator>>=(const int s) {
    //      return operator<<=(-s);
    // }

    Poly operator-() const {
        auto t = *this;
        for(auto &e: t.a) e = -e;
        return t;
    }

    Poly operator+(const Poly &rhs) const { return Poly(*this) += rhs; }
    Poly operator-(const Poly &rhs) const { return Poly(*this) -= rhs; }
    Poly operator*(const Poly &rhs) const { return Poly(*this) *= rhs; }
    Poly operator*(const R &rhs)    const { return Poly(*this) *= rhs; }
    Poly operator/(const Poly &rhs) const { return Poly(*this) /= rhs; }
    Poly operator/(const R &rhs)    const { return Poly(*this) /= rhs; }
    Poly operator%(const Poly &rhs) const { return Poly(*this) %= rhs; }

    // mutable reference of the coefficient
    R& coef(size_t index) { return this->a[index]; }

    bool operator == (const Poly &rhs) const { return a == rhs.a; }
    bool operator != (const Poly &rhs) const { return a != rhs.a; }


    // Computes the derivative
    Poly deriv(int k = 1) const {
        if(deg() + 1 < k) return Poly(R(0));
        std::vector<R> res(deg() - k +  1, R(0));
        if(k == 1) {
            for(int i = 1; i <= deg(); i++)
                res[i - 1] = get(i) * R(i);
        } else {
            for(int i = k; i <= deg(); i++)
                res[i - k] = get(i) * internal::fact<R>(i) / internal::fact<R>(i - k);
        }
        return Poly(res);
    }

    /// @return The integral of the polynomial (with constant 0)
    Poly integr() const {
        std::vector<R> res(deg() + 2);
        for(int i = 1; i <= deg() + 1; i++) {
            res[i] = get(i-1) / i;
        }
        return Poly(res);
    }

    /// @return The (natural) logarithm of the polynomial (mod x^m)
    Poly log(int m) const {
        return (Poly(*this).deriv().mod_xk(m) * Poly(*this).inv(m)).mod_xk(m); 
    }

    Poly exp(int m) const {
        assert(get(0) == 0);
        std::cerr << "T0" << std::endl;
        Poly Q(1), P(*this);
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

    Poly sqrt(int k = 2) const {

    }

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
    static Poly<R> build_poly_tree(std::vector<Poly<R>> &tree, 
        int v, typename std::vector<R>::iterator l, typename std::vector<R>::iterator r) {
            if(r - l == 1) {
                return tree[v] = Poly({-*l, 1});
            } else {
                auto m = l + (r - l) / 2;
                return tree[v] = build_poly_tree(tree, 2*v, l, m) 
                    * build_poly_tree(tree, 2*v+1, m, r);
            }
        }

    std::vector<R> eval(std::vector<Poly<R>> &tree, 
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
        std::vector<Poly<R>> tree(4*n);
        build_poly_tree(tree, 1, begin(x), end(x));
        return eval(tree, 1, begin(x), end(x));
    }


    static Poly<R> _interpolate_tree(std::vector<Poly<R>> const &tree, 
        int v, typename std::vector<R>::iterator l, typename std::vector<R>::iterator r) {
            if(r - l == 1) {
                return tree[v] = Poly<R>(*l);
            } else {
                auto m = l + (r - l) / 2;
                Poly<R> A = _interpolate_tree(tree, 2*v, l, m);
                Poly<R> B = _interpolate_tree(tree, 2*v+1, m, r);
                return tree[v] = A * tree[2*v+1] + tree[2*v] * B;
            }
    }


    static Poly<R> interpolate(std::vector<R> const &x, std::vector<R> const &y) {
        assert(x.size() == y.size());
        std::vector<Poly<R>> tree(4*x.size());
        std::vector<R> u = Poly<R>::build_poly_tree(tree, 1, begin(x), end(x)).deriv().eval(x);
        for(int i = 0; i < y.size(); i++)
            u[i] = y[i] / u[i];
        return _interpolate_tree(tree, 1, begin(u), end(u));
    }


    /// @return The Borel transformation of the polynomial (a_k x^k -> a_k / k! x^k)
    Poly<R> borel() const {
        Poly<R> ret = *this;
        for(int k = 0; k <= deg(); k++)
            ret->a[k] /= internal::fact(k);
        return ret;
    }

    /// @return The inverse Borel (Laplace) transformation of the polynomial (a_k x^k -> k!a_k x^k)
    Poly<R> invborel() const {
        Poly<R> ret = *this;
        for(int k = 0; k <= deg(); k++)
            ret->a[k] *= internal::fact(k);
        return ret;
    }

};

template<typename R>
Poly<R> operator*(const R &lhs, const Poly<R> &rhs) {
    return Poly(rhs) *= lhs;
}

template<typename R>
Poly<R> operator/(const R &lhs, const Poly<R> &rhs) {
    return Poly(lhs) /= rhs;
}


template<typename R>
std::ostream& operator<<(std::ostream& os, const Poly<R> &P) {
    os << "[ ";
    for(int i = 0; i <= P.deg(); i++)
        os << P.get(i) << " ";
    os << "] ";
    return os;
}

namespace polynomial
{

template<class K>
void euclid_division(Poly<K> const &A, Poly<K> const &B, Poly<K> &Q, Poly<K> &R) {
    R = A, Q = Poly<K>(0);
    while(R.deg() >= B.deg()) {
        K s = R.leadcoef() / B.leadcoef();
        int sd = R.deg() - B.deg();

        Q.set(sd, Q.get(sd) + s);
        R = R - s * B.shift(sd);
    }
}

template<class K>
std::pair<Poly<K>, Poly<K>> euclid_division(Poly<K> const &A, Poly<K> const &B) {
    Poly<K> Q, R;
    polynomial::euclid_division(A, B, Q, R);
    return {Q, R};
}

template<class K>
Poly<K> gcd(Poly<K> const &A, Poly<K> const &B) {
    if(B == Poly<K>(0)) return A;
    Poly<K> Q, R;
    polynomial::euclid_division(A, B, Q, R);
    return polynomial::gcd(B, R);
}

template<class K>
Poly<K> extended_gcd(Poly<K> const &A, Poly<K> const &B, Poly<K> &U, Poly<K> &V) {
    Poly<K> D(A), V1 = Poly<K>(0), V3 = Poly<K>(B), T, Q, R;
    U = Poly<K>(1);
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