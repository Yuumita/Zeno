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

/**
 * @brief The Polynomial class
 * 
 * @tparam R The type of coefficients
 */
template<typename R>
class Poly {
private:
    /**
     * @brief Coefficients of the polynomial
     */
    std::vector<R> a;
public:
    /**
     * @brief Reduce the 0 coefficients of the Polynomial (i.e remove all trailing zeroes)
     */
    void normalize(){
        while(!this->a.empty() && this->a.back() == R(0))
            this->a.pop_back();
    }

    /**
     * @brief Construct a new Polynomial object with degree -1 (No coefficients) 
     */
    Poly(): a({}) {}

    /**
     * @brief Construct a new Polynomial object with coefficients from the list coefs 
     * 
     * @param coefs The coefficient list
     * @param n The size of the list
     */
    Poly(const int n, const R *coefs) {
        this->a.resize(n);
        for(int i = 0; i < n; i++) a[i] = coefs[i];
        this->normalize();
    }

    /**
     * @brief Construct a new Polynomial object with coefficients from the array coefs 
     * 
     * @param coefs The coefficient array
     */
    Poly(const std::vector<R> &coefs): a(coefs) {
        this->normalize();
    }

    /**
     * @brief Construct a new (constant, i.e of 0 degree) Polynomial object which is just v
     * 
     * @param v The value of the constant polynomial
     */
    Poly(const R &v) {
        this->a.assign(1, v);
        this->normalize();
    }

    /**
     * @brief Construct a new Polynomial object of degree n-1 with initial values v 
     * 
     * @param n The coefficient count
     * @param v The value to initialize the coefficients
     */
    Poly(const int n, const R &v) {
        this->a.resize(n);
        for(int i = 0; i < n; i++) a[i] = v;
    }

    /**
     * @brief Return the degree of the polynomial (-1 if there aren't any coefficients)
     * 
     * @return int The degree of the polynomial
     */
    int deg() const {
        return this->a.size() - 1;
    }

    int is_zero() const {
        return this->a.empty();
    }

    /**
     * @brief Get the coefficient of x^i
     * 
     * @param i The index of the coefficient
     * @return R The coefficient
     */
    R get(const int i) const {
        return (0 <= i && i < this->a.size() ? this->a[i] : R(0));
    }

    /**
     * @brief Set the coefficient of x^i
     * 
     * @param i The index of the coefficient
     * @param c The value of the new coefficient
     */
    void set(const int i, const R c) {
        assert(0 <= i);
        if(this->a.size() <= i) this->a.resize(i+1);
        this->a[i] = c;
    }

    /**
     * @brief Addition of two polynomials (*this and P)
     * 
     * @param P The other polynomial
     * @return Poly& The resulting polynomial
     */
    Poly plus(const Poly& P) const {
        return *this + P;
    }


    /**
     * @brief Multiplication of two polynomials (*this and P)
     * @return Poly& The resulting polynomial
     */
    Poly times(const Poly& P) const {
        return *this * P;
    }

    /**
     * @brief Multiplication of this polynomial with the scalar v
     * @return Poly The resulting polynomial
     */
    Poly times(const R& v) const {
        return *this * v;
    }

    /**
     * @brief Subtraction of two polynomials (*this - P)
     * @return Poly The resulting polynomial
     */
    Poly minus(const Poly& P) const {
        return *this - P;
    }

    /**
     * @return Poly A copy of this polynomial mod x^k (getting the first k coefficients)
     */
    Poly mod_xk(int k) const {
        return Poly(std::vector<R>(a.begin(), a.begin() + std::min(k, (int)a.size())));
    }

    /**
     * @return Poly A copy of this polynomial times x^k (coefficients shifted k steps to the right)
     */
    Poly times_xk(int k) const {
        std::vector<int> coefs = this->a;
        coefs.insert(coefs.begin(), k, 0);
        return Poly(coefs);
    }

    /**
     * @return Poly A copy of this polynomial with coefficients shifted k steps to the left
     */
    Poly floordiv_xk(int k) const {
        return Poly(std::vector<R>(a.begin() + k, a.end()));
    }


    /**
     * @return Poly 
     */
    Poly sqrt_x() const {
        Poly P(*this);
        for(int i = 0; i <= P.deg(); i++) {
            P.coef(i) = P[2*i];
        }
        P.normalize();
        return P;
    }

    /**
     * @return Poly 
     */
    Poly sqr_x() const {
        Poly P(*this);
        P.a.resize(2*P.deg() + 1);
        for(int i = 0; i <= 2*P.deg(); i += 2) {
            P.set(i, P[i/2]);
        }
        P.normalize();
        return P;
    }

    /**
     * @return Poly The multiplicative inverse of this polynomial mod x^k
     */
    Poly inverse(int k) const {
        Poly A = this->mod_xk(k);
        if(k == 1) return Poly(1 / A[0]);
        Poly Am = A.scale_x(R(-1));
        Poly B = Am * A;
        B = B.sqrt_x().inverse(k / 2).sqr_x();
        return Am * B;
    }

    /**
     * @return Poly  A copy of this polynomial with coefficients reversed
     */
    Poly reverse() const {
        return Poly(std::vector<R>(a.rbegin(), a.rend()));
    }

    /**
     * @brief Euclidean division of this polynomial with the polynomial P
     * @return pair<Poly, Poly> the quotiend (first) and remainder (second) polynomials
     */
    std::pair<Poly, Poly> div(const Poly& P) const {
        assert(!P.is_zero());
        if(P.deg() > this->deg()) {
            return {Poly(0), Poly(*this)};
        }
        Poly A_R = this->reverse(), B_R = P.reverse();
        int d = this->deg() - P.deg();
        Poly D = ((A_R).mod_xk(d + 1) * P.inverse(d + 1)).mod_xk(d + 1).reverse();
        return {D, *this - D * P};
    }

    /**
     * @brief component-wise multiplication with v^k
     * @return Poly A copy of the polynomial with the x value scaled by v
     */ 
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
        this->a = fft::convolution(this->a, rhs.a);
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

    Poly operator-() const {
        auto t = *this;
        for(auto &e: t.a) e = -e;
        return t;
    }

    Poly operator+(const Poly &rhs) const { return Poly(*this) += rhs; }
    Poly operator-(const Poly &rhs) const { return Poly(*this) -= rhs; }
    Poly operator*(const Poly &rhs) const { return Poly(*this) *= rhs; }
    Poly operator*(const R &rhs) const { return Poly(*this) *= rhs; }
    Poly operator/(const Poly &rhs) const { return Poly(*this) /= rhs; }
    Poly operator/(const R &rhs) const { return Poly(*this) /= rhs; }
    Poly operator%(const Poly &rhs) const { return Poly(*this) %= rhs; }

    R operator[](const int index) const { return this->get(index); }

    // mutable reference of the coefficient
    R& coef(size_t index) { return this->a[index]; }

    bool operator == (const Poly &rhs) const { return a == rhs.a; }
    bool operator != (const Poly &rhs) const { return a != rhs.a; }


    // Computes the derivative
    Poly deriv(int k = 1) {
        if(deg() + 1 < k) return Poly(R(0));
        std::vector<R> res(deg() - k +  1);
        if(k == 1) {
            for(int i = 1; i <= deg(); i++) {
                res[i - 1] = a[i] * i;
            }
        } else {
            for(int i = k; i <= deg(); i++) {
                res[i - k] = a[i] * internal::fact<R>(i) / internal::fact<R>(i - k);
            }
        }
        return Poly(res);
    }

    // Compute integral with constant = 0
    Poly integr() {
        std::vector<R> res(deg() + 2);
        for(int i = 1; i <= deg() + 1; i++) {
            res[i] = a[i-1] / i;
        }
        return Poly(res);
    }


    /**
     * @brief Evaluate the polynomial in point x0
     * 
     * @param x0 The point to evaluate
     * @return R The value of the polynomial at point x0 
     */
    R eval(R x0) const {
        R y(0);
        for(int i = this->deg(); i >= 0; i--) {
            y = (x0 * y) + a[i];
        }
        return y;
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

} // namespace zeno

#endif /* End of ZENO_POLY_HPP*/