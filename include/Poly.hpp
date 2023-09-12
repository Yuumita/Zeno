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

#include "fft.hpp"
#include "internal_math.hpp"


namespace zeno {

/**
 * @brief The Polynomial class
 * 
 * @tparam T The type of coefficients
 */
template<typename T>
class Poly {
private:
    /**
     * @brief Coefficients of the polynomial
     */
    std::vector<T> a;
public:
    /**
     * @brief Reduce the 0 coefficients of the Polynomial (i.e remove all trailing zeroes)
     */
    void normalize(){
        while(!this->a.empty() && this->a.back() == T(0)) {
            this->a.pop_back();
        }
        if(this->a.size() == 0) this->a = {0};
    }

    /**
     * @brief Construct a new Polynomial object with degree -1 (No coefficients) 
     */
    Poly(): a({0}){}

    /**
     * @brief Construct a new Polynomial object with degree n (but all coefficients 0) 
     */
    Poly(const int n) {
        this->a.resize(n);
    }

    /**
     * @brief Construct a new Polynomial object with coefficients from the list coefs 
     * 
     * @param coefs The coefficient list
     * @param n The size of the list
     */
    Poly(const int n, const T *coefs) {
        this->a.resize(n);
        for(int i = 0; i < n; i++) a[i] = coefs[i];
        this->normalize();
    }

    /**
     * @brief Construct a new Polynomial object with coefficients from the array coefs 
     * 
     * @param coefs The coefficient array
     */
    Poly(const std::vector<T> &coefs): a(coefs) {
        this->normalize();
    }

    /**
     * @brief Construct a new (constant, i.e of 0 degree) Polynomial object which is just v
     * 
     * @param v The value of the constant polynomial
     */
    Poly(const T &v) {
        this->a.assign(std::vector<T>({v}));
        this->normalize();
    }

    /**
     * @brief Construct a new Polynomial object of degree n-1 with initial values v 
     * 
     * @param n The coefficient count
     * @param v The value to initialize the coefficients
     */
    Poly(const int n, const T &v) {
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

    /**
     * @brief Get the coefficient of x^i
     * 
     * @param i The index of the coefficient
     * @return T The coefficient
     */
    T get(const int i) const {
        return (0 <= i && i < this->a.size() ? this->a[i] : T(0));
    }

    /**
     * @brief Set the coefficient of x^i
     * 
     * @param i The index of the coefficient
     * @param c The value of the new coefficient
     */
    void set(const int i, const T c) {
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
     * 
     * @param P The other polynomial
     * @return Poly& The resulting polynomial
     */
    Poly times(const Poly& P) const {
        return *this * P;
    }

    /**
     * @brief Multiplication of this polynomial with the scalar v
     * 
     * @param v The scalar 
     * @return Poly& The resulting polynomial
     */
    Poly times(const T& v) const {
        return *this * v;
    }

    /**
     * @brief Subtraction of two polynomials (*this - P)
     * 
     * @param P The other polynomial
     * @return Poly& The resulting polynomial
     */
    Poly minus(const Poly& P) const {
        return *this - P;
    }

    /**
     * @brief Division of two polynomials (*this / P)
     * 
     * @param P The other polynomial
     * @return Poly& The resulting polynomial
     */
    Poly div(const Poly& P) const {
        return *this / P;
    }

    /**
     * @brief Division of this polynomial with the scalar v
     * 
     * @param v The scalar 
     * @return Poly& The resulting polynomial
     */
    Poly div(const T& v) const {
        return *this / v;
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

    /* TODO: Multiplication */
    Poly& operator*=(const Poly &rhs) {
        this->a = fft::convolution(this->a, rhs.a);
        this->normalize();
        return *this;
    }

    Poly& operator*=(const T &rhs) {
        for(auto &e: this->a) e *= rhs;
        this->normalize();
        return *this;
    }

    /* TODO: Division */
    Poly& operator/=(const Poly &rhs) {
        return *this;
    }

    Poly& operator/=(const T &rhs) {
        for(auto &e: this->a) e /= rhs;
        return *this;
    }

    Poly operator-() const {
        auto t = *this;
        for(auto &e: t.a) e = -e;
        return t;
    }

    Poly operator+(const Poly &rhs) const { return Poly(*this) += rhs; }
    Poly operator-(const Poly &rhs) const { return Poly(*this) -= rhs; }
    Poly operator*(const Poly &rhs) const { return Poly(*this) *= rhs; }
    Poly operator*(const T &rhs) const { return Poly(*this) *= rhs; }
    Poly operator/(const Poly &rhs) const { return Poly(*this) *= rhs; }
    Poly operator/(const T &rhs) const { return Poly(*this) *= rhs; }

    T operator[](const int index) const { return this->get(index); }

    // mutable reference of the coefficient
    T& coef(size_t index) { return this->a[index]; }

    bool operator == (const Poly &rhs) const { return a == rhs.a; }
    bool operator != (const Poly &rhs) const { return a != rhs.a; }


    // Computes the derivative
    Poly deriv(int k = 1) {
        if(deg() + 1 < k) return Poly(T(0));
        vector<T> res(deg() - k +  1);
        if(k == 1) {
            for(int i = 1; i <= deg(); i++) {
                res[i - 1] = a[i] * i;
            }
        } else {
            for(int i = k; i <= deg(); i++) {
                res[i - k] = a[i] * internal::fact<T>(i) * internal::invfact(i - k);
            }
        }
        return Poly(res);
    }

    // Compute integral with constant = 0
    Poly integr() {
        vector<T> res(deg() + 2);
        for(int i = 1; i <= deg() + 1; i++) {
            res[i] = a[i-1] / i;
        }
        return Poly(res);
    }


    /**
     * @brief Evaluate the polynomial in point x0
     * 
     * @param x0 The point to evaluate
     * @return T The value of the polynomial at point x0 
     */
    T eval(T x0) const {
        T y(0);
        for(int i = this->deg(); i >= 0; i--) {
            y = (x0 * y) + a[i];
        }
        return y;
    }
};

template<typename T>
Poly<T> operator*(const T &lhs, const Poly<T> &rhs) {
    return Poly(rhs) *= lhs;
}

template<typename T>
Poly<T> operator/(const T &lhs, const Poly<T> &rhs) {
    return Poly(rhs) /= lhs;
}


template<typename T>
std::ostream& operator<<(std::ostream& os, const Poly<T> &P) {
    os << "[ ";
    for(int i = 0; i <= P.deg(); i++)
        os << P.get(i) << " ";
    os << "] ";
    return os;
}

} // namespace zeno

#endif /* End of ZENO_POLY_HPP*/