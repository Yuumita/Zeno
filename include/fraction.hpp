#ifndef ZENO_STRUCTURES_HPP
#define ZENO_STRUCTURES_HPP

#include <iostream>
#include <utility> // ???
#include <vector>
#include "euclid.hpp"
// #include <math.h> // acos
#include <algorithm> // std::swap

namespace zeno {

template<class Z = int64_t>
class fraction {
private:
    Z num, den;  // numerator, denominator
    void normalize() {
        Z g = zeno::gcd<Z>(num, den);
        num /= g, den /= g;
    }
public:
    fraction() : num(Z(0)), den(Z(1)) {}
    fraction(Z num_): num(num_), den(1) {}

    template<typename T>
    fraction(T n): num(Z(n)), den(1) {}

    fraction(Z num_, Z den_): num(num_), den(den_) {
        normalize();
    }
    

    fraction& operator+=(fraction const &rhs) { 
        num = (num * rhs.den) + (rhs.num * den) ;
        den = den * rhs.den;
        normalize();
        return *this;
    }

    fraction& operator++() { 
        num += den;
        normalize();
        return *this;
    }

    fraction operator++(int) { 
        fraction tmp(*this);
        operator++();
        return tmp;
    }

    fraction operator-() const { 
        return fraction(0) -= *this;
    }

    fraction operator-(fraction f) { 
        return fraction(*this) -= f;
    }

    fraction& operator-=(fraction const &rhs) { 
        num = (num * rhs.den) - (rhs.num * den) ;
        den = den * rhs.den;
        normalize();
        return *this;
    }

    fraction& operator--() { 
        num -= den;
        normalize();
        return *this;
    }

    fraction operator--(int) { 
        fraction tmp(*this);
        operator--();
        return tmp;
    }

    fraction& operator*=(fraction const &rhs) { 
        num = num * rhs.num;
        den = den * rhs.den;
        normalize();
        return *this;
    }

    fraction& operator/=(fraction const &rhs) { 
        num = num * rhs.den;
        den = den * rhs.num;
        normalize();
        return *this;
    }


    /* casting */
    template<class T> fraction  operator+(T x) const { return fraction(*this) += x; }
    template<class T> fraction  operator-(T x) const { return fraction(*this) -= x; }
    template<class T> fraction  operator*(T x) const { return fraction(*this) *= x; }
    template<class T> fraction  operator/(T x) const { return fraction(*this) /= x; }
    template<class T> fraction& operator+=(T x)      { return operator+=(fraction(x)); }
    template<class T> fraction& operator-=(T x)      { return operator-=(fraction(x)); }
    template<class T> fraction& operator*=(T x)      { return operator*=(fraction(x)); }
    template<class T> fraction& operator/=(T x)      { return operator/=(fraction(x)); }

    friend std::ostream &operator<<(std::ostream &os, const fraction &f) { return os << f.num << "/" << f.den; }
};
using rational = fraction<int64_t>;






} // namespace zeno

#endif /* ZENO_STRUCTURES_HPP */