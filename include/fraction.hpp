#pragma once

#include <iostream>
#include <utility> // ???
#include <vector>
#include "euclid.hpp"
// #include <math.h> // acos
#include <algorithm> // std::swap

namespace zeno {

template<class Z = int64_t, class longZ = __int128_t>
class FractionLong {
    using fraction = FractionLong;
private:
    Z num, den;  // numerator, denominator
    void normalize() {
        Z g;
        if(std::is_integral_v<Z>) {
            g = zeno::_gcd<Z>(num, den);
        } else {
            g = zeno::gcd<Z>(num, den);
        }
        if(g) num /= g, den /= g;
        if(den < 0) num = -num, den = -den;
    }
public:
    fraction() : num(Z(0)), den(Z(1)) {}
    fraction(Z num_): num(num_), den(1) {}

    template<typename T>
    fraction(T n): num(Z(n)), den(1) {}

    fraction(Z num_, Z den_): num(num_), den(den_) {
        assert(den != 0);
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

    fraction operator+() const { return fraction(*this); }
    fraction operator-() const { return fraction(0) -= *this; }

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
   friend fraction  operator+(const  fraction &x) const { return fraction(*this) += x; }
   friend fraction  operator-(const  fraction &x) const { return fraction(*this) -= x; }
   friend fraction  operator*(const  fraction &x) const { return fraction(*this) *= x; }
   friend fraction  operator/(const  fraction &x) const { return fraction(*this) /= x; }
   template<typename T> friend fraction& operator+=(const T &x) { return operator+=(fraction(x)); }
   template<typename T> friend fraction& operator-=(const T &x) { return operator-=(fraction(x)); }
   template<typename T> friend fraction& operator*=(const T &x) { return operator*=(fraction(x)); }
   template<typename T> friend fraction& operator/=(const T &x) { return operator/=(fraction(x)); }



    friend bool operator==(const R& l, const R& r) { return l.x == r.x && l.y == r.y; };
    friend bool operator!=(const R& l, const R& r) { return l.x != r.x || l.y != r.y; };

    friend bool operator<(const R& l, const R& r) {
        return longZ(l.num) * r.den < longZ(l.den) * r.num;
    };
    friend bool operator<=(const R& l, const R& r) { return l < r || l == r; }

    friend bool operator>(const R& l, const R& r) {
        return longZ(l.num) * r.den > longZ(l.den) * r.num;
    };
    friend bool operator>=(const R& l, const R& r) { return l > r || l == r; }

    friend ostream& operator<<(ostream& os, const R& r) {
        os << r.num;
        if (r.den != 0 && r.den != 1) os << "/" << r.den;
        return os;
    }

};

using rational   = FractionLong<int64_t, __int128_t>;
using rational32 = FractionLong<int32_t, int64_t>;
using rational64 = FractionLong<int64_t, __int128_t>;

template<typename Z>
using fraction = FractionLong<Z, Z>;


} // namespace zeno