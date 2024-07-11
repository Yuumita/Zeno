#ifndef ZENO_STRUCTURES_HPP
#define ZENO_STRUCTURES_HPP

#include <iostream>
#include <utility> // ???
#include <vector>
// #include <math.h> // acos
#include <algorithm> // std::swap

namespace zeno {


template<typename R>
class Fraction {

private:
    R num, den;  // numerator, denominator
    void normalize() {

    }
public:
    Fraction<R>() : num(R(0)), den(R(1)) {}
    Fraction<R>(R num_, R den_): num(num_), den(den_) {
        normalize();
    }
    Fraction<R>(R num_): num(num_), den(1) {}
    Fraction<R>(long long n): num(R(n)), den(1) {}


    Fraction& operator+=(Fraction const &rhs) { 
        num = (num * rhs.den) + (rhs.num * den) ;
        den = den * rhs.den;
        normalize();
    }
    Fraction& operator++() { 
        return Fraction(*this) += R(1);
    }

};




template <int MOD = 998244353, bool is_prime = true>
struct modular {
    int v;
    modular(): v(0) {}
    modular(int vv) : v(vv < 0 ? vv % MOD + MOD : vv % MOD) {}
    modular(long long vv) : v(static_cast(vv < 0 ? vv % MOD + MOD : vv % MOD)) {}
    modular(int vv, std::nullptr_t) : v(vv) {}
    operator int() const { return v; }

    modular& operator+=(modular x) { v += x.v; if(v >= MOD) v -= MOD; return *this; }
    modular& operator++() { if(v == MOD - 1) v = 0; else v++; return *this; }
    modular operator++(int) { modular ans(*this); operator++(); return ans; }
    modular operator-() const { return modular(0) -= *this; }
    modular operator-(modular x) const { return modular(*this) -= x; }
    modular& operator-=(modular x) { if(v < x.v) v += MOD; v -= x.v; return *this; }
    modular& operator--() { if(v == 0) v = MOD - 1; else v--; return *this; }
    modular operator--(int) { modular ans(*this); operator--(); return ans; }
    modular& operator*=(modular x) { v = 1ll * v * x.v % MOD; return *this; }
    modular& operator/=(modular x) { return operator*=(x.inverse()); }

    /* casting */
    template<class T> modular operator+(T x) const { return modular(*this) += x; }
    template<class T> modular& operator+=(T x) { return operator+=(modular(x)); }
    template<class T> modular operator-(T x) const { return modular(*this) -= x; }
    template<class T> modular& operator-=(T x) { return operator-=(modular(x)); }
    template<class T> modular operator*(T x) const { return modular(*this) *= x; }
    template<class T> modular& operator*=(T x) { return operator*=(modular(x)); }
    template<class T> modular operator/(T x) const { return modular(*this) /= x; }
    template<class T> modular& operator/=(T x) { return operator/=(modular(x)); }

    modular pow(int n) const { /* change it so it uses zeno::pow */
        modular ret(1), mul(v);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul, n /= 2;
        }
        return ret;
    }

    modular inv() const {

        // Fermat's Little Theorem
        if(is_prime) return pow(MOD - 2);

        int a = v, b = MOD, u = 1, v = 0, t;
        while (b > 0) {
            t = a / b;
            std::swap(a -= t * b, b);
            std::swap(u -= t * v, v);
        }
        return modular(u);
    }

    modular inverse() const { return inv(); }

    int get_mod() { return MOD; }

    friend std::ostream &operator<<(std::ostream &os, const modular &p) { return os << p.v; }
    friend std::istream &operator>>(std::istream &is, modular &a) {
        int t; is >> t;
        a = modular(t);
        return (is);
    }
};

template<int M>
using modint = modular<M, false>; 

using modint998244353 = modular<998244353, true>; 

} // namespace zeno

#endif /* ZENO_STRUCTURES_HPP */