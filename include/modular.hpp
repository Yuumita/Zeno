#ifndef ZENO_MODULAR_HPP
#define ZENO_MODULAR_HPP

#include <iostream>
#include <vector>
#include <math.h> // acos
#include <algorithm> // std::swap

namespace zeno {

template <int MOD = 998244353, bool is_prime = true>
struct modular {
    int v;
    modular(): v(0) {}
    modular(int vv) : v(vv < 0 ? vv % MOD + MOD : vv % MOD) {}
    modular(long long vv) : v(vv < 0 ? vv % MOD + MOD : vv % MOD) {}
    modular(int vv, nullptr_t) : v(vv) {}
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

    template<class T> modular operator+(T x) const { return modular(*this) += x; }
    template<class T> modular& operator+=(T x) { return operator+=(modular(x)); }
    template<class T> modular operator-(T x) const { return modular(*this) -= x; }
    template<class T> modular& operator-=(T x) { return operator-=(modular(x)); }
    template<class T> modular operator*(T x) const { return modular(*this) *= x; }
    template<class T> modular& operator*=(T x) { return operator*=(modular(x)); }
    template<class T> modular operator/(T x) const { return modular(*this) /= x; }
    template<class T> modular& operator/=(T x) { return operator/=(modular(x)); }

    modular pow(int n) const {
        modular ret(1), mul(v);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul, n /= 2;
        }
        return ret;
    }

    modular inverse() const {
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

    friend std::ostream &operator<<(std::ostream &os, const modular &p) { return os << p.v; }
    friend std::istream &operator>>(std::istream &is, modular &a) {
        int t; is >> t;
        a = modular(t);
        return (is);
    }
};

} // namespace zeno

#endif // ZENO_MODULAR_HPP