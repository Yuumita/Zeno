#pragma once

#include <cstdint>
#include <stdlib.h>

#include "modular.hpp" // is_modular class

namespace zeno {
    

/// TODO: check correctness.
/// @brief Implementation of Montgomery space. A number X converted to montgometry(X): (x*r) mod n.
template <typename Z, typename uZ, typename longZ, typename ulongZ, uZ MOD>
class Montgomery {
    using M = Montgomery;

private:
    Z logsize;
    uZ modinv; // modinv = mod^{-1} mod r
    uZ val;

    uZ reduce(ulongZ const &x) const {
        uZ q = uZ(x) * modinv;               // q = xn' mod r
        uZ m = (ulongZ(q) * MOD) >> logsize; // m = qn
        return (x >> logsize) + MOD - m;
    }

public:

    static constexpr uZ mod() { return MOD; }

    Montgomery() : val(0) {
        modinv = MOD;
        while(modinv * MOD != 1) modinv *= uZ(2) - MOD * modinv;
    }

    Montgomery(uZ val_) : Montgomery() {
        val = val_;
    }

    template <typename T, std::enable_if_t<std::is_integral<T>::value>* = nullptr>
    Montgomery(T val_) {
        longZ x = static_cast<longZ>(val_ % static_cast<longZ>(MOD));
        if (x < 0) x += MOD;
        val = static_cast<uZ>(x);
    }

    template <typename T, std::enable_if_t<!std::is_integral<T>::value>* = nullptr>
    Montgomery(T val_): Montgomery(static_cast<longZ>(val_)) { }

    M &operator+=(M const &b) {
        val += b.val;
        while(val >= MOD) val -= MOD;
        return *this;
    }

    M &operator-=(M const &b) {
        val += MOD - b.val;
        while(val >= MOD) val -= MOD;
        return *this;
    }

    M &operator*=(M const &b) {
        val = reduce(ulongZ(val) * b.val);
        return *this;
    }

    M &operator/=(M const &b) {
        *this *= b.inv();
        return *this;
    }

    M operator+(const M &b) const { return M(*this) += b; }
    M operator-(const M &b) const { return M(*this) -= b; }
    M operator*(const M &b) const { return M(*this) *= b; }
    M operator/(const M &b) const { return M(*this) /= b; }

    M operator+() const { return M(*this); }
    M operator-() const { return M(0) - M(*this); }


    friend bool operator==(const M &A, const M &B) { return A.retrieve() == B.retrieve(); }
    friend bool operator!=(const M &A, const M &B) { return A.retrieve() != B.retrieve(); }

    friend std::ostream &operator<<(std::ostream &os, const M &x) { return os << x.retrieve(); }
    friend std::istream &operator>>(std::istream &is, M &x) {
        longZ t; is >> t;
        x = M(t);
        return (is);
    }

    M pow(Z n) const { /* change it so it uses zeno::pow ??? */
        M ret(1), mul(*this);
        if(n < 0) mul = mul.inv(), n = -n;
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul, n >>= 1;
        }
        return ret;
    }

    /// TODO: prime check for fermat's theorem ???
    M inv() const {
        Z x = retrieve(), y = MOD, u = 1, v = 0;
        while (y > 0) {
            Z t = x / y;
            std::swap(x -= t * y, y);
            std::swap(u -= t * v, v);
        }
        return M(u);
    }

    operator int() const { return retrieve(); }
    uZ retrieve() const {
        uZ ret = reduce(val);
        while(ret >= MOD) ret -= MOD;
        return ret;
    }

};


template <typename Z, typename uZ, typename longZ, typename ulongZ, uZ mod>
struct is_modular<zeno::Montgomery<Z, uZ, longZ, ulongZ, mod>> : std::true_type {};

template<uint32_t mod>
using Montgomery32 = Montgomery<int32_t, uint32_t, int64_t, uint64_t, mod>;



} // namespace zeno
