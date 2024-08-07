#pragma once

#include <cstdint>
#include <stdlib.h>

namespace zeno {
    

/// TODO: check correctness.
/// @brief Implementation of Montgomery space. A number X converted to montgometry(X): (x*r) mod n.
template <typename Z, typename uZ, typename longZ, typename ulongZ, uZ mod>
class Montgomery {
    using M = Montgomery;

private:
    Z logsize;
    uZ modinv; // modinv = mod^{-1} mod r
    uZ val;

    uZ reduce(ulongZ const &x) {
        uZ q = uZ(x) * modinv;               // q = xn' mod r
        uZ m = (ulongZ(q) * mod) >> logsize; // m = qn
        return (x >> logsize) + mod - m;
    }

public:

    Montgomery() : val(0) {
        modinv = mod;
        while(modinv * mod != 1) modinv *= uZ(2) - mod * modinv;
    }

    Montgomery(uZ val_) : val(val_) {
        Montgometry();
    }

    M &operator+=(M const &b) {
        val += b.val;
        while(val >= mod) val -= mod;
        return *this;
    }

    M &operator-=(M const &b) {
        val += mod - b.val;
        while(val >= mod) val -= mod;
        return *this;
    }

    M &operator*=(M const &b) {
        val = reduce(ulongZ(val) * b.val);
        return *this;
    }

    M &operator/=(M const &b) {
        *this *= b.inverse();
        return *this;
    }

    M operator+(const M &b) const { return M(*this) += b; }
    M operator-(const M &b) const { return M(*this) -= b; }
    M operator*(const M &b) const { return M(*this) *= b; }
    M operator/(const M &b) const { return M(*this) /= b; }

    M operator+() const { return M(*this); }
    M operator-() const { return M(0) - M(*this); }


    friend ostream &operator<<(ostream &os, const M &x) {
        return os << x.retrieve();
    }

    friend istream &operator>>(istream &is, M &x) {
        longZ t; is >> t;
        x = ArbitraryLazyMontgomeryModIntBase(t);
        return is;
    }

    /// TODO: prime check for fermat's theorem
    M inverse() const {
        Z x = retrieve(), y = mod, u = 1, v = 0;
        while (y > 0) {
            Z t = x / y;
            swap(x -= t * y, y);
            swap(u -= t * v, v);
        }
        return M(u);
    }

    uZ retrieve() {
        uZ ret = reduce(val);
        while(ret >= mod) ret -= mod;
        return ret;
    }
};

template<uint32_t mod>
using Montgomery32 = Montgomery<int32_t, uint32_t, int64_t, uint64_t, mod>;




} // namespace zeno
