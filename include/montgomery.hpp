#pragma once

#include <cstdint>
#include <stdlib.h>

#include "modular.hpp" // importing is_modular class

namespace zeno {
    


/// @brief Implementation of Montgomery space. Numbers X are mapped to (x*r) mod m where r = 2^logsize.
/// @tparam Z Integer ring.
/// @tparam uZ The unsigned integer domain of Z.
/// @tparam longZ Integer ring holding numbers at least as large as the squares of numbers in Z.
/// @tparam longZ Integer ring holding numbers at least as large as the squares of numbers in Z.
/// @tparam ulongZ The unsigned integer domain of longZ.
/// @tparam MOD The modulus of the class —needs to be coprime with r = 2^logsize = 2^{number of bits in uZ}.
///         The modulus of each Montgomery object can be changed in runtime.
template <typename Z, typename uZ, typename longZ, typename ulongZ, uZ MOD = 0>
class Montgomery {
    using M = Montgomery;
protected:
    static constexpr Z logsize = sizeof(uZ) * 8;
//  uZ r = 2^{logsize}
    uZ _mod;   // The, perhaps changed, modulus of the current Montgomery object.
    uZ modinv; // m^{-1} mod r.
    uZ val;    // the value of X in the (X * r) mod m representation.

    /// @return (x * r^{-1}) mod m.
    uZ reduce(ulongZ const &x) const {
        uZ q = uZ(x) * modinv;               // q = xn' mod r
        uZ m = (ulongZ(q) * mod()) >> logsize; // m = qn
        return (x >> logsize) + mod() - m;
    }

    /// @return (x * r) mod m.
    uZ transform(uZ const &x) const {
        return (ulongZ(x) << logsize) % mod();
    }

public:
    /// @return The template parameter MOD.
    static uZ mod()   { return MOD; }
    /// @return The, perhaps changed, mod of this Montgomery class.
    uZ get_mod() { return _mod; }

    void set_mod(uZ m) {
        _mod = m;
        modinv = _mod;
        assert(m % 2 == 1);
        while(modinv * _mod != 1) modinv *= uZ(2) - _mod * modinv;
    }

    Montgomery() : val(0) {
        set_mod(MOD);
    }

    Montgomery(longZ val_) : Montgomery() {
        val = transform(val_ % _mod + _mod);
    }

    template <typename T>
    Montgomery(T val_): Montgomery(static_cast<longZ>(val_)) { }

    M &operator+=(M const &b) {
        val += b.val;
        while(val >= _mod) val -= _mod;
        return *this;
    }

    M &operator-=(M const &b) {
        val += _mod - b.val;
        while(val >= _mod) val -= _mod;
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


    // If A, B have the same modulus then equality can be checked faster from their values "inside" the space.
    friend bool operator==(const M &A, const M &B) { return A.get() == B.get(); }
    friend bool operator!=(const M &A, const M &B) { return A.get() != B.get(); }

    friend std::ostream &operator<<(std::ostream &os, const M &x) { return os << x.get(); }
    friend std::istream &operator>>(std::istream &is, M &x) {
        Z t; is >> t;
        x = M(longZ(t));
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
        Z x = get(), y = _mod, u = 1, v = 0;
        while (y > 0) {
            Z t = x / y;
            std::swap(x -= t * y, y);
            std::swap(u -= t * v, v);
        }
        return M(u);
    }

    operator int() const { return get(); }
    /// @return The "real" value of the object, i.e. extract X from (X * r) mod m.
    uZ get() const {
        uZ ret = reduce(val);
        while(ret >= _mod) ret -= _mod;
        return ret;
    }

};


template <typename Z, typename uZ, typename longZ, typename ulongZ, uZ mod>
struct is_modular<zeno::Montgomery<Z, uZ, longZ, ulongZ, mod>> : std::true_type {};

template<uint32_t mod = 0>
using Montgomery32 = Montgomery<int32_t, uint32_t, int64_t, uint64_t, mod>;

template<uint64_t mod = 0>
using Montgomery64 = Montgomery<int64_t, uint64_t, __int128_t, __uint128_t, mod>;



} // namespace zeno
