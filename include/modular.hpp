#pragma once

#include "primality.hpp"
#include <cassert>

namespace zeno 
{


/// @brief Static (i.e. modulus given at compile time) implementation of the Z/mZ ring.
/// @tparam Z Integer ring holding numbers at least as large as MOD.
/// @tparam longZ Integer ring holding numbers at least as large as the square of a number in Z.
/// @tparam MOD The modulo m of the Z/mZ ring.
template <typename Z, Z MOD, typename longZ = Z>
class static_modular {
    static_assert(1 <= MOD, "Template parameter MOD must be greater than or equal to 1.");
    using modular = static_modular<Z, MOD, longZ>;
    using uZ = std::make_unsigned_t<Z>;
    using longuZ = std::make_unsigned_t<longZ>;
private: 
    uZ v;
    static constexpr uZ umod() { return MOD; }
    static constexpr bool is_prime = internal::is_prime_constexpr<MOD>;
public:
    static constexpr Z mod() { return MOD; }
    static Z cardinality() { return mod(); }
    static Z size() { return mod(); }

    static_modular(): v(0) {}
    template <typename T, std::enable_if_t<std::is_integral<T>::value>* = nullptr>
    static_modular(T vv) {
        longZ x = static_cast<longZ>(vv % static_cast<longZ>(mod()));
        if (x < 0) x += mod();
        v = static_cast<uZ>(x);
    }

    template <typename T, std::enable_if_t<!std::is_integral<T>::value>* = nullptr>
    static_modular(T vv): static_modular(static_cast<longZ>(vv)) { }

    operator int() const { return v; }
    uZ val() { return v; }

    modular& operator+=(const modular &x) { v += x.v; if(v >= umod()) v -= umod(); return *this; }
    modular& operator-=(const modular &x) { if(v < x.v) v += umod() - x.v; else v -= x.v; return *this; }
    modular& operator*=(const modular &x) { v = (uZ)((longuZ)(v) * x.v % umod()); return *this; }
    modular& operator/=(const modular &x) { return operator*=(x.inv()); }

    modular& operator++() { if(v == umod() - 1) v = 0; else v++; return *this; }
    modular operator++(int) { modular ret(*this); operator++(); return ret; }

    modular& operator--() { if(v == 0) v = umod() - 1; else v--; return *this; }
    modular operator--(int) { modular ret(*this); operator--(); return ret; }

    modular operator+() const { return *this; }
    modular operator-() const { return modular(0) - *this; }

    //modular operator-(modular x) const { return modular(*this) -= x; }


    /* casting (maybe avoid it ???) */
    template<typename T, std::enable_if_t<std::is_integral_v<T>>* = nullptr> 
    modular& operator+=(const T &x) { return operator+=(modular(x)); }
    
    template<typename T, std::enable_if_t<std::is_integral_v<T>>* = nullptr> 
    modular& operator-=(const T &x) { return operator-=(modular(x)); }
    
    template<typename T, std::enable_if_t<std::is_integral_v<T>>* = nullptr> 
    modular& operator*=(const T &x) { return operator*=(modular(x)); }
    
    template<typename T, std::enable_if_t<std::is_integral_v<T>>* = nullptr> 
    modular& operator/=(const T &x) { return operator/=(modular(x)); }
    

    template<typename T, std::enable_if_t<std::is_integral_v<T>>* = nullptr> 
    friend modular operator+(const modular &M, const T &x) { return modular(M) += x; }

    template<typename T, std::enable_if_t<std::is_integral_v<T>>* = nullptr> 
    friend modular operator-(const modular &M, const T &x) { return modular(M) -= x; }

    template<typename T, std::enable_if_t<std::is_integral_v<T>>* = nullptr> 
    friend modular operator*(const modular &M, const T &x) { return modular(M) *= x; }

    template<typename T, std::enable_if_t<std::is_integral_v<T>>* = nullptr> 
    friend modular operator/(const modular &M, const T &x) { return modular(M) /= x; }

    template<typename T> friend bool operator==(const modular &A, const T &B) { return A.v == modular(B).v; }
    template<typename T> friend bool operator!=(const modular &A, const T &B) { return A.v != modular(B).v; }

    friend modular operator+(const modular &lhs, const modular &rhs) { return modular(lhs) += rhs; }
    friend modular operator-(const modular &lhs, const modular &rhs) { return modular(lhs) -= rhs; }
    friend modular operator*(const modular &lhs, const modular &rhs) { return modular(lhs) *= rhs; }
    friend modular operator/(const modular &lhs, const modular &rhs) { return modular(lhs) /= rhs; }

    friend bool operator==(const modular &A, const modular &B) { return A.v == B.v; }
    friend bool operator!=(const modular &A, const modular &B) { return A.v != B.v; }


    modular pow(Z n) const { /* change it so it uses zeno::pow ??? */
        modular ret(1), mul(v);
        if(n < 0) mul = mul.inv(), n = -n;
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul, n >>= 1;
        }
        return ret;
    }

    modular inv() const {
        // Fermat's Little Theorem
        if(is_prime)
            return pow(umod() - 2);

        Z a = v, b = mod(), u = 1, v = 0, t;
        while (b > 0) {
            t = a / b;
            std::swap(a -= t * b, b);
            std::swap(u -= t * v, v);
        }
        assert(a*u + b*v == 1); // (b*v == 0) TODO: test this 
        return modular(u);
    }
    modular inverse() const { return inv(); }

    friend std::ostream &operator<<(std::ostream &os, const modular &p) { return os << p.v; }
    friend std::istream &operator>>(std::istream &is, modular &a) {
        Z t; is >> t;
        a = modular(t);
        return (is);
    }
};


// /// @brief Dynamic (modulo given at runtime time) implementation of the Z/mZ ring.
// /// @tparam Z Integer ring holding numbers at least as large as MOD.
// /// @tparam longZ Integer ring holding numbers at least as large as the square of a number in Z.
// template <typename Z, typename longZ = int64_t>
// class dynamic_modular {
//     static_assert(std::is_integral<Z>::value, "Template parameter Z must be integral.");
//     using modular = dynamic_modular;
//     using uZ = std::make_unsigned_t<Z>;
//     using longuZ = std::make_unsigned_t<longZ>;
// private: 
//     uZ v;
//     Z MOD;
//     uZ umod() { return MOD; }
//     // static constexpr is_prime = internal::is_prime_constexpr<MOD> ???
// public:
//     Z mod() { return MOD; }
//     void set_mod(int m) { 
//         assert(1 <= m);
//         MOD = m;
//     }
//     Z cardinality() { return mod(); }
//     Z size() { return mod(); }
// 
//     dynamic_modular(): v(0), MOD(2) {}
//     dynamic_modular(Z m) : MOD(m) {}
//     template <typename T, std::enable_if_t<std::is_integral<T>::value>* = nullptr>
//     dynamic_modular(T vv, Z m) : MOD(m) {
//         longZ x = static_cast<longZ>(vv % static_cast<longZ>(mod()));
//         if (x < 0) x += mod();
//         v = static_cast<uZ>(x);
//     }
//     template <typename T, std::enable_if_t<!std::is_integral<T>::value>* = nullptr>
//     dynamic_modular(T vv): dynamic_modular(static_cast<longZ>(vv)) { }
// 
//     operator int() const { return v; }
//     uZ val() { return v; }
// 
//     modular& operator+=(const modular &x) { v += x.v; if(v >= umod()) v -= umod(); return *this; }
//     modular& operator-=(const modular &x) { if(v < x.v) v += umod() - x.v; else v -= x.v; return *this; }
//     modular& operator*=(const modular &x) { v = (uZ)((longuZ)(v) * x.v )% umod(); return *this; }
//     modular& operator/=(const modular &x) { return operator*=(x.inv()); }
// 
//     modular& operator++() { if(v == umod() - 1) v = 0; else v++; return *this; }
//     modular operator++(int) { modular ret(*this); operator++(); return ret; }
// 
//     modular& operator--() { if(v == 0) v = umod() - 1; else v--; return *this; }
//     modular operator--(int) { modular ret(*this); operator--(); return ret; }
// 
//     modular operator+() const { return *this; }
//     modular operator-() const { return modular(0) - *this; }
// 
//     //modular operator-(modular x) const { return modular(*this) -= x; }
// 
// 
//     /* casting (maybe avoid it ???) */
//     template<typename T> modular& operator+=(const T &x) { return operator+=(modular(x)); }
//     template<typename T> modular& operator-=(const T &x) { return operator-=(modular(x)); }
//     template<typename T> modular& operator*=(const T &x) { return operator*=(modular(x)); }
//     template<typename T> modular& operator/=(const T &x) { return operator/=(modular(x)); }
// 
//     template<typename T, std::enable_if_t<std::is_integral_v<T>>* = nullptr> 
//     friend modular operator+(const modular &M, const T &x) { return modular(M) += x; }
// 
//     template<typename T, std::enable_if_t<std::is_integral_v<T>>* = nullptr> 
//     friend modular operator-(const modular &M, const T &x) { return modular(M) -= x; }
// 
//     template<typename T, std::enable_if_t<std::is_integral_v<T>>* = nullptr> 
//     friend modular operator*(const modular &M, const T &x) { return modular(M) *= x; }
// 
//     template<typename T, std::enable_if_t<std::is_integral_v<T>>* = nullptr> 
//     friend modular operator/(const modular &M, const T &x) { return modular(M) /= x; }
// 
//     friend bool operator==(const modular &A, const modular &B) { return A.v == B.v; }
//     friend bool operator!=(const modular &A, const modular &B) { return A.v != B.v; }
// 
// 
//     modular pow(Z n) const { /* change it so it uses zeno::pow ??? */
//         modular ret(1), mul(v);
//         if(n < 0) mul = mul.inv(), n = -n;
//         while (n > 0) {
//             if (n & 1) ret *= mul;
//             mul *= mul, n >>= 1;
//         }
//         return ret;
//     }
// 
//     modular inv() const {
//         assert(v != 0);
//         
//         // TODO: check dynamically if prime
// 
//         Z a = v, b = mod(), u = 1, v = 0, t;
//         while (b > 0) {
//             t = a / b;
//             std::swap(a -= t * b, b);
//             std::swap(u -= t * v, v);
//         }
//         assert(a*u + b*v == 1); // (b*v == 0) TODO: test this 
//         return modular(u);
//     }
// 
//     modular inverse() const { return inv(); }
// 
//     friend std::ostream &operator<<(std::ostream &os, const modular &p) { return os << p.v; }
//     friend std::istream &operator>>(std::istream &is, modular &a) {
//         Z t; is >> t;
//         a = modular(t);
//         return (is);
//     }
// };


template<typename T>
struct is_modular : std::false_type {};

template<typename Z, Z m, typename longZ>
struct is_modular<zeno::static_modular<Z, m, longZ>> : std::true_type {};

// template<typename Z, typename longZ>
// struct is_modular<zeno::dynamic_modular<Z, longZ>> : std::true_type {};

template<typename T>
inline constexpr bool is_modular_v = is_modular<T>::value;


template<int M>
using modint = static_modular<int, M, long long>; 

template<int64_t M>
using modlong = static_modular<int64_t, M, __int128_t>; 

template<int M>
using modint32 = static_modular<int32_t, M, int64_t>; 

template<int M>
using modint64 = static_modular<int64_t, M, __int128_t>; 

using modint998244353 = modint<998244353>; 

// using dmint = dynamic_modular<int32_t>;
// using dmint32 = dynamic_modular<int32_t>;
// using dmint64 = dynamic_modular<int64_t>;

}; // namespace zeno