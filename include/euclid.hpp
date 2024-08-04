#pragma once

#include <vector>
#include <cstdint>
#include "internal.hpp"

namespace zeno {


/// @brief Euclid's gcd algorithm.
/// @tparam Z An integral domain.
/// @return gcd(a, b)
template<typename Z = int64_t>
Z gcd(Z a, Z b) { 
    return (b == Z(0) ? a : gcd(b, a % b));
}

/// @tparam Z An integral domain.
/// @return lcm(a, b)
template<typename Z = int64_t>
Z lcm(Z a, Z b) {
    return (a / gcd(a, b)) * b;
}

/// @brief Fast implementation of Euclid's gcd algorithm.
/// @tparam Z An integral domain.
/// @return gcd(a, b)
template<typename Z = int64_t>
Z binary_gcd(Z a, Z b) { /* fast implementation of Euclid's algorithm */
    if(a < b) std::swap(a, b);
    if(b == 0) 
        return a;

    Z r = a % b; 
    a = b; b = r; /* one euclid step */

    if(b == 0) 
        return a;

    Z k = 0;
    while((a & 1) && (b & 1))
        k += 1, a >>= 1, b >>= 1;

    while(a & 1) a >>= 1;
    while(b & 1) b >>= 1;

    do {
        Z t = (a - b) >> 1;
        if(t == 0) 
            return a << k;
            
        while(t & 1) t >>= 1;

        if(t > 0) 
            a = t;
        else
            b = -t;

    } while(b != 0);

    return a << k;
}

template<typename Z = int64_t>
Z _gcd(Z a, Z b) {
    if constexpr (std::is_integral_v<Z>)
        return binary_gcd(a, b);
    return gcd(a,b);
}

template<typename Z = int64_t>
Z _lcm(Z a, Z b) {
    return (a / _gcd(a, b)) * b;
}

/// @brief Extended Euclid's algorithm. x, y are set to those s.t. ax + by = gcd(a, b)
/// @tparam Z An integral domain.
/// @return gcd(a, b)
template <typename Z = int64_t>
Z extended_gcd(Z a, Z b, Z &x, Z &y) {
    if (b == 0) {
        x = 1, y = 0;
        return a;
    }
    Z x1, y1;
    Z g = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return g;
}

// fast implementation of extended_gcd
//template<typename Z = int64_t>
//Z extended_binary_gcd(Z a, Z b, Z &x, Z &y) {
//    // TODO: complete
//}




/// @brief Compute a solution to the linear diophantine equation ax + by = c. If a solution
///        exists then the variables x, y are set to one solution and g to gcd(a, b).
/// @tparam Z An integral domain.
/// @return true iff ax + by = c has at least one solution.
template<typename Z = int64_t>
bool solve_linear_diophantine(Z a, Z b, Z c, Z &x, Z &y, Z &g) {
    if (a == 0 && b == 0) {
        if (c == 0) {
            x = y = g = 0;
            return true;
        }
        return false;
    }

    if (a == 0) {
        if (c % b == 0) {
            x = 0;
            y = c / b;
            g = (b < 0 ? -b : b);
            return true;
        }
        return false;
    }

    if (b == 0) {
        if (c % a == 0) {
            x = c / a;
            y = 0;
            g = (a < 0 ? -a : a);
            return true;
        }
        return false;
    }

    g = extended_gcd(a, b, x, y);
    if (c % g != 0)
        return false;
    Z dx = c / a;
    c -= dx * a;
    Z dy = c / b;
    c -= dy * b;
    x = dx + internal::safe_mul_mod(x, c / g, b);
    y = dy + internal::safe_mul_mod(y, c / g, a);
    g = (g < 0 ? -g: g);
    return true;
}


/// @brief Apply the inductive CRT to solve for the number X such that X = x[i] (mod m[i]). 
///        The moduli m[i] should be pairwise coprime.
/// @ref https://en.wikipedia.org/wiki/Chinese_remainder_theorem#Existence_(constructive_proof)
/// @tparam Z An integral domain.
/// @return An integer X such that X = x[i] (mod m[i])
template<typename Z = int64_t>
Z chinese_remainder_theorem(std::vector<Z> const &m, std::vector<Z> const x) {
    assert(m.size() == x.size());

    Z M = m[0], X = x[0];
    Z u, v, d;
    for(Z i = 1; i < (Z)m.size(); i++) {
        d = extended_gcd(M, m[i], u, v);
        assert(d == 1);

        X = u * M * x[i] + v * m[i] * X;
        M *= m[i];
        X %= M;
    }

    return x;
}

} // namespace zeno