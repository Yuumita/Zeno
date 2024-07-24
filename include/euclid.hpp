#pragma once

#include <vector>
#include <cstdint>

namespace zeno {


// Here, Z is any Euclid domain

template<typename Z = int64_t>
Z gcd(Z a, Z b) { /* Euclid's algorithm */
    return (b == Z(0) ? a : gcd(b, a % b));
}

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
Z _gcd(Z a, Z b) { /* Euclid's algorithm */
    static_assert(std::is_integral_v<Z>, "Template parameter Z should be integral for binary_gcd");
    return binary_gcd(a, b);
}

// Extended Euclid's algorithm, returns gcd(a, b) and sets x, y to those s.t. ax + by = gcd(a, b)
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


/// @brief Safely (i.e. in a way that avoid overflows) compute (a * b) mod n
template<typename Z = int64_t>
Z safe_mul_mod(Z a, Z b, Z n) {
    bool sign = true;
    for(Z &x: {a, b}) {
        if(x < 0) {
            x = -x;
            sign = !sign;
        }
    }

    a %= n, b %= n;
    Z ret = 0;
    while(b > 0) {
        if(b % 2) {
            ret = (ret + a) % n;
        }
        a = (a + a) % n, b /= 2;
    }

    return (sign ? ret : (-ret) % n);
}

/// @brief Compute a solution to the diophantine equation ax + by = c. 
/// @return true iff ax + by = c has at least one solution.
template<typename Z = int64_t>
bool solve_diophantine(Z a, Z b, Z c, Z &x, Z &y, Z &g) {
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
    x = dx + safe_mul_mod(x, c / g, b);
    y = dy + safe_mul_mod(y, c / g, a);
    g = (g < 0 ? -g: g);
    return true;
}


template<typename Z = int64_t>
Z inductive_chinese_remainder_theorem(std::vector<Z> const &m, std::vector<Z> const x) {
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