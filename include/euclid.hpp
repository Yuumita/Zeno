#ifndef ZENO_EUCLID_HPP
#define ZENO_EUCLID_HPP

#include <vector>

namespace zeno {

namespace internal {

template<typename Z = int64_t>
Z inv_gcd(Z a, Z b) { /* Euclid's algorithm */
    return (b == 0 ? a : gcd(b, a % b));
}
    
} // namespace internal
    

// Here, Z is any Euclid domain

template<typename Z = int64_t>
Z gcd(Z a, Z b) { /* Euclid's algorithm */
    return (b == 0 ? a : gcd(b, a % b));
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
    return binary_gcd(a, b);
}

// Extended Euclid's algorithm, returns gcd(a, b) and sets x, y to those s.t. ax + by = gcd(a, b)
template<typename Z = int64_t>
Z extended_gcd(Z a, Z b, Z &x, Z &y) {
    if (b == 0) {
        x = 1, y = 0;
        return a;
    }
    Z x1, y1;
    Z d = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}

// fast implementation of extended_gcd
template<typename Z = int64_t>
Z extended_binary_gcd(Z a, Z b, Z &x, Z &y) {
    Z y1, y3, t1, t3;
    bool f2 = false;
    if(a < b) { 
        std::swap(a, b);
        std::swap(x, y);
    }
    if(b == 0) {
        x = 1, y = 0;
        return a;
    }
    Z q = a / b, r = a % b;
    a = b;
    b = r;
    if(b == 0) {
        x = 1, y = 0;
        return a;
    }

    Z k = 0;
    while((a & 1) && (b & 1))
        k += 1, a >>= 1, b >>= 1;
    
    if(~b & 1) {
        std::swap(a, b);
        std::swap(x, y);
    }

    Z d = a;
    x = 1, y1 = y3 = b;

    if(a & 1) 
        t1 = 0, t3 = -b;
    else {

    }
    
    // TODO: complete
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


#endif /* ZENO_EUCLID_HPP */