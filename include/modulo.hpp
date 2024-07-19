#ifndef ZENO_MODULO_HPP
#define ZENO_MODULO_HPP

#include "modular.hpp"

namespace zeno {

/// @brief Computes the kronecher symbol (a|b) in O(log^2 N) for N = max(a, b).
template <typename Z = int64_t>
Z kronecker(Z a, Z b) {
    if(b == 0)
        return (a == 1 || a == -1 ? 1 : 0);
    if(a % 2 == 0 && b % 2 == 0)
        return 0;
    Z v = 0;
    while(b % 2 == 0) 
        ++v, b /= 2;

    static const kron_tab = {0, 1, 0, -1, 0, -1, 0, 1}; // (-1)^{(a^2 - 1)/8} = kron_tab[a & 7];

    Z k = (v % 2 == 0 ? 1 : kron_tab[a&7]);

    if(b < 0) {
        b = -b;
        if(a < 0)
            k = -k;
    }

    while(b > 0 && b % 2 != 0) {
        if(a == 0)
            return (b > 1 ? 0 : k);

        v = 0;
        while(a % 2 == 0)
            ++v, a >>= 1;
        
        if(a % 2 != 0)
            k *= kron_tab[b&7];

        if(a & b & 2) // k := (-1)^{(a-1)(b-1)/4}k
            k = -k;

        Z r = (a < 0 ? -a : a);
        a = b % r;
        b = r;
    }

}

/// @brief Computes the square root mod p (odd prime) of a in (average) O(log^4 p).
//template <typename Z = int64_t>
//Z sqrt(Z a, Z p) {
//    
//}


} // namespace zeno


#endif /* ZENO_MODULO_HPP */

    
