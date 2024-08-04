#pragma once

#include "modular.hpp"
#include "primality.hpp"
#include "internal.hpp"
#include <vector>

namespace zeno 
{

template<class G, class Z = int64_t>
Z get_order(G g) {
    Z h = g.cardinality();
    std::vector<Z> p, v; 
    primality::standard_factorize(p, v);
    Z e = h, i = 0;
    G g1;
    for(int i = 0; i < p.size(); i++) {
        e = e / internal::pow(p, v);
        g1 = internal::pow(g, e);
        while(g1 != G(1)) {
            g1 = internal::pow(g1, p[i]);
            e = e * p[i];
        }
    }
    return e;
}


/// @brief Given an odd prime p, return a primitive root modulo p.
template<class Z = int64_t>
Z primitive_root(Z p) {
    std::vector<Z> fp, fv;
    primality::standard_factorize(p-1, fp, fv);

    for(Z i = 2; i < p; i++) {
        bool ok = true;
        for(auto &f: fp) 
            if(internal::binary_exponentiation_mod<Z>(i, (p - 1) / f, p) == 1) {
                ok = false;
                break;
            }
        if(ok) return i;
    }

    return 0;
}

}; // namespace zeno