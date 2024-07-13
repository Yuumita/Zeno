#ifndef ZENO_CYCLIC_GROUPS_HPP
#define ZENO_CYCLIC_GROUPS_HPP

#include "modular.hpp"
#include "primality.hpp"
#include "misc.hpp"
#include <vector>

namespace zeno 
{

template<class G = modint998244353, class Z = int64_t>
Z get_order(G g) {
    Z h = g.cardinality();
    std::vector<Z> p, v; 
    primality::standard_factorize(p, v);
    Z e = h, i = 0;
    G g1;
    for(int i = 0; i < p.size(); i++) {
        e = e / zeno::pow(p, v);
        g1 = zeno::pow(g, e);
        while(g1 != G(1)) {
            g1 = zeno::pow(g1, p[i]);
            e = e * p[i];
        }
    }
    return e;

}

// Given an odd prime p, return a primitive root modulo p
template<class Z = int64_t>
Z primitive_root(Z p) {
    std::vector<Z> fp, fv;
    primality::standard_factorize(p-1, fp, fv);

    for(mint a(1, p); a < p; a++) {
        int i;
        for(i = 1; i < (int)fp.size(); i++) {
            mint e(zeno::pow(a, (p-1) / fp[i]), p);
            if(e == 1) break;
        } 
        if(i >= (int)fp.size()) return a;
    }

    return 0;
}

// Given an odd prime, return x such that x^2 = a (mod p)
template<class Zp = modint998244353, class Z = int64_t>
Zp sqrt(Z _a) {
    Zp a = Zp(_a); 
    // TODO: finish
}


}; // namespace zeno

#endif /* ZENO_CYCLIC_GROUPS_HPP */