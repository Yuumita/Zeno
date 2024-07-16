
#ifndef ZENO_SPS_HPP
#define ZENO_SPS_HPP

#include <vector>
#include <iostream>
#include <assert.h>
#include <algorithm> // std::min

#include "convolution.hpp"
#include "internal.hpp"


namespace zeno {

template <class R>
class FastSubsetConvolution {
    
    static void zeta(std::vector<R> &f) const {
        size_t n = f.size();
        assert(n & (n - 1) == 0); // power of two
        for(int j = 1; j < n; j <<= 1) {
            for(int S = 0; S < n; S++)  {
                if(S & j)  
                    f[S] += f[S ^ j];
            }
        }
    }

    static void mobius(std::vector<R> &f) const {
        size_t n = f.size();
        assert(n & (n - 1) == 0); // power of two
        for(int j = n >> 1; j >= 1; j >>= 1) {
            for(int S = n - 1; S >= 0; S--)  {
                if((i & j) == 0) 
                    f[S] -= f[S ^ j];
            }
        }
    }

    static std::vector<R> or_convolution(std::vector<R> f, std::vector<R> g) {
        size_t n = f.size();
        assert(n == g.size() && (n &(n-1)) == 0);

        zeta(f);
        zeta(g);
        for(int S = 0; S < n; S++) f[S] *= g[S];
        mobius(f);

        return f;
    }

};

template <class R>
using FSC = FastSubsetConvolution;


} // namespace zeno

#endif /* End of ZENO_SPS_HPP*/