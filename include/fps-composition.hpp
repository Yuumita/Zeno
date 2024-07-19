#ifndef ZENO_FPS_COMPOSITION_HPP
#define ZENO_FPS_COMPOSITION_HPP

#include "fps.hpp"

namespace zeno {



/// @brief Computes f(G(x)) mod x^m in O(nlog^2n).
template <class R>
FPS<R> composition(const FPS<R> &f, const FPS<R> &G, int m) {
    if(m <= 0) return FPS({});
    if(G.is_zero())
        return FPS(f.get(0));
    
}



} // namespace zeno


#endif /* ZENO_FPS_COMPOSITION_HPP */