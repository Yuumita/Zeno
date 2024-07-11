#ifndef ZENO_MISC_HPP
#define ZENO_MISC_HPP


#include <cstdint>

namespace zeno 
{ 

/* G should form a group, Z should be an integer struct */
template<class G = int64_t, class Z = int64_t>
G binary_exponentiation(G g, Z n) {
    if(n < 0) g = 1/g, n = -n;

    G x = g, y = G(1);
    while(n > 0) {
        if(n % 2 == 1) y *= x;
        n /= 2, x *= x;
    }

    return y
}

template<class G = int64_t, class Z = int64_t>
G pow(G g, Z n) {
    return binary_exponentiation(g, n);
}



} // namespace zeno

#endif /* ZENO_MISC_HPP */