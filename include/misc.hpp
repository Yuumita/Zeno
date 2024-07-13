#ifndef ZENO_MISC_HPP
#define ZENO_MISC_HPP


#include <cstdint>

namespace zeno 
{ 


template<class G = int64_t, class Z = int64_t>
G pow(G g, Z n) {
    return binary_exponentiation(g, n);
}



} // namespace zeno

#endif /* ZENO_MISC_HPP */