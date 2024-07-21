#ifndef ZENO_MFPS_HPP
#define ZENO_MFPS_HPP

#include <stdlib.h>
#include <array>

namespace zeno {

/// @brief Implementation of the ring R[[x_0, ...,x_{N-1}]] / <x_0^n[0], ..., x_{N-1}^n[N-1]>
template<class R, size_t N>
class mfps {
private:
    fps F;
    std::array<int, N> n;
public:
    mfps();
};
    
    
} // namespace zeno


#endif /* ZENO_MFPS_HPP */