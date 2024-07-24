#pragma once

#include <vector>
namespace zeno
{
    


namespace internal {

/// @brief Garner's algorithm
const int m0 = 167772161;
const int m1 = 469762049;
const int m2 = 754974721;
using mint0 = modint<m0>;
using mint1 = modint<m1>;
using mint2 = modint<m2>;
const int r01 = mint1(m0).inv().val();
const int r02 = mint2(m0).inv().val();
const int r12 = mint2(m1).inv().val();


} // namespace internal


using int128 = __int128_t;
template <typename T>
std::vector<int128> convolution_int128(std::vector<T> const &a, std::vector<T> const &b) {
    if (std::min(a.size(), b.size()) <= fft::magic_number) 
        return convolution_naive<int64_t>(a, b);

    int 

}


} // namespace zeno