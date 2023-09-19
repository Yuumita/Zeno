
#ifndef ZENO_INTERNAL_MATH_HPP
#define ZENO_INTERNAL_MATH_HPP


#include <vector>

namespace zeno {

namespace internal {


    // Computes factorial
    template<typename T>
    T fact(int n) {
        static std::vector<T> F({1});
        if(F.size() <= n) {
            F.push_back(F.back() * F.size());
        }
        return F[n];
    }

    // Computes inverse factorial
    template<typename T>
    T invfact(int n) {
        static std::vector<T> iF({1});
        if(iF.size() <= n) {
            iF.push_back(iF.back() * (T(1) / fact<T>(iF.size())));
        }
        return iF[n];
    }

}

}

#endif