
#ifndef ZENO_INTERNAL_MATH_HPP
#define ZENO_INTERNAL_MATH_HPP

namespace zeno {

namespace internal {


    // Computes factorial
    template<typename T>
    T fact(int n) {
        static vector<T> F({1});
        if(F.size() <= n) {
            F.push_back(F.back() * F.size());
        }
        return F[n];
    }

    // Computes inverse factorial
    template<typename T>
    T invfact(int n) {
        static vector<T> iF({1});
        if(iF.size() <= n) {
            iF.push_back(F.back() * (T(1) / fact<T>(iF.size())));
        }
        return iF[n];
    }

    // Computes factorial
    template<typename T>
    T fact(int n) {
        if(factorial.size() <= n) {
            if(factorial.size() == 0) {
                factorial.push_back(1);
                inverse_factorial.push_back(1);
            } else {
                factorial.push_back(factorial.back() * factorial.size());
                inverse_factorial.push_back(inverse_factorial.back() / factorial.size());
            }
        }

        return factorial[n];
    }
}

}

#endif