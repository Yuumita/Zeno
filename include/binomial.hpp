#pragma once

#include "convolution.hpp"
#include <vector>
#include <algorithm>

namespace zeno {

/// TODO: perhaps support negative binomial coefficients etc.
/// @brief The Binomial class used for computing factorials and binomial coefficients.
template <typename Tp>
class Binomial {
private:
    static std::vector<Tp> _factorial;
    static std::vector<Tp> _ifactorial;

    /// @brief Precomputes the factorial and inverse factorial up to N.
    static void precompute(size_t N = 0) {
        size_t n = _factorial.size();
        if(N == 0) N = n * 2;
        if(N <= n) return;


        _factorial.resize(N + 1), _ifactorial.resize(N + 1);

        for(size_t i = n; i <= N; ++i) 
            _factorial[i] = _factorial[i - 1] * Tp(i);

        _ifactorial[N] = Tp(1) / _factorial[N];
        for(size_t i = N - 1; i >= n; --i)
            _ifactorial[i] = _ifactorial[i + 1] * Tp(i + 1);

    }

public:

    Binomial(size_t n = 0) {
        if(n) precompute();
    }

    /// @return The factorial of n, n!.
    static Tp factorial(int n) {
        if(n < 0) return Tp(0);
        if(n >= _factorial.size()) precompute();
        return _factorial[n];
    }

    /// @return The inverse of the factorial of n, 1 / n!.
    static Tp ifactorial(int n) {
        if(n < 0) return Tp(0);
        if(n >= _ifactorial.size()) precompute();
        return _ifactorial[n];
    }

    /// @return The binomial coefficient C(n, k) = n! / (k! * (n - k)!).
    static Tp binom(int n, int k) {
        if(n < 0 || k < 0 || n < k) return Tp(0); 
        return factorial(n) * ifactorial(k) * ifactorial(n - k);
    }
    inline Tp operator()(int n, int m) { return binom(n, m); }

    /// TODO: the following
    // Tp int_binom(int n, int m) const {
    //     if(n < 0 || m < 0 || n < m) return Tp(0); 
    // }

    /// @return 1 / n.
    static Tp inv(int n) {
        if(n <= 0) return n;
        return factorial(n - 1) * ifactorial(n);
    }

    /// @return The falling factorial n! / (n - k)! = n * (n - 1) * ... * (n - k + 1).
    static Tp falling(int n, int k) {
        if(k < 0) return Tp(0);
        return factorial(n) * ifactorial(n - k);
    }

    /// @return The rising factorial (n + k - 1)! / (n - 1)! = n * (n + 1) * ... * (n + k - 1).
    static Tp rising(int n, int k) {
        if(k < 0) return Tp(0);
        return factorial(n + k - 1) * ifactorial(n - 1);
    }

    /// @return The multinomial coefficient (m1 + ... + mk) / (m1! * .... * mk!).
    template <typename Z, std::enable_if_t<!std::is_integral<Z>::value>* = nullptr>
    static Tp multinomial(std::vector<Z> const &m) {
        Z n = Z(0);
        for(Z &k: m) {
            if(k < 0) return Tp(0);
            n += k;
        }
        Tp ret = factorial(n);
        for(Z &k: m)
            ret *= ifactorial(k);
        return ret;
    }
    template <typename Z, std::enable_if_t<!std::is_integral<Z>::value>* = nullptr>
    inline Tp operator()(std::vector<Z> const &m) { return multinomial(m); }
};

template <typename Tp>
std::vector<Tp> Binomial<Tp>::_factorial = {Tp(1)};

template <typename Tp>
std::vector<Tp> Binomial<Tp>::_ifactorial = {Tp(1)};


} // namespace zeno
