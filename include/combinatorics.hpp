#pragma once

#include "convolution.hpp"
#include <vector>
#include <algorithm>

namespace zeno {

template <typename Tp>
class Binomial {
private:
    static std::vector<Tp> _factorial;
    static std::vector<Tp> _ifactorial;

    static void precompute(size_t N = 0) {
        size_t n = _factorial.size();
        if(N <= n) return;
        if(N == 0) N = n * 2;

        _factorial.resize(N), _ifactorial.resize(N);

        for(size_t i = n; i < N; ++i) 
            _factorial[i] = _factorial[i - 1] * Tp(i);

        _ifactorial.back() = Tp(1) / _factorial.back();
        for(size_t i = N - 2; i >= n; --i)
            _ifactorial[i] = _ifactorial[i + 1] * Tp(i + 1);

    }

public:

    Binomial(size_t n = 0) {
        if(n) precompute();
    }

    static Tp factorial(int n) {
        if(n < 0) return Tp(0);
        if(n >= _factorial.size()) precompute();
        return _factorial[n];
    }

    static Tp ifactorial(int n) {
        if(n < 0) return Tp(0);
        if(n >= _ifactorial.size()) precompute();
        return _ifactorial[n];
    }

    static Tp binom(int n, int m) {
        if(n < 0 || m < 0 || n < m) return Tp(0); 
        return factorial(n) * ifactorial(m) * ifactorial(n - m);
    }

    inline Tp operator()(int n, int m) { return binom(n, m); }

    /// TODO: the following
    // Tp int_binom(int n, int m) const {
    //     if(n < 0 || m < 0 || n < m) return Tp(0); 
    // }

    static Tp inv(int n) {
        if(n <= 0) return n;
        return factorial(n - 1) * ifactorial(n);
    }

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


/// TODO: check edge cases (e.g. S(0, 0) = 1, s(0, 0) = 1)
/// @brief
template <typename Tp>
class Stirling {
private:
    static std::vector<std::vector<Tp>> _stirling1;
    static std::vector<std::vector<Tp>> _stirling2;

    /// @brief Precompute the stirling numbers of the second kind S(N, *) for the given N
    /// @ref https://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind#Explicit_formula
    static void precompute_stirling2(size_t N) {
        if(_stirling2.size() <= N) _stirling2.resize(N, {});
        if(_stirling2[N].size() >= N + 1) return;

        // S(N, *) = convolution(a, b) where a[i] = (-1)^i / i!, b[i] = i^N / i!
        std::vector<Tp> a(N + 1), b(N + 1);
        for(size_t i = 0; i <= N; i++) {
            a[i] = Binomial<Tp>::ifactorial(i) * (i % 2 == 0 ? +1 : -1);
            b[i] = internal::pow<Tp>(i, N) * Binomial<Tp>::ifactorial(i);
        }
        _stirling2[N] = zeno::convolution(a, b);
        _stirling2[N].resize(N + 1);
    }

public:

    static Tp stirling1(int n, int k) {

    }

    static Tp stirling2(int n, int k) {
        if(n < 0 || k < 0 || k > n) return 0;
        if(_stirling2.size() <= n)
            precompute_stirling2(n);
        return _stirling2[n][k];
    }


};

template <typename Tp>
std::vector<std::vector<Tp>> Stirling<Tp>::_stirling1 = {{Tp(1)}};

template <typename Tp>
std::vector<std::vector<Tp>> Stirling<Tp>::_stirling2 = {{Tp(1)}};


} // namespace zeno
