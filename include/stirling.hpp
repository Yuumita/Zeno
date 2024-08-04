#pragma once

#include "fps.hpp"
#include "convolution.hpp"
#include <vector>
#include <algorithm>

namespace zeno {

/// TODO: check edge cases (e.g. S(0, 0) = 1, s(0, 0) = 1)
/// @brief The Stirling class used for computing the stirling numbers.
template <typename Tp>
class Stirling {

    using FPS = zeno::FormalPowerSeries;

private:
    static std::vector<std::vector<Tp>> _stirling1;
    static std::vector<std::vector<Tp>> _stirling2;

    /// @brief Precompute the stirling numbers of the first kind s(N, *) for the given N.
    /// @ref https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind#Definitions
    static void precompute_stirling1(size_t N) {
        if(_stirling1.size() <= N) _stirling1.resize(N + 1, {});
        if(_stirling1[N].size() > 0) return;
        if(N == 0) { _stirling1[N] = {1}; return; }

        precompute_stirling1(N >> 1);

        /// falling(x, n) = sum_{k=0}^n s(n, k) x^k
        /// falling(x, n) = falling(x, n/2) * falling(x, n/2) if n is even etc.
        FPS<Tp> P = FPS<Tp>(_stirling1[N >> 1]);
        P *= P.taylor_shift(Tp(-(N >> 1)));
        if(N & 1) P *= FPS<Tp>({-Tp(N - 1), 1});

        _stirling1[N].resize(N + 1);
        for(size_t i = 0; i <= N; i++)
            _stirling1[N][i] = P[i];
    }

    /// @brief Precompute the stirling numbers of the second kind S(N, *) for the given N
    /// @ref https://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind#Explicit_formula
    static void precompute_stirling2(size_t N) {
        if(_stirling2.size() <= N) _stirling2.resize(N + 1, {});
        if(_stirling2[N].size() > 0) return;

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

    Stirling() {}

    /// @return s(n, k) where s(*, *) is the stirling number of the first kind.
    /// @ref https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind
    static Tp stirling1(int n, int k) {
        if(n < 0 || k < 0 || k > n) return 0;
        precompute_stirling1(n);
        return _stirling1[n][k];

    }

    /// @return S(n, k) where S(*, *) is the stirling number of the second kind.
    /// @ref https://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind
    static Tp stirling2(int n, int k) {
        if(n < 0 || k < 0 || k > n) return 0;
        precompute_stirling2(n);
        return _stirling2[n][k];
    }

};

template <typename Tp>
std::vector<std::vector<Tp>> Stirling<Tp>::_stirling1 = {{Tp(1)}};

template <typename Tp>
std::vector<std::vector<Tp>> Stirling<Tp>::_stirling2 = {{Tp(1)}};


} // namespace zeno
