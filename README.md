# Zeno

*Zeno* is a C++ template library for computational algebra.

## Overview

*Zeno* is a C++ template library for computational algebra which implements many modern algorithms and practices for formal power series, arbitrary precision arithmetic, modular arithmetic and other algebraic structures. Zeno is designed to be efficient and flexible, with a clear focus on simplicity and elegancy.

## Key Features

- Formal power series ($R[[x]]$).
- Sparse multivariate polynomials ($R[x_{1},\dots,x_{n}]$).
- Set power series operations ($R[[x_{1},\dots,x_{n}]]/\left< x_{1}^{2}-q(x_{1}),\dots,x_{n}^{2}-q(x_{n}) \right>$ where $q \in R[x]$ with $\mathrm{deg}(q)\leq 1$).
- Arbitrary precision arithmetic (efficient for integers $\leq 10^{10^6}$).
- Vector space operations on bitsets ($\mathbb{F}_{2}^d$).
- Misc.: Ring of integers modulo $m$ ($\mathbb{Z}/m\mathbb{Z}$). convolution via FFT, NTT etc., elementary combinatorics, elementary number theory, field of fractions ($\mathrm{Frac}(R)$, e.g. $\mathbb{Q}=\mathrm{Frac}(\mathbb{Z})$), complex numbers ($\mathbb{C}$).


## Some Examples

```cpp
#include <iostream>
#include "include/modular.hpp" // import modint32
#include "include/fps.hpp"     // import FormalPowerSeries
using mint = zeno::modint32<998244353>;    // mint is now the Z_{998244353} field
using FPS = zeno::FormalPowerSeries<mint>; // FPS is now the Z_p[[x]] ring where p = 998244353

int main() {
    FPS f(3, 1);             // 1 + 2x + 3x^2
    FPS g = f * FPS({1, 1}); // (1 + 2x + 3x^2) * (1 + x)
    std::cout << g.eval(1) << std::endl; // evaluate g(x) at x = 1
    std::cout << (f + g).eval(1) << std::endl; // evaluate (f + g)(x) at x = 1

    int N = 4;
    FPS h = f.inv(N);    // h = 1 / f mod x^N
    assert((h * f).mod_xk(N) == FPS(1)); // (h * f) mod x^N == 1 --> true

    f[0] = 0; // [x^0]f(x) is set to 0 (so that exp(f) is well defined in the field of mint)
    FPS exp_f = f.exp(N); // exp_f = {1 + f + f^2/2 + f^3/6}, i.e. f.exp(N) = exp(f) mod x^N
    f[0] = 1; // [x^0]f(x) is set to 1 (so that log(f) is well defined in the field of mint)
    FPS log_f = f.log(N); // (1 + x).log(N) = log(1 + x) mod x^N

    std::vector<mint> Y = f.eval({1, 2, 3, 4, 5}); // evaluate f at x = 1, ..., 5

    mint c = 8;
    FPS fts = f.taylor_shift(c); // fts(x) = f(x + c).
}
```


```cpp
#include <iostream>
#include "include/bigint.hpp"

using MPI = zeno::MultiPrecisionInteger;

int main() {

    MPI a(123456789);
    MPI b; std::cin >> b; // e.g. b = 98765432143243204378204578

    MPI c("-123456789123456789123456789");

    // Basic Arithmetic Operations
    MPI addition       = a + b;
    MPI subtraction    = a - b;
    MPI multiplication = a * b;
    MPI division       = c / a; 
    MPI remainder      = c % a;

    // Increment and Decrement
    ++a, a--;

    // Comparison Operations and outputs
    std::cout << "a + b: " << a + b << std::endl;
    std::cout << "Is a equal to b? " << (a == b ? "Yes" : "No") << std::endl;
    std::cout << "Is a less than b? " << (a < b ? "Yes" : "No") << std::endl;

    return 0;
}
```
