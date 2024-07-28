#include <bits/stdc++.h>


#include "../include/modular.hpp"
#include "../include/sps.hpp"

using namespace std;
#define rep(i, N) for (long long i = 0; i < (long long)(N); i++)

template <typename T>
void test(int n) {
    
    using SPS = zeno::SetPowerSeries<T>;

    assert((n & (n - 1)) == 0);
    vector<T> a(n), b(n), c(n), d;

    rep(i, n) {
        a[i] = rand() & 0xFFFF;
        b[i] = rand() & 0xFFFF;
    }

    {
        d = a;
        SPS::zeta(d);
        SPS::mobius(d);
        assert(a == d);
    }

    {
        d = a;
        SPS::super_zeta(d);
        SPS::super_mobius(d);
        assert(a == d);
    }

    {
        d = a;
        SPS::walsh_hadamard_transform(d);
        SPS::walsh_hadamard_transform(d, true);
        assert(a == d);
    }

    // and convolution
    {
        d = SPS::and_convolution(a, b);
        c.assign(n, 0);
        rep(i, n) rep(j, n) c[i & j] += a[i] * b[j];
        assert(c == d);
    }

    // or convolution
    {
        d = SPS::or_convolution(a, b);
        c.assign(n, 0);
        rep(i, n) rep(j, n) c[i | j] += a[i] * b[j];
        assert(c == d);
    }

    // xor convolution
    {
        d = SPS::xor_convolution(a, b);
        c.assign(n, 0);
        rep(i, n) rep(j, n) c[i ^ j] += a[i] * b[j];
        assert(c == d);
    }
}

int main() {
    srand(time(0));
    for(int N = 1; N <= 4096; N *= 2) { // n = logN <= 10
        test<long long>(N);
        test<zeno::modint998244353>(N);
    }
}
