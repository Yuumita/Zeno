#include <bits/stdc++.h>



#include "../include/modular.hpp"
#include "../include/fps.hpp"

using namespace std;
#define rep(i, N) for (long long i = 0; i < (long long)(N); i++)

template <typename T>
void test(int n) {

    using FPS = zeno::FormalPowerSeries<T>;
    
    vector<T> a(n), b(n), x(n), y(n);

    rep(i, n) {
        a[i] = rand() & 0xFFFF;
        b[i] = rand() & 0xFFFF;
        x[i] = rand() & 0xFFFF;
        y[i] = rand() & 0xFFFF;
    }

    // inv
    {
        FPS f = FPS(a);
        FPS g = f.inv();
        assert((f * g) == FPS(T(1)));
    }

    // exp/log
    {
        FPS f = FPS(a);
        FPS g = f.exp();
        g = g.log();
        assert(f == g);
    }

    // multipoint evaluation
    {
        /// TODO: ...
    }

    // interpolation
    {
        FPS f = FPS::interpolate(x, y);
        for(int i = 0; i < n; i++) {
            assert(f.eval(x[i]) == y[i]);
        }
    }
}


int main() {
    srand(time(0));
    for(int N = 1; N <= 200; N++) { // n = logN <= 10
        test<long long>(N);
        test<zeno::modint998244353>(N);
    }
}
