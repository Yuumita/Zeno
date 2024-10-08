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


    // simple operations
    {
        FPS f = FPS(a), g = FPS(b);

        assert(f + g - f == g);
        assert((f - f).is_zero());

        assert(((f + g) - f) + g == T(2)*g);

        vector<T> c = zeno::convolution_naive(a, b);
        assert( (f*g) == FPS(c));
    }

    // inv
    {
        FPS f = FPS(a).mod_xk(n+1);
        if(f[0] == 0) f[0] = T(1);
        FPS g = f.inv(n+1);
        assert((f * g).mod_xk(n+1) == FPS(T(1)));
    }

    // exp/log
    {
        FPS f = FPS(a);
        f[0] = T(0);
        assert(f.exp(n+1).log(n+1) == f.mod_xk(n+1));
    }

    // pow
    {
        FPS f = FPS(a), g = FPS(1);
        int k = rand() & 0xFF;
        for(int i = 0; i < k; i++) g = (f*g).mod_xk(n+1);
        assert(f.pow(k, n+1) == g);
    }

    // borel transform
    {
        FPS f = FPS(a);
        assert(f.borel().inv_borel() == f);
    }

    // Taylor shift
    {
        FPS f = FPS(a), g = FPS(a);
        T c = T(rand() & 0xFF);
        f = f.taylor_shift(c);
        for(int i = 0; i < n; i++) 
            assert(f.eval(T(i)) == g.eval(T(i) + c));
    }

    // multipoint evaluation
    {
        FPS f = FPS(a);
        vector<T> v = f.eval(x);
        for(int i = 0; i < n; i++)
            assert(v[i] == f.eval(x[i]));
            
    }

    bool ok = true;
    for(int i = 0; i < n; i++) 
        for(int j = i + 1; j < n; j++) 
            if(x[i] == x[j]) ok = false;
    if(ok) // interpolation
    {
        FPS f = FPS::interpolate(x, y);
        for(int i = 0; i < n; i++)
            assert(f.eval(x[i]) == y[i]);
    }
}


int main() {
    srand(time(0));
    zeno::FPS<zeno::modint998244353> f;
    f = {0, 0, 0, 1};
    cout << f << " -> " << f.taylor_shift(3) << endl;;
    
    for(int N = 200; N <= 800; N *= 2) { // n = logN <= 10
        cerr << "N = " << N << endl;
        test<zeno::modint998244353>(N);
    }
}
