#include <bits/stdc++.h>

#include "../include/modular.hpp"
#include "../include/polynomial.hpp"
#include "../include/combinatorics.hpp"

using namespace std;
#define rep(i, N) for (long long i = 0; i < (long long)(N); i++)

template <typename T>
void test(int n) {
}


template <class R, size_t N>
using SMPoly = zeno::SparseMultivariatePolynomial<R, N>;

using mint = zeno::modint998244353;

template <class R, size_t N>
using mmial = zeno::MononomialTerm<R, N>;

int main() {
    srand(time(0));
    

    SMPoly<mint, 2> P;
    P += mmial<mint, 2>(mint(2), array<int, 2>{1, 1});
    P += mmial<mint, 2>(mint(13), array<int, 2>{0, 2});
    P += mmial<mint, 2>(mint(5), array<int, 2>{1, 1});

    cout << P << endl;

}
