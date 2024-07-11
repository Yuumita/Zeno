#ifndef ZENO_PRIMALITY_HPP
#define ZENO_PRIMALITY_HPP

#include <vector>
#include <cstdint>
#include <stdlib.h> /* abs */
#include "structures.hpp"

namespace zeno 
{

namespace primality
{

template<class Z = int64_t>
Z find_factor(Z n);


std::vector<int> primes;
// O(1) prime check and O(logA) factorization in [1, N] with O(N) precomputation with spf (Smallest prime factor)
const int MAX_ESIEVE_NUM = 2e6;
int spf[MAX_ESIEVE_NUM + 1];
bool has_precomputed_esieve = false;
void prime_sieve() {
    for(int i = 1; i <= MAX_ESIEVE_NUM; i++) spf[i] = i; 
    for(int i = 2; i <= MAX_ESIEVE_NUM; i += 2) spf[i] = 2;
    primes.push_back(2);
    for(int i = 3; i <= MAX_ESIEVE_NUM; i += 2) if(spf[i] == i) {
        primes.push_back(i);
        for(long long j = (long long)i*i; j <= MAX_ESIEVE_NUM; j += i){
            if(spf[j] == j)
                spf[j] = i;
        }
    }
    has_precomputed_esieve = true;
}

// Pollard's rho polynomial g(x) = (x^2 + c) mod m
template<class Z = int64_t>
Z _G(Z x, Z c, Z m) {
    return ((__int128_t) x * x % m + c) % m;
} 

template<class Z = int64_t>
Z pollard_rho(Z N, Z x0 = 2, Z c = 1) {
    Z x = x0, y = x0;
    Z g = 1;
    while (g == 1) { // floyd's cycle detection
        x = _G(x, c, N);
        y = _G(_G(y, c, N), c, N);
        g = zeno::gcd((x < y ? y - x : x - y), N);
    }
    return g;
}

template<class Z = int64_t>
Z find_factor(Z n) {
    if(n <= MAX_ESIEVE_NUM) 
        return spf[n];
    Z d = pollard_rho(n, 2, 1);
    if(d == n)
        d = pollard_rho(n, 1, 2);

    if(d == n) {
        // TODO: check if n is prime, else find another factor
        for(d = 2; d*d <= n; d++) {
            if(n % d == 0) 
                return d;
        }
        return n;
    }
    return d; 
}

template<class Z = int64_t>
std::vector<Z> factorize(Z n) {

    if(!has_precomputed_esieve) prime_sieve();

    std::vector<Z> ret;
    do {
        Z d = find_factor(n);
        ret.push_back(d);
        n /= d;
    } while(n > 1);

    return ret;
}

template<class Z = int64_t>
void standard_factorize(Z n, std::vector<Z> &p, std::vector<Z> &v) {
    std::vector<Z> &f = factorize(n);
    if(f.size() == 0) return;
    for(int i = 0, j; i < f.size(); i++) {

        for(j = i; j+1 < f.size() && f[i] == f[j+1]; j++);

        p.push_back(f[i]);
        v.push_back(j - i + 1);
    }
    return;
}

} // namespace primality


}; // namespace zeno

#endif /* ZENO_PRIMALITY_HPP */