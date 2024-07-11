#ifndef ZENO_INTERNAL_HPP
#define ZENO_INTERNAL_HPP


#include <cstdint>

namespace zeno
{

namespace internal
{

const int DEFAULT_MOD = 998244353;

// Computes factorial
template<typename Z = int64_t>
Z fact(int n) {
    static std::vector<Z> _factorial({1});
    while(_factorial.size() <= n)
        _factorial.push_back(_factorial.back() * _factorial.size());
    return _factorial[n];
}


int64_t ceil_log2(int64_t n) {
    int64_t x = 0;
    while((1u << x) < (uint64_t)(n)) x++;    
    return x;
}

}; // namespace internal

}; // namespace zeno


#endif /* ZENO_INTERNAL_HPP */