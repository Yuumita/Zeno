#ifndef INTERNAL_ZENO_BIT_HPP
#define INTERNAL_ZENO_BIT_HPP

namespace zeno {

namespace internal {

    /**
    * @param n A non-negative integer
    * @return ceil(log(n)) (i.e minimum `h` s.t. `n <= 2**h`)
    */
    int ceil_log2(int n) {
        int x = 0;
        while((1u << x) < (unsigned int)(n)) x++;    
        return x;
    }

} // namespace internal

} // namespace zeno

#endif // INTERNAL_ZENO_BIT_HPP 