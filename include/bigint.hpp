#ifndef ZENO_BIGINT_HPP
#define ZENO_BIGINT_HPP

#include <vector>
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <string>

namespace zeno {


/// TODO: change 'int' to 'uint32_t'
class MultiPrecisionInteger {
    using MPI = MultiPrecisionInteger;

private:

    static const int B = 1000'000'000; // base
    static const int logB = 9;

    bool neg;
    std::vector<int> digits;


public:

    MultiPrecisionInteger(): neg(false), digits({}) {};

    MultiPrecisionInteger(bool neg_, std::vector<int> const &digits_)
        : neg(neg_), digits(digits_) {};

    template <typename Z>
    MultiPrecisionInteger(Z x) : neg(false) {
        static_assert(std::is_integral_v<Z>, "Template parameter Z is not integral.");
        if (std::is_signed_v<Z>)
            if (x < Z(0)) neg = true, x = -x;
        while (x > Z(0))  dat.push_back(x % B), x /= B;
    }

    template<typename Iterator>
    MultiPrecisionInteger(Iterator begin, Iterator end) : neg(false) {
        using T = typename std::iterator_traits<Iterator>::value_type;
        static_assert(std::is_same<T, char>::value, "Iterator value type must be of type char");
        if(*begin == '+' || *begin == '-') neg = (*begin == '-'), ++begin;
        while(end - begin > 0) {
            int seg_len = std::min(end - begin, logB);
            Iterator seg_begin = std::prev(end, seg_len);

            int x = 0;
            for(Iterator it = std::prev(end, seg_len); it != end; ++it)
                x = 10 * x + (*it - '0');
            
            digits.push_back(x);
            end = seg_begin;
        }
        _shrink();
    }

    template<typename Iterator>
    MultiPrecisionInteger(Iterator begin) : neg(false) {
        using T = typename std::iterator_traits<Iterator>::value_type;
        static_assert(std::is_same<T, char>::value, "Iterator value type must be of type char");
        Iterator it = begin;
        while(it && *it != '\0') ++it;
        MultiPrecisionInteger(begin, end);
    }

    MultiPrecisionInteger(const std::string &s) : neg(false) {
        MultiPrecisionInteger(s.begin(), s.end());
    }

    MPI& operator+=(const MPI& x) { 

    }

    MPI& operator-=(const MPI& rhs) { 

    }

    MPI& operator*=(const MPI& rhs) { 

    }

    MPI& operator/=(const MPI& rhs) { 

    }

    MPI& operator%=(const MPI& rhs) { 

    }



    MPI& operator++() { 

    }
    MPI operator++(int) { MPI ret(*this); operator++(); return ret; }

    MPI& operator--() { 

    }
    MPI operator--(int) { MPI ret(*this); operator--(); return ret; }

    MPI operator+() const { return MPI(*this); }
    MPI operator-() const { return MPI(!neg, digits); }

    // TODO: check for equivalences between representations
    friend bool operator==(const MPI &A, const MPI &B) { 
        return A.neg == B.neg && A.digits == B.digits; 
    }
    friend bool operator!=(const MPI &A, const MPI &B) { 
        return A.neg != B.neg || A.digits != B.digits;
    }
    friend bool operator<(const MPI& lhs, const MPI& rhs) {
        return _cmp(lhs, rhs) < 0;
    }
    friend bool operator<=(const MPI& lhs, const MPI& rhs) {
        return _cmp(lhs, rhs) <= 0;
    }
    friend bool operator>(const MPI& lhs, const MPI& rhs) {
        return _cmp(lhs, rhs) > 0;
    }
    friend bool operator>=(const MPI& lhs, const MPI& rhs) {
        return _cmp(lhs, rhs) >= 0;
    }

private:

    void _shrink() {
        while(!digits.empty() && digits.back() == 0) digits.pop_back();
    }

    static void _shrink(std::vector<int> &a) {
        while(!a.empty() && a.back() == 0) a.pop_back();
    }

    /// @return -1, 0, +1 if lhs is <, ==, > rhs respectively.
    static int _cmp(const MPI &lhs, const MPI &rhs) {
        if(lhs.neg != rhs.neg)
            return (lhs.neg ? +1 : -1);
        std::vector<int> const &l = lhs.digits, &r = rhs.digits;
        if(l.size() != r.size()) 
            return (l.size() < r.size() ? -1: +1);
        for(int i = l.size() - 1; i >= 0; i--) 
            if(l[i] != r[i]) return (l[i] < r[i] ? -1: +1);
        return 0;
    }

    /// @brief Add the unsigned MPI vectors (a + b) and store the result to a
    static void _add(std::vector<int> &a, std::vector<int> const &b) {
        int carry = 0;
        for(int i = 0; i < a.size() || i < b.size() || carry; i++) {
            if(i == a.size()) a.push_back(0);
            a[i] += carry + (i < b.size() ? b[i] : 0);
            carry = (a[i] >= B);
            if(carry) a[i] -= B;
        }
    }

    /// @brief Subtract the two unsigned MPI vectors (a - b) and store the result to a
    /// @pre a >= b
    static void _sub(std::vector<int> &a, std::vector<int> const &b) {
        int carry = 0;
        for(int i = 0; i < b.size() || carry; i++) {
            a[i] -= carry + (i < b.size() ? b[i] : 0);
            carry = (a[i] < 0);
            if(carry) a[i] += B;
        }
        _shrink(a);
    }
    
    
    static const int magic_number = 60;

    /// @brief Multiply the unsigned MPI vector a with b and store the result to a
    /// @pre b < B
    static void _mul(std::vector<int> &a, std::vector<int> const &b) {
        // Ideally a >= b
        if(a.size() <= magic_number && b.size() <= magic_number) {

        }
    }
};

using MPI    = MultiPrecisionInteger;
using BigInt = MultiPrecisionInteger;
using bigint = MultiPrecisionInteger;
    
} // namespace zeno
#endif /* ZENO_BIGINT_HPP */