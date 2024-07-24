#pragma once

#include <vector>
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <iostream>
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
        while (x > Z(0))  digits.push_back(x % B), x /= B;
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


    friend MPI operator+(const MPI& lhs, const MPI& rhs) {
        if (lhs.neg == rhs.neg) 
            return MPI(lhs.neg, _add(lhs.digits, rhs.digits));
        std::vector<int> c;
        bool nneg;
        if (_cmp(lhs.digits, rhs.digits) < 0) { // |l| <= |r|
            c = _sub(rhs.digits, lhs.digits);
            nneg = MPI::is_zero(c) ? false : rhs.neg;
        } else {
            c = _sub(lhs.digits, rhs.digits);
            nneg = MPI::is_zero(c) ? false : lhs.neg;
        }
        return MPI(nneg, std::move(c));
    }

    friend MPI operator-(const MPI& lhs, const MPI& rhs) { return lhs + (-rhs); }

    friend MPI operator*(const MPI& lhs, const MPI& rhs) {
        std::vector<int> c = _mul(lhs.digits, rhs.digits);
        bool nneg = MPI::is_zero(c) ? false : (lhs.neg ^ rhs.neg);
        return MPI(nneg, std::move(c));
    }

    friend std::pair<MPI, MPI> divmod(const MPI& lhs, const MPI& rhs) {
        auto dm = _divmod_newton(lhs.digits, rhs.digits);
        bool dn = _is_zero(dm.first) ? false : lhs.neg != rhs.neg;
        bool mn = _is_zero(dm.second) ? false : lhs.neg;
        return {MPI(dn, std::move(dm.first)), MPI(mn, std::move(dm.second))};
    }

    friend MPI operator/(const MPI& lhs, const MPI& rhs) {
        return divmod(lhs, rhs).first;
    }

    friend MPI operator%(const MPI& lhs, const MPI& rhs) {
        return divmod(lhs, rhs).second;
    }

    MPI& operator+=(const MPI& rhs) { return (*this) = (*this) + rhs; }
    MPI& operator-=(const MPI& rhs) { return (*this) = (*this) - rhs; }
    MPI& operator*=(const MPI& rhs) { return (*this) = (*this) * rhs; }
    MPI& operator/=(const MPI& rhs) { return (*this) = (*this) / rhs; }
    MPI& operator%=(const MPI& rhs) { return (*this) = (*this) % rhs; }


    MPI& operator++() { return (*this) = (*this) + MPI(1); }
    MPI operator++(int) { MPI ret(*this); operator++(); return ret; }

    MPI& operator--() { return (*this) = (*this) - MPI(1); }
    MPI operator--(int) { MPI ret(*this); operator--(); return ret; }

    MPI operator+() const { return MPI(*this); }
    MPI operator-() const { return MPI((digits.empty() ? false: !neg), digits); }

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
    static int _cmp(std::vector<int> const &l, std::vector<int> const &r) {
        if(l.size() != r.size()) 
            return (l.size() < r.size() ? -1: +1);
        for(int i = l.size() - 1; i >= 0; i--) 
            if(l[i] != r[i]) return (l[i] < r[i] ? -1: +1);
        return 0;
    }

    /// @return -1, 0, +1 if lhs is <, ==, > rhs respectively.
    static int _cmp(const MPI &lhs, const MPI &rhs) {
        if(lhs.neg != rhs.neg)
            return (lhs.neg ? +1 : -1);
        const std::vector<int> &l = lhs.digits, &r = rhs.digits;
        return _cmp(l, r);
    }


    /// @brief Add the unsigned MPI vectors (a + b) and store the result to a
    static void _add_op(std::vector<int> &a, std::vector<int> const &b) {
        int carry = 0;
        for(int i = 0; i < a.size() || i < b.size() || carry; i++) {
            if(i == a.size()) a.push_back(0);
            a[i] += carry + (i < b.size() ? b[i] : 0);
            carry = (a[i] >= B);
            if(carry) a[i] -= B;
        }
    }

    static std::vector<int> _add(std::vector<int> const &a, std::vector<int> const &b) {
        std::vector<int> ret = a;
        _add_op(ret, b);
        return ret;
    }

    /// @brief Subtract the two unsigned MPI vectors (a - b) and store the result to a
    /// @pre a >= b
    static void _sub_op(std::vector<int> &a, std::vector<int> const &b) {
        int carry = 0;
        for(int i = 0; i < b.size() || carry; i++) {
            a[i] -= carry + (i < b.size() ? b[i] : 0);
            carry = (a[i] < 0);
            if(carry) a[i] += B;
        }
        _shrink(a);
    }

    static std::vector<int> _sub(std::vector<int> const &a, std::vector<int> const &b) {
        std::vector<int> ret = a;
        _sub_op(ret, b);
        return ret;
    }

    static const int magic_number = 60;
    static const long long BB = (long long)(B) * B;

    static std::vector<int> _mul(std::vector<int> const &a, int const &b0) {
        std::vector<int> c(a.size() + 1);
        long long carry = 0;
        for(int i = 0; i < a.size() || carry; ++i) {
            carry += (long long)(a[i]) * b0;
            c[i]   = (int)(carry % B);
            carry /= B;
        }
        _shrink(c);
    }

    /// @return The multiplication, a * b, result of the (unsigned) MPIs a, b.
    static std::vector<int> _mul(std::vector<int> const &a, std::vector<int> const &b) {
        if(a.empty() || b.empty()) return {};
        if(b.size() == 1) return _mul(a, b[0]);
        if(a.size() <= magic_number || b.size() <= magic_number) { // naive convolution
            std::vector<long long> c(a.size() + b.size());
            for(int i = 0; i < a.size(); i++) {
                for(int j = 0; j < b.size(); j++) {
                    c[i + j] += (long long)(a[i]) * (long long)b[j]; 
                    // a[i] * b[i] < 10^18 => c[i + j] < 10^18 + 10^
                    if(c[i + j] >= BB) {
                        c[i + j] -= BB;
                        c[i + j + 1]  += (long long)(B);
                    }
                }
            }
            std::vector<int> ret(c.size() + 3);
            long long carry = 0;
            int i;
            for(i = 0; i < c.size(); i++) {
                carry += c[i];
                ret[i] = carry % B;
                carry /= B;
            }
            while(carry) {
                ret[i++] = carry % B;
                carry /= B;
            }
            _shrink(ret);
            return ret;
        } 

        // ntt convolution 
        /// TODO: check sizes for overflows
        std::vector<__int128_t> c = BigNTT::multiply_int128(a, b);
        std::vector<int> ret(c.size() + 3);
        __int128_t carry = 0;
        int i;
        for(i = 0; i < c.size(); i++) {
            ret.push_back(carry % B);
            carry /= B;
        }
        while(carry) {
            ret[i++] = carry % B;
            carry /= B;
        }
        _shrink(ret);
        return ret;
    }
    
    
    /// @brief Multiply the unsigned MPI vector a with b and store the result to a
    /// @pre b < B
    static void _mul_op(std::vector<int> &a, std::vector<int> const &b) {
        a = _mul(a, b);
    }

    static std::pair<std::vector<int>, std::vector<int>> _div(std::vector<int> const &a, int const &b0) {
        if(b0 == 0) {
            std::cerr << "Division by zero" << std::endl;
            exit(1);
        }
        if(b0 == 1) return {a, {}};
        std::vector<int> q(a.size()), r;
        long long carry = 0;
        for(int i = a.size() - 1; i >= 0; --i) {
            carry = carry * B + a[i];
            q[i]  = (int)(carry / b0);
            carry %= b0;
        }
        _shrink(q);
        if(carry) r.push_back(carry);
        return {q, r};
    }

    /// @ref Knuth's ACP vol. 2, Seminumerical Algorithms
    static std::pair<std::vector<int>, std::vector<int>> 
    _div_long(std::vector<int> const &a, std::vector<int> const &b) {
        if(MPI::is_zero(b)) {
            std::cerr << "Division by zero" << std::endl;
            exit(1);
        }
        if(b.size() == 1) return _div(a, b[0]);
        if(_cmp(a, b) < 0) return {{}, a};

        int K = (B - 1) / b.back(); // normalization [a//b = (Ka)//(Kb) = u//v]
        std::vector<int> u = _mul(a, K);
        std::vector<int> v = _mul(b, K);
        // now v.back() >= B/2 so that 2 * v.back() >= B

        std::vector<int> q(u.size() - v.size() + 1, 0);
        std::vector<int> r(u.end() - v.size(), u.end());
        for(int i = q.size() - 1; i >= 0; --i) {
            if(r.size() > v.size()) {
                int qq = ((long long)(r[r.size() - 1]) * B + r[r.size() - 2]) / v.back();
                // qq - 2 <= Q <= qq where Q = r // v (from normalization)
                std::vector<int> vqq = _mul(v, qq);

                while(_cmp(r, vqq) < 0) qq--, vqq = _sub(vqq, v);
                r = _sub(r, vqq);
                while(_cmp(v, r) <= 0) qq++, r = _sub(r, v);
                q[i] = qq;
            } else if(r.size() == v.size()) { // from normalization qq is then at most 1
                if(_cmp(v, r) < 0) { 
                    q[i] = 1, r = _sub(r, v);
                }
            }
            if(i > 0) r.insert(r.begin(), u[i - 1]);
        }
        _shrink(q), _shrink(r);
        auto [qT, rT] = _div(r, K); // denormalization
        return {q, qT};
    }

    /// @ref https://en.wikipedia.org/wiki/Division_algorithm#Newton–Raphson_division
    static std::pair<std::vector<int>, std::vector<int>> 
    _div_newton_raphson (std::vector<int> const &a, std::vector<int> const &b) {
        if(MPI::is_zero(b)) {
            std::cerr << "Division by zero" << std::endl;
            exit(1);
        }
        if(b.size() <= magic_number) return _div_long(a, b);
        if(_cmp(a, b) < 0) return {{}, a};

        int K = (B - 1) / b.back(); // normalization [a//b = (Ka)//(Kb) = u//v]
        std::vector<int> u = _mul(a, K);
        std::vector<int> v = _mul(b, K);
        // now v.back() >= B/2 so that 2 * v.back() >= B

        std::vector<int> q(u.size() - v.size() + 1, 0);
        std::vector<int> r(u.end() - v.size(), u.end());
        for(int i = q.size() - 1; i >= 0; --i) {
            if(r.size() > v.size()) {
                int qq = ((long long)(r[r.size() - 1]) * B + r[r.size() - 2]) / v.back();
                // qq - 2 <= Q <= qq where Q = r // v (from normalization)
                std::vector<int> vqq = _mul(v, qq);

                while(_cmp(r, vqq) < 0) qq--, vqq = _sub(vqq, v);
                r = _sub(r, vqq);
                while(_cmp(v, r) <= 0) qq++, r = _sub(r, v);
                q[i] = qq;
            } else if(r.size() == v.size()) { // from normalization qq is then at most 1
                if(_cmp(v, r) < 0) { 
                    q[i] = 1, r = _sub(r, v);
                }
            }
            if(i > 0) r.insert(r.begin(), u[i - 1]);
        }
        _shrink(q), _shrink(r);
        auto [qT, rT] = _div(r, K); // denormalization
        return {q, qT};
    }

public:

    static MPI abs(MPI const &m) { return MPI(false, m.digits); }
    bool is_zero() const { return digits.empty(); }
    static bool is_zero(std::vector<int> const &a) { return a.empty(); }

    static std::string to_string(int value, bool padding = false) {
        std::string ret = std::to_string(value);
        if(padding) {
            std::string p;
            for(int i = ret.size(); i < 9; i++) p.push_back('0');
            return p + ret;
        }
        return ret;
    }

    std::string to_string() const {
        if(is_zero()) return "0";
        std::string ret;
        if(neg) ret.push_back('-');
        for(int i = digits.size() - 1; i >= 0; --i)
            ret += MPI::to_string(digits[i], i != digits.size() - 1);
        return ret;
    }

    friend std::istream& operator>>(std::istream& is, MPI &m) {
        std::string s;
        is >> s;
        m = MPI(s);
        return is;
    }
    friend std::ostream& operator<<(std::ostream& os, MPI const &m) {
        return os << m.to_string();
    }
};

using MPI    = MultiPrecisionInteger;
using BigInt = MultiPrecisionInteger;
using bigint = MultiPrecisionInteger;
    
} // namespace zeno