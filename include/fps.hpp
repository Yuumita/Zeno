#ifndef ZENO_FPS_HPP
#define ZENO_FPS_HPP

#include <vector>
#include <iostream>
#include <assert.h>
#include <algorithm> // std::min

#include "convolution.hpp"
#include "internal.hpp"


namespace zeno {

/// @tparam R The type of coefficients
template<typename R>
class FormalPowerSeries : public std::vector<R> {
    using FPS = FormalPowerSeries;
public:

    /// @brief Reduce the 0 coefficients of the FPS (i.e remove all trailing zeroes)
    void normalize(){
        while(!this->empty() && this->back() == R(0))
            this->pop_back();
    }

    /// @brief Constructs a new FPS object with degree -1 (no coefficients) 
    FormalPowerSeries(): (*this)({}) {}

    /// @brief Constructs a new FPS object with coefficients from the array coefs 
    FormalPowerSeries(const std::vector<R> &coefs): a(coefs) {
        this->normalize();
    }

    /// @brief Constructs a new (constant, i.e of 0 degree) FPS object which is just v
    FormalPowerSeries(const R &v) {
        this->assign(1, v);
        this->normalize();
        normalize();
    }

    /// @brief Constructss a new FPS object of degree n-1 with initial values v 
    FormalPowerSeries(const int n, const R &v) {
        this->resize(n);
        for(int i = 0; i < n; i++) a[i] = v;
        normalize();
    }

    FormalPowerSeries(std::vector<R>::iterator begin, std::vector<R>::iterator end) {
        this->assign(begin, end);
        normalize();
    }

    /// @brief Returns the degree of the FPS (-1 if the FPS is zero)
    int deg() const {
        return this->size() - 1;
    }

    /// @brief Returns the order of the FPS (-1 if the FPS is zero)
    int ord() const {
        for(int i = 0; i < a.size(); i++) 
            if(this->at(i) != R(0)) return i;
        return -1;
    }

    int is_zero() const {
        return this->empty();
    }

    /// @brief Returns a copy of the coefficient of x^i
    R get(const int i) const {
        return (0 <= i && i <= deg() ? this->at(i) : R(0));
    }

    /// @brief Sets the coefficient of x^i
    void set(const int i, const R c) {
        if(c == R(0)) return;
        assert(0 <= i);
        if(this->size() <= i) this->resize(i+1, R(0));
        this->at(i) = c;
    }

    /// @brief Returns a copy of the leading coefficient
    R leadcoef() const {
        if(this->empty()) return 0;
        return this->back();
    }


    /// @brief Addition of two FPSs (*this and P)
    FPS plus(const FPS &P) const {
        return *this + P;
    }

    /// @brief Multiplication of two FPSs (*this and P)
    FPS times(const FPS &P) const {
        return *this * P;
    }

    /// @brief Multiplication of this FPS with the scalar v
    FPS times(const R &v) const {
        return *this * v;
    }

    /// @brief Subtraction of two FPSs (*this - P)
    FPS minus(const FPS &P) const {
        return *this - P;
    }

    FPS& operator+=(const FPS &rhs) {
        int d = (this->deg() > rhs.deg() ? this->deg() : rhs.deg());
        this->resize(d);
        for(int i = 0; i < (int)rhs.size(); i++){
            this->at(i) += rhs[i];
        }
        this->normalize();
        return *this;
    }

    FPS& operator-=(const FPS &rhs) {
        int d = (this->deg() > rhs.deg() ? this->deg() : rhs.deg());
        this->resize(d + 1);
        for(int i = 0; i < (int)rhs.size(); i++){
            this->at(i) -= rhs[i];
        }
        this->normalize();
        return *this;
    }

    FPS& operator*=(const FPS &rhs) {
        this->a = fft::convolution<R>(this->a, rhs.a);
        this->normalize();
        return *this;
    }

    FPS& operator*=(const R &rhs) {
        for(auto &e: this) e *= rhs;
        this->normalize();
        return *this;
    }

    FPS& operator/=(const FPS &rhs) {
        return *this = div(rhs).first;
    }

    FPS& operator/=(const R &rhs) {
        for(auto &e: this) e /= rhs;
        return *this;
    }

    FPS& operator%=(const FPS &rhs) {
        return *this = div(rhs).second;
    }

    // TODO: >>=, <<= shifts
    // FPS& operator<<=(const int s) {
    //     return *this = (*this).shift(s);
    // }
    // FPS& operator>>=(const int s) {
    //      return operator<<=(-s);
    // }

    FPS operator-() const {
        FPS t = FPS(this->size());
        for(int i = 0; i < this->size(); i++) t
        for(auto &e: t.a) e = -e;
        return t;
    }

    FPS operator+(const FPS &rhs) const { return FPS(*this) += rhs; }
    FPS operator-(const FPS &rhs) const { return FPS(*this) -= rhs; }
    FPS operator*(const FPS &rhs) const { return FPS(*this) *= rhs; }
    FPS operator*(const R &rhs)   const { return FPS(*this) *= rhs; }
    FPS operator/(const FPS &rhs) const { return FPS(*this) /= rhs; }
    FPS operator/(const R &rhs)   const { return FPS(*this) /= rhs; }
    FPS operator%(const FPS &rhs) const { return FPS(*this) %= rhs; }

    R coef(size_t index) { return this->get(index); }

    // bool operator == (const FPS &rhs) const { normalize(); return *this == rhs; }
    // bool operator != (const FPS &rhs) const { normalize(); return *this != rhs; }

    /// @return FPS A copy of this FPS mod x^k (getting the first k coefficients)
    FPS mod_xk(int k) const {
        return FPS(std::vector<R>(this->begin(), this->begin() + std::min(k, (int)this->size())));
    }
    inline static FPS mod_xk(FPS &f, int k) const { return f.mod_xk(k); }

    /// @return FPS A copy of this FPS times x^k (coefficients shifted k steps to the right)
    FPS times_xk(int k) const {
        std::vector<R> coefs = *this;
        coefs.insert(coefs.begin(), k, R(0));
        return FPS(coefs);
    }
    inline static FPS times_xk(FPS &f, int k) const { return f.times_xk(k); }

    /// @return FPS A copy of this FPS with coefficients shifted k steps to the left
    FPS floordiv_xk(int k) const {
        return FPS(std::vector<R>(this->begin() + k, this->end()));
    }
    inline static FPS floordiv_xk(FPS &f, int k) const { return f.floordiv_xk(k); }

    FPS shift(int k) const {
        if(k > 0)
            return times_xk(k);
        if(k < 0)
            return floordiv_xk(-k);
        return FPS(*this);
    }
    inline static FPS shift(FPS &f, int k) const { return f.shift(k); }


    /// @return FPS P(x^{1/2})
    FPS sqrt_x() const {
        std::vector<R> ret(this->size(), R(0));
        for(int i = 0; i < this->size(); i++)
            ret[i] = get(2*i);
        return FPS(ret);
    }
    inline static FPS sqrt_x(FPS &f) const { return f.sqrt_x(); }

    /// @return FPS P(x^2)
    FPS sqr_x() const {
        std::vector<R> ret(2 * this->size(), R(0));
        for(int i = 0; i < this->size(); i++)
            ret[2*i] = this->a[i];
        return FPS(ret);
    }
    inline static FPS sqr_x(FPS &f) const { return f.sqr_x(); }

    /// @return FPS The multiplicative inverse of this FPS mod x^m
    FPS inv(int m = deg() + 1) const {
        assert(get(0) != R(0));
        if(m == 0) return FPS(0);
        FPS Q(R(1) / get(0));
        for(int i = 1; i < m; i *= 2)
            Q = (Q * R(2) - Q * Q * this->mod_xk(2*i)).mod_xk(2*i);
        return Q.mod_xk(m);
    }
    inline static FPS inv(FPS &f, int m = deg() + 1) const { return f.inv(m); }

    /// @return FPS  A copy of this FPS with coefficients reversed
    FPS reverse() const {
        return FPS(std::vector<R>(this->rbegin(), this->rend()));
    }
    inline static FPS reverse(FPS &f) const { return f.reverse(Q); }

    /// @brief Euclidean division of this FPS with the FPS P
    /// @return pair<FPS, FPS> the quotiend (first) and remainder (second) FPSs
    std::pair<FPS, FPS> div(const FPS &P) const {
        assert(!P.is_zero());
        if(P.deg() > this->deg())
            return {FPS(0), FPS(*this)};
        FPS A_R = this->reverse(), B_R = P.reverse();
        int d = this->deg() - P.deg();
        FPS D = ((A_R).mod_xk(d + 1) * P.inv(d + 1)).mod_xk(d + 1).reverse();
        return {D, *this - D * P};
    }
    inline static std::pair<FPS, FPS> div(FPS &f, const FPS &Q) const { return f.div(Q); }

    /// @return FPS A copy of the FPS with the x value scaled by v
    FPS scale_x(const R &v) const {
        R c = 1;
        FPS P(*this);
        for(int i = 0; i <= P.deg(); i++, c *= v) { 
            P[i] *= c;
        }
        return  P;
    }
    inline static FPS scale_x(FPS &f, const R &v) const { return f.scale_x(v); }



    /// @return The derivative of the FPS
    FPS deriv(int k = 1) const {
        if(deg() + 1 < k) return FPS(R(0));
        std::vector<R> res(deg() - k +  1, R(0));
        if(k == 1) {
            for(int i = 1; i <= deg(); i++)
                res[i - 1] = this->at(i) * R(i);
        } else {
            for(int i = k; i <= deg(); i++)
                res[i - k] = this->at(i) * internal::fact<R>(i) * internal::inv_fact<R>(i - k);
        }
        return FPS(res);
    }

    inline static FPS deriv(FPS &f, int k = 1) const { return f.deriv(k); }


    /// @return The integral of the FPS (with constant 0)
    FPS integr() const {
        std::vector<R> res(deg() + 2);
        for(int i = 1; i <= deg() + 1; i++) {
            res[i] = this->at(i-1) / i;
        }
        return FPS(res);
    }
    inline static FPS integr(FPS &f) const { return f.integr(k); }

    /// @return The (natural) logarithm of the FPS (mod x^m)
    FPS log(int m) const { return (FPS(*this).deriv().mod_xk(m) * FPS(*this).inv(m)).mod_xk(m);  }
    inline static FPS log(FPS &f, int m) const { return f.log(m); }

    FPS exp(int m) const {
        assert(get(0) == 0);
        FPS Q(1), P(*this);
        P[0] += 1;
        for(int i = 1; i < m; i *= 2) { // at the end of each step: F = e^P (mod x^{2^{i+1}})
            // Q = (Q * (P - Q.log(2*i))).mod_xk(2*i);
            Q = (Q * (P - Q.mod_xk(2*i))).mod_xk(2*i);
        }
        return Q.mod_xk(m);
    }
    inline static FPS exp(FPS &f, int m) const { return f.exp(m); }

    // FPS sqrt(int k = 2) const {

    // }
    // inline static FPS sqrt(FPS &f, int m) const { return f.sqrt(m); }

    /// @return The k-th power of the FPS
    FPS pow(int k, int m = -1) const {
        if(m == -1) m = deg() * k;
        if(is_zero()) {
            assert(k > 0);
            return (k == 0 ? FPS(1) : FPS(0));
        }

        if(k < 0)
            return pow(-k, m).inv();

        int ord = this->ord();
        R a = this->at(ord);
        FPS T(this->begin() + ord, this->end());
        T /= a; // P(x) = ax^ord T(x); T(0) = 1
        
        // P(x)^k = a^k x^{k*ord} T(x)^k = a^k x^{k*ord} exp[klnT(x)]
        return (T.log(m) * k).exp(m).shift(k * ord) * zeno::pow(a, k);
    }

    /// @brief Evaluate the FPS at point x0
    /// @return The value of the FPS at point x0 
    R eval(R x0) const {
        R y(0);
        for(int i = this->deg(); i >= 0; i--) {
            y = (x0 * y) + this->at(i);
        }
        return y;
    }

    /// @brief Builds the FPS evaluation tree.
    /// @return The FPS P(x) = (x - x_L) (x - x_{L+1}) ... (x - x_R)
    static FPS<R> build_poly_tree(std::vector<FPS<R>> &tree, 
        int v, std::vector<R>::iterator l, std::vector<R>::iterator r) {
            if(r - l == 1) {
                return tree[v] = FPS({-*l, 1});
            } else {
                auto m = l + (r - l) / 2;
                return tree[v] = build_poly_tree(tree, 2*v, l, m) 
                    * build_poly_tree(tree, 2*v+1, m, r);
            }
        }

    std::vector<R> _eval(std::vector<FPS<R>> &tree, 
        int v, std::vector<R>::iterator l, std::vector<R>::iterator r) {
            if(r - l == 1) {
                return {eval(*l)};
            } else {
                auto m = l + (r - l) / 2;
                std::vector<R> A = (*this % tree[2*v]).eval(tree, 2*v, l, m);
                std::vector<R> B = (*this % tree[2*v+1]).eval(tree, 2*v+1, m, r);
                A.insert(end(A), begin(B), end(B));
                return A;
            }
        }

    /// @brief Evaluate the FPS in points x0, x1, ..., x_{n-1}
    /// @return The value of the FPS at points x0, x1, ..., x_{n-1}
    std::vector<R> eval(std::vector<R> const &x) const {
        int n = x.size();
        if(is_zero()) return std::vector<R>(n, R(0));
        std::vector<FPS<R>> tree(4*n);
        build_poly_tree(tree, 1, begin(x), end(x));
        return _eval(tree, 1, begin(x), end(x));
    }


    static FPS<R> _interpolate_tree(std::vector<FPS<R>> const &tree, 
        int v, typename std::vector<R>::iterator l, typename std::vector<R>::iterator r) {
            if(r - l == 1) {
                return tree[v] = FPS<R>(*l);
            } else {
                auto m = l + (r - l) / 2;
                FPS<R> A = _interpolate_tree(tree, 2*v, l, m);
                FPS<R> B = _interpolate_tree(tree, 2*v+1, m, r);
                return tree[v] = A * tree[2*v+1] + tree[2*v] * B;
            }
    }


    static FPS<R> interpolate(std::vector<R> const &x, std::vector<R> const &y) {
        assert(x.size() == y.size());
        std::vector<FPS<R>> tree(4*x.size());
        std::vector<R> u = FPS<R>::build_poly_tree(tree, 1, begin(x), end(x)).deriv().eval(x);
        for(int i = 0; i < y.size(); i++)
            u[i] = y[i] / u[i];
        return _interpolate_tree(tree, 1, begin(u), end(u));
    }

    /// @return The FPS P(x) = Prod_{i=0}^n f_i(x)
    static FPS<R> eval_product(std::vector<R> &f) {
        std::deque<FPS<R>> q;
        for(auto &e: x) q.emplace_back(FPS(e));
        for(int i = 0; i < x.size() - 1; i++) {
            FPS<R> A = q.front();
            q.pop_front();
            FPS<R> B = q.front();
            q.pop_front();
            q.emplace_back(A * B);
        }
        return q.front();
    }


    /// @return The Borel transformation of the FPS (a_k x^k -> a_k / k! x^k)
    FPS<R> borel() const {
        FPS<R> ret = FPS(*this);
        for(int k = 0; k <= deg(); k++)
            ret[k] *= internal::inv_fact(k);
        return ret;
    }

    /// @return The inverse Borel (Laplace) transformation of the FPS (a_k x^k -> k!a_k x^k)
    FPS<R> inv_borel() const {
        FPS<R> ret = FPS(*this);
        for(int k = 0; k <= deg(); k++)
            ret[k] *= internal::fact(k);
        return ret;
    }


    static std::pair<FPS, FPS> euclid_division(FPS const &A, FPS const &B) {
        return A.div(B);
    }

    static void euclid_division(FPS const &A, FPS const &B, FPS &Qf, FPS &Rf) {
        auto [Qf, Rf] = A.div(B);
    }

    static FPS gcd(FPS const &A, FPS const &B) { // general gcd also works
        return (B.is_zero() ? A : FPS::gcd(B, A % B));
    }

    static FPS extended_gcd(FPS const &A, FPS const &B, FPS &U, FPS &V) {
        FPS D(A), V1 = FPS(0), V3 = FPS(B), T, Q, R;
        U = FPS(1);
        while(!V3.is_zero()) {
            FPS::euclid_division(D, V3, Q, R);
            T = U - V1 * Q;
            U = V1, D = V3;
            V1 = T, V3 = R;
        }
        V = (D - A*U) / B;
    }


    /// @return P(x + c) where P is this FPS.
    static FPS taylor_shift(R c) const {
        size_t n = this->deg();
        FPS a, b; a.resize(n), b.resize(n);
        R cc = c;
        for(int i = 0; i < n; i++) a[i] = internal::fact(i) * this->at(i);
        for(int i = 0; i < n; i++) b[i] = cc * internal::inv_fact(i), cc *= c;
        return (a * b).borel();
    }

    static FPS hadamard_product(FPS const &f) {
        FPS ret({});
        for(int i = 0; i < this->size() && i < f.size(); i++)
            ret.push_back(this->at(i) * f[i]);
        ret.normalize();
        return ret;
    }

};

template<class R>
using FPS = FormalPowerSeries<R>;

template<class R>
using Poly = FormalPowerSeries<R>;

template<class R>
FormalPowerSeries<R>  operator*(const R &lhs, const FormalPowerSeries<R> &rhs) {
    return FPS<R>(rhs) *= lhs;
}

template<class R>
FormalPowerSeries<R> operator/(const R &lhs, const FormalPowerSeries<R> &rhs) {
    return FPS(lhs) /= rhs;
}


template<class R>
std::ostream& operator<<(std::ostream& os, const FormalPowerSeries<R> &P) {
    os << "[ ";
    for(int i = 0; i <= P.deg(); i++)
        os << P.get(i) << " ";
    os << "] ";
    return os;
}

} // namespace zeno

#endif /* End of ZENO_FPS_HPP */