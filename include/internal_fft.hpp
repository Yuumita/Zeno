 #ifndef INTERNAL_ZENO_FFT_HPP
 #define INTERNAL_ZENO_FFT_HPP

#include <math.h>

namespace zeno {

namespace fft {

    size_t compute_convolution_size(size_t n, size_t m) {
        return (size_t(1) << internal::ceil_log2(n + m - 1));
    }

    /// @brief Complex numbers struct
    template<typename T>
    struct complex { 
        T x, y;

        T Re() { return x; }
        T Im() { return y; }

        complex(): x(0), y(0) {}
        complex(T x0): x(x0), y(0) {}
        complex(T x0, T y0): x(x0), y(y0) {}

        static complex polar(T rho, T angle) {
            return complex(rho * cos(angle), rho * sin(angle));
        }

        complex conj() const {
            return complex(x, -y);
        }

        T norm() const {
            return x*x + y*y;
        };

        complex& operator+=(const complex &rhs) {
            this->x += rhs.x; 
            this->y += rhs.y;
            return *this;
        }
        complex& operator-=(const complex &rhs) {
            this->x -= rhs.x; 
            this->y -= rhs.y;
            return *this;
        }
        complex& operator*=(const complex &rhs) {
            return *this = complex(
                this->x * rhs.x - this->y * rhs.y,
                this->x * rhs.y + this->y * rhs.x
            ); 
        }
        complex& operator/=(const complex &rhs) {
            return *this = complex(
                (x * rhs.x + y * rhs.y) / rhs.norm(),
                (y * rhs.x - x * rhs.y) / rhs.norm() 
            ); 
        }

        complex operator-() const {
            return complex(-x, -y);
        }
        complex operator+(const complex &rhs) const {
            return complex(*this) += rhs;
        }
        complex operator-(const complex &rhs) const {
            return complex(*this) -= rhs;
        }
        complex operator*(const complex &rhs) const {
            return complex(*this) *= rhs;
        }
        complex operator/(const complex &rhs) const {
            return complex(*this) /= rhs;
        }

        friend std::ostream &operator<<(std::ostream &os, const complex &z) { return os << z.x << " + " << z.y << "i"; }
    };


}  // namespace fft

} // namespace zeno

 #endif // INTERNAL_ZENO_FFT_HPP