#pragma once

#include <stdexcept>

template<typename T>
struct is_imaginary_extension {};

template <typename R> struct ImaginaryExtension {

  R real{};
  R imag{};

  constexpr ImaginaryExtension(R r, R i) : real{r}, imag{i} {};
  constexpr ImaginaryExtension(R r) : real{r}, imag{0} {};

  constexpr R norm() const { return real * real + imag * imag; }
  constexpr ImaginaryExtension conjugate() const {
    return ImaginaryExtension{real, 0 - imag};
  }

  constexpr bool operator==(ImaginaryExtension const &b) const {
    return (real == b.real) && (imag == b.imag);
  }
  constexpr bool operator!=(ImaginaryExtension const &b) const {
    return (real != b.real) || (imag != b.imag);
  }

  constexpr ImaginaryExtension operator+(ImaginaryExtension const &b) const {
    return ImaginaryExtension{real + b.real, imag + b.imag};
  }
  constexpr ImaginaryExtension operator-(ImaginaryExtension const &b) const {
    return ImaginaryExtension{real - b.real, imag - b.imag};
  }
  ImaginaryExtension &operator+=(ImaginaryExtension const &b) {
    real += b.real;
    imag += b.imag;
    return *this;
  }
  ImaginaryExtension &operator-=(ImaginaryExtension const &b) {
    real -= b.real;
    imag -= b.imag;
    return *this;
  }
  template<typename S>
  constexpr auto operator*(ImaginaryExtension<S> const &b) const -> ImaginaryExtension<decltype(real * b.real)> {
    return ImaginaryExtension<decltype(real * b.real)>{
        real * b.real - imag * b.imag, real * b.imag + imag * b.real};
  }
  ImaginaryExtension &operator*=(ImaginaryExtension const &b) {
    auto oldr = real;
    real = oldr * b.real - imag * b.imag;
    imag = oldr * b.imag + imag * b.real;
    return *this;
  }
  constexpr ImaginaryExtension operator-() const {
    return ImaginaryExtension{0 - real, 0 - imag};  // Stops MSVC complaining about unary minus for unsigned types
  }

  template <typename S> constexpr operator ImaginaryExtension<S>() const {
    return ImaginaryExtension<S>{static_cast<S>(real), static_cast<S>(imag)};
  }
};

template <typename T> struct is_imaginary_extension<ImaginaryExtension<T>> {
  static const bool value = true;
};

template <typename R> void swap(R &a, R &b) {
  auto c = a;
  a = b;
  b = a;
}

template <typename R, typename S, std::enable_if_t<!is_imaginary_extension<R>::value>>
auto operator*(R r, ImaginaryExtension<S> const &a) -> ImaginaryExtension<decltype(r * a.real)> {
  return ImaginaryExtension<R>{a.real * r, a.imag * r};
}

template <typename OS, typename R>
OS &operator<<(OS &ostr, ImaginaryExtension<R> const &a) {
  ostr << a.real << " + " << a.imag << "i";
  return ostr;
}

template <typename R> static inline R norm(ImaginaryExtension<R> const &a) {
  return a.norm();
}

template <typename Z> struct GaussianInteger : public ImaginaryExtension<Z> {

  constexpr GaussianInteger(ImaginaryExtension<Z> &&e)
      : ImaginaryExtension<Z>{std::move(e)} {}
  constexpr GaussianInteger(Z a) : ImaginaryExtension<Z>{a} {}
  constexpr GaussianInteger(Z a, Z b) : ImaginaryExtension<Z>{a, b} {}

  GaussianInteger operator/(Z d) const {

    if (d == 0)
      throw std::overflow_error{"Division by zero integer"};

    Z re = this->real;
    Z im = this->imag;

    if (d < 0) {
      re *= -1;
      im *= -1;
      d *= -1;
    }

    auto quot = [](Z a, Z b) { 
      int sgn = (a < 0) ? -1 : 1;
      a = (a < 0) ? -a : a;

      Z q = a / b, r = a % b;
      return sgn *((2 * r <= b) ? q : q + 1);
    };
    
    return GaussianInteger{quot(re, d), quot(im, d)};
  }

  GaussianInteger operator/(GaussianInteger const &d) const {
    auto m = d.norm();
    if (m == 0)
      throw std::overflow_error{"Division by zero Gaussian integer"};

    GaussianInteger c = *this * d.conjugate();

    return c / m;
  }

  GaussianInteger operator%(GaussianInteger const &d) const {
    return *this - (*this / d) * d;
  }
};

template <typename Re> struct ComplexNumber : public ImaginaryExtension<Re> {
  ComplexNumber(ImaginaryExtension<Re> &&e)
      : ImaginaryExtension<Re>{std::move(e)} {}
  ComplexNumber(Re r) : ImaginaryExtension<Re>{r} {}

  ComplexNumber operator/(ComplexNumber const &d) const {
    return ComplexNumber{(*this * d.conjugate()) / d.norm()};
  }
};

template <typename T> struct EuclidsAlgorithm {
  struct Result {
    T gcd;
    T m;
    T n;
  };
  static Result bezouts_identity(T a, T b) {

    T c{1}, d{0}, e{0}, f{1};
    if (norm(a) < norm(b)) {
      c = 0;
      d = 1;
      e = 1;
      f = 0;
      std::swap(a, b);
    }

    while (norm(b) > 0) {
      T q = a / b;
      a -= q * b;
      c -= q * e;
      d -= q * f;
      std::swap(a, b);
      std::swap(c, e);
      std::swap(d, f);
    }
    return Result{a, c, d};
  }
};
