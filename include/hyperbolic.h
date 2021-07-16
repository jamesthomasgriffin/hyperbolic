#pragma once

#include <complex>

#include "glm/glm.hpp"
#include "glm/ext/matrix_clip_space.hpp"
#include "glm/ext/matrix_transform.hpp"

#include "glm_lorentz.h"

#include "imaginary_extension.h"

namespace glm {
namespace hyperbolic {



using T = vec4::value_type;

inline T distance(vec4 const &p, vec4 const &q) {
  return std::acosh(lorentz::dot(p, q));
}

namespace detail {

// A class that represents an angle as a pair of cos/sin or cosh/sinh values.
template <int SIGN> class angle {
  T c, s;
  static_assert(SIGN * SIGN == 1, "SIGN must be +/-1");

public:
  angle() : c{1}, s{0} {}
  angle(T theta)
      : c{(SIGN == 1) ? cosh(theta) : cos(theta)}, s{(SIGN == 1) ? sinh(theta)
                                                                 : sin(theta)} {
  }

  T get_angle() const { return (SIGN == 1) ? acosh(c) : atan2(s, c); }
  T get_c() const { return c; }
  T get_s() const { return s; }

  operator T() const { return get_angle(); };

  angle operator+(angle const &a) const {
    return angle{c * a.c + SIGN * s * a.s, c * a.s + s * a.c};
  }

  angle& operator+=(angle const &a) const {
    T old_c = c;
    c = c * a.c + SIGN * s * a.s;
    s = old_c * a.s + s * a.c;
    return *this;
  }

  angle operator-(angle const &a) const {
    return angle{c * a.c + s * a.s, SIGN * c * a.s + s * a.c};
  }

  angle operator-() const { return angle{c, -s}; }

  static angle unsafe_bypass_checks(T _c, T _s) { return angle{_c, _s}; }

private:
  angle(T _c, T _s) : c{_c}, s{_s} {};
};

} // namespace detail
using spherical_angle = detail::angle<-1>;
using hyperbolic_angle = detail::angle<1>;

inline T sin(spherical_angle const &angle) { return angle.get_s(); }
inline T cos(spherical_angle const &angle) { return angle.get_c(); }
inline T sinh(hyperbolic_angle const &angle) { return angle.get_s(); }
inline T cosh(hyperbolic_angle const &angle) { return angle.get_c(); }

class point {
private:
  vec4 m_v;

public:
  point(vec4 const &p) : m_v{lorentz::normalize(p, -1)} {}

  vec4 const &vec() const { return m_v; }

  static point unsafe_bypass_checks(vec4 const &p) { return point{p, 1}; }

private:
  point(vec4 const &p, int) : m_v{p} {}
};

class boundary_point {
private:
  vec4 m_v;

public:
  boundary_point(vec4 const &p) : m_v{p} {
    assert(lorentz::dot(p, p) < 1e-12);
  }

  vec4 const &vec() const { return m_v; }

  static boundary_point unsafe_bypass_checks(vec4 const &p) {
    return boundary_point{p, 1};
  }

private:
  boundary_point(vec4 const &p, int) : m_v{p} {}
};

class unit_tangent {
private:
  vec4 m_v;

public:
  unit_tangent(vec4 const &n) : m_v{lorentz::normalize(n, 1)} {}

  vec4 const &get_normal() const { return m_v; }

  T operator()(point const &p) const {
    return lorentz::dot(get_normal(), p.vec());
  }
  T operator()(boundary_point const &p) const {
    return lorentz::dot(get_normal(), p.vec());
  }

  static unit_tangent unsafe_bypass_checks(vec4 const &n) { return unit_tangent{n, 1}; }

private:
  unit_tangent(vec4 const &n, int) : m_v{n} {}
};


inline point reflect(unit_tangent const &p, point const &v) {
  return point::unsafe_bypass_checks(v.vec() -
                                     2 * p(v) * p.get_normal());
}

inline unit_tangent reflect(unit_tangent const &p1, unit_tangent const &p2) {
  vec4 const &v = p1.get_normal();
  vec4 const &w = p2.get_normal();
  return unit_tangent::unsafe_bypass_checks(w - 2 * lorentz::dot(v, w) * v);
}

inline boundary_point reflect(unit_tangent const &p, boundary_point const &w) {
  return boundary_point::unsafe_bypass_checks(w.vec() -
                                              2 * p(w) * p.get_normal());
}

inline point midpoint(point const &p, point const &q) {
  return point{p.vec() + q.vec()};
}

inline vec4 transport_unit_tangent_to(vec4 const &u, vec4 const &p) {
  T pu = lorentz::dot(p, p);
  return inversesqrt(1 + pu * pu) * (u + pu * p);
}

inline mat4 transport_frame_to(mat4 const &frame, point const &p) {
  mat4 result{};
  result[0] = transport_unit_tangent_to(frame[0], p.vec());
  result[1] = transport_unit_tangent_to(frame[1], p.vec());
  result[2] = transport_unit_tangent_to(frame[2], p.vec());
  result[3] = p.vec();
  return result;
}

inline vec4 transport_tangent_to(vec4 const &v, vec4 const &p) {
  T vv = lorentz::dot(v, v);
  T pv = lorentz::dot(p, p);
  return sqrt(vv / (vv + pv * pv)) * (v + pv * p);
}

inline vec4 closest_point(unit_tangent const &pl) {
  vec4 const &n = pl.get_normal();
  vec4 const v = vec4{0, 0, 0, 1} - n[3] * n;
  return inversesqrt(1 + n[3] * n[3]) * v;
}

inline point closest_point(point const &p, unit_tangent const &pl) {
  T a = pl(p);
  vec4 v = p.vec() - a * pl.get_normal();
  return point::unsafe_bypass_checks(inversesqrt(1 + a * a) * v);
}

inline mat4 horographic(T fov, T aspect) {
  return glm::perspective(fov, aspect, 0.00001f, 2.0f) *
         glm::lookAt(glm::vec3{0, 0, 1}, glm::vec3{0, 0, -1},
                     glm::vec3{0, 1, 0});
}

template <typename T>
mat<4, 4, T> exponentiate_from_tangent(vec<3, T> const &e, vec<3, T> const &f) {

  std::complex<T> const detA{glm::dot(e, e) / 4 - glm::dot(f, f) / 4, glm::dot(e, f) / 2};

  std::complex<T> const a = std::sqrt(detA);
  
  // Below values are independent of chosen root above
  std::complex<T> const sinc_a = (std::norm(a) < 1e-12) ? a : sin(a) / a;
  std::complex<T> const cos_a = cos(a);

  ImaginaryExtension<mat<2, 2, T>> A{
      mat<2, 2, T>{f[0], f[1] - e[2], f[1] - e[2], -f[0]},
      mat<2, 2, T>{e[0], e[1] - f[2], e[1] - f[2], -f[0]}};

  ImaginaryExtension<mat<2, 2, T>> expA =
      ImaginaryExtension<mat<2, 2, T>>{mat<2, 2, T>(cos_a.real()),
                                       mat<2, 2, T>(cos_a.imag())} +
      ImaginaryExtension<T>{sinc_a.real(), sinc_a.imag()} * A;

  return lorentz::detail::action_via_hermitian_matrix(expA);

  // mat<4, 4, T> const Z{0,     -e[0], e[1], f[0], e[2], 0,    -e[0], f[1],
  //               -e[1], e[0],  0,    f[2], f[0], f[1], f[2],  0};
  //
  // mat<4, 4, T> const iZ{0,     -f[0], f[1], -e[0], f[2],  0,     -f[0],
  // -e[1],
  //                -f[1], f[0],  0,    -e[2], -e[0], -e[1], -e[2], 0};
  //
  // mat<4, 4, T> bZ = b.real() * Z - b.imag() * iZ;
  //
  // vec<4, T> result{0, 0, 0, cos_norm_a};
  // result += bZ[3];
  // result += (sinc_norm_a / 2) * vec<4, T>{cross(f, e), norm_a};
  //
  // T test_norm = lorentz::dot(result, result);

  // return result;
}

} // namespace hyperbolic

} // namespace glm