#pragma once

#include <complex>

#include "glm/glm.hpp"
#include "glm/ext/matrix_clip_space.hpp"
#include "glm/ext/matrix_transform.hpp"

#include "glm_lorentz.h"
#include "lorentz_lie_algebra.h"

namespace hyperbolic {

namespace lorentz = glm::lorentz;

template <typename T, qualifier Q>
inline typename T distance(glm::vec<4, T, Q> const &p,
                           glm::vec<4, T, Q> const &q);

template <typename T, qualifier Q>
inline glm::vec<4, T, Q> move_d_from_p_to_q(T d, glm::vec<4, T, Q> const &p,
                                            glm::vec<4, T, Q> const &q);

template <typename T, qualifier Q>
inline glm::vec<4, T, Q> midpoint(glm::vec<4, T, Q> const &p,
                                  glm::vec<4, T, Q> const &q);

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q>
affine_transform_to_hyperbolic(glm::mat<4, 4, T, Q> const &p);

namespace detail {

template <typename T, typename U>
inline void assert_close(T a, U b, T epsilon = static_cast<T>(1e-6)) {
  assert(abs(a - b) < epsilon);
}

template <typename T, qualifier Q>
inline void assert_unit_tangent(glm::vec<4, T, Q> const &v) {
  assert_close(lorentz::dot(v, v), 1);
}

template <typename T, qualifier Q>
inline void assert_positive(glm::vec<4, T, Q> const &v) {
  assert(lorentz::dot(v, v) > 0);
}

template <typename T, qualifier Q>
inline void assert_negative(glm::vec<4, T, Q> const &v) {
  assert(lorentz::dot(v, v) < 0);
}

template <typename T, qualifier Q>
inline void assert_point(glm::vec<4, T, Q> const &v) {
  assert_close(lorentz::dot(v, v), -1);
}

template <typename T, qualifier Q>
inline void assert_null(glm::vec<4, T, Q> const &v) {
  assert_close(lorentz::dot(v, v), 0);
}

// A class that represents an angle as a pair of cos/sin or cosh/sinh values.
// This can avoid calls to inverse (hyperbolic) trig fns.
template <int SIGN, typename T> class angle {
  T c, s;
  static_assert(SIGN * SIGN == 1, "SIGN must be +/-1");

public:
  angle() : c{1}, s{0} {}
  explicit angle(T theta)
      : c{(SIGN == 1) ? std::cosh(theta) : std::cos(theta)},
        s{(SIGN == 1) ? std::sinh(theta) : std::sin(theta)} {}

  T get_angle() const { return (SIGN == 1) ? asinh(s) : atan2(s, c); }
  T get_c() const { return c; }
  T get_s() const { return s; }

  operator T() const { return get_angle(); };

  angle operator+(angle const &a) const {
    return angle{c * a.c + SIGN * s * a.s, c * a.s + s * a.c};
  }

  angle &operator+=(angle const &a) {
    T old_c = c;
    c = c * a.c + SIGN * s * a.s;
    s = old_c * a.s + s * a.c;
    return *this;
  }

  angle operator-(angle const &a) const {
    return angle{c * a.c - SIGN * s * a.s, -c * a.s + s * a.c};
  }

  angle &operator-=(angle const &a) {
    T old_c = c;
    c = c * a.c - SIGN * s * a.s;
    s = -old_c * a.s + s * a.c;
    return *this;
  }

  angle operator-() const { return angle{c, -s}; }

  static angle unsafe_bypass_checks(T _c, T _s) { return angle{_c, _s}; }
  template <length_t L, qualifier Q>
  static angle between_points(glm::vec<L, T, Q> const &p,
                              glm::vec<L, T, Q> const &q) {
    T c = (SIGN == 1) ? -lorentz::dot(p, q) : glm::dot(p, q);
    T s = glm::sqrt(SIGN * (c * c - 1));
    return angle{c, s};
  }

private:
  angle(T _c, T _s) : c{_c}, s{_s} {};
};

template <typename T, qualifier Q>
inline glm::mat<3, 3, T, Q>
log_special_orthogonal_matrix(glm::mat<3, 3, T, Q> const &m) {
  using vec3_t = glm::vec<3, T, Q>;
  using mat3_t = glm::mat<3, 3, T, Q>;

  // NB angle is +ve, sign is encoded in K
  mat3_t const sinc_angle_K = (m - glm::transpose(m)) * static_cast<T>(0.5);

  T const cos_angle = (glm::trace(m) - 1) / 2;
  T const angle = glm::acos(cos_angle); // In [0, pi]
  if (angle < static_cast<T>(1e-6))     // sinc is approx 1 for small angles
    return sinc_angle_K;

  T const sinc_angle = glm::sin(angle) / angle;
  if (sinc_angle > static_cast<T>(1e-6))
    return sinc_angle_K / sinc_angle;

  // Angle is close to pi, so find kernel of m - I, which is orthogonal to the
  // span (column space)
  mat3_t const n = m - mat3_t{1};

  // Both of these candidates cannot be zero
  vec3_t const cand1 = glm::cross(n[0], n[1]);
  vec3_t const cand2 = glm::cross(n[0], n[2]);
  vec3_t const cand3 = glm::cross(n[1], n[2]);
  T const cc1 = glm::dot(cand1, cand1);
  T const cc2 = glm::dot(cand2, cand2);
  T const cc3 = glm::dot(cand2, cand2);
  vec3_t const cand = (cc1 > static_cast<T>(0.01))   ? cand1
                      : (cc2 > static_cast<T>(0.01)) ? cand2
                                                     : cand3;
  vec3_t const w = static_cast<T>(3.1415926) * glm::normalize(cand);

  return mat3_t{0, w[2], -w[1], -w[2], 0, w[0], w[1], -w[0], 0};
}

// m is assumed to be of the form of a 3x3 special orthogonal matrix and a
// 4-vector with 4th component 1
template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q>
log_special_affine_matrix(glm::mat<4, 4, T, Q> const &m) {
  using vec3_t = glm::vec<3, T, Q>;
  using vec4_t = glm::vec<4, T, Q>;
  using mat3_t = glm::mat<3, 3, T, Q>;
  using mat4_t = glm::mat<4, 4, T, Q>;

  mat3_t S = log_special_orthogonal_matrix(mat3_t{m});

  vec3_t p = vec3_t{m[3]};

  // p = Bv = S^-1(exp(S) - 1)v, interpreted as a power series (S is not
  // invertible).  Using S^3 = -a^2 S, we evaluate B,
  // B = I + (cos(a) - 1) / a^2) S + (sin(a) - a) / a^3 S^2

  T const a2 = -0.5f * glm::trace(S * S);
  T const a = glm::sqrt(a2);
  T const cos_a = glm::cos(a);
  T const sin_a = glm::sin(a);
  T const coeff1 = (a < 1e-5) ? (a2 - 12) / 24 : (cos_a - 1) / a2;
  T const coeff2 = (a < 1e-4) ? (a2 - 20) / 120 : (sin_a - a) / (a * a2);
  mat3_t const B = mat3_t{1} + coeff1 * S + coeff2 * S * S;
  vec3_t const v = glm::inverse(B) * p;

  mat4_t logm{S};
  logm[3] = vec4_t{v, 0};
  return logm;
}

} // namespace detail

template <typename T> using spherical_angle = detail::angle<-1, T>;
template <typename T> using hyperbolic_angle = detail::angle<1, T>;

template <typename T> inline T sin(spherical_angle<T> const &angle) {
  return angle.get_s();
}
template <typename T> inline T cos(spherical_angle<T> const &angle) {
  return angle.get_c();
}
template <typename T> inline T sinh(hyperbolic_angle<T> const &angle) {
  return angle.get_s();
}
template <typename T> inline T cosh(hyperbolic_angle<T> const &angle) {
  return angle.get_c();
}

template <typename T, qualifier Q>
T distance(glm::vec<4, T, Q> const &p, glm::vec<4, T, Q> const &q) {
  return std::acosh(-lorentz::dot(p, q));
}

template <typename T, qualifier Q>
inline glm::vec<4, T, Q> move_d_from_p_to_q(hyperbolic_angle<T> d,
                                            glm::vec<4, T, Q> const &p,
                                            glm::vec<4, T, Q> const &q) {
  using hyp_angle = hyperbolic_angle<T>;
  hyp_angle const pq = hyp_angle::between_points(p, q);
  return (sinh(pq - d) * p + sinh(d) * q) / sinh(pq);
}

template <typename T, qualifier Q>
inline glm::vec<4, T, Q> reflect(glm::vec<4, T, Q> const &plane,
                                 glm::vec<4, T, Q> const &v) {
  detail::assert_unit_tangent(plane);
  return v - 2 * lorentz::dot(plane, v) * plane;
}

template <typename T, qualifier Q>
inline glm::vec<4, T, Q> midpoint(glm::vec<4, T, Q> const &p,
                                  glm::vec<4, T, Q> const &q) {
  return lorentz::normalize(p + q);
}

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q>
affine_transform_to_hyperbolic(glm::mat<4, 4, T, Q> const &m) {
  using vec3_t = glm::vec<3, T, Q>;
  using lie_algebra_t = lorentz_lie_algebra_t<vec3_t>;

  glm::mat<4, 4, T, Q> log_m = detail::log_special_affine_matrix(m);

  // NB from_matrix looks only at the final column of log_m to extract the boost
  // direction which is exactly what we want
  lie_algebra_t S = lie_algebra_t::from_matrix(log_m);
  return exponential(S);
}

template <typename T, qualifier Q>
inline glm::vec<4, T, Q> transport_unit_tangent_to(glm::vec<4, T, Q> const &u,
                                                   glm::vec<4, T, Q> const &p) {
  detail::assert_point(p);
  detail::assert_unit_tangent(u);

  T pu = lorentz::dot(p, u);
  return glm::inversesqrt(1 + pu * pu) * (u + pu * p);
}

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q>
transport_frame_to(glm::mat<4, 4, T, Q> const &frame,
                   glm::vec<4, T, Q> const &p) {
  detail::assert_point(p);

  glm::mat<4, 4, T, Q> result{};
  result[0] = transport_unit_tangent_to(frame[0], p);
  result[1] = transport_unit_tangent_to(frame[1], p);
  result[2] = transport_unit_tangent_to(frame[2], p);
  result[3] = p;
  return result;
}

template <typename T, qualifier Q>
inline glm::vec<4, T, Q> transport_tangent_to(glm::vec<4, T, Q> const &v,
                                              glm::vec<4, T, Q> const &p) {
  detail::assert_point(p);
  detail::assert_positive(v);

  T vv = lorentz::dot(v, v);
  T pv = lorentz::dot(p, v);
  return glm::sqrt(vv / (vv + pv * pv)) * (v + pv * p);
}

template <typename T, qualifier Q>
inline glm::vec<4, T, Q>
closest_point_to_origin_in_plane(glm::vec<4, T, Q> const &plane) {
  using vec4_t = glm::vec<4, T, Q>;
  detail::assert_unit_tangent(plane);

  vec4_t const v = vec4_t{0, 0, 0, 1} - plane[3] * plane;
  return glm::inversesqrt(1 + plane[3] * plane[3]) * v;
}

template <typename T, qualifier Q>
inline glm::vec<4, T, Q>
closest_point_in_plane(glm::vec<4, T, Q> const &p,
                       glm::vec<4, T, Q> const &plane) {
  detail::assert_point(p);
  detail::assert_unit_tangent(plane);

  T const a = lorentz::dot(p, plane);
  glm::vec<4, T, Q> v = p - a * plane;
  return glm::inversesqrt(1 + a * a) * v;
}

template <typename T>
inline glm::mat<4, 4, T, glm::defaultp> horographic(T fov, T aspect) {
  return glm::perspective(fov, aspect, static_cast<T>(0.00001), static_cast<T>(2.0)) *
         glm::lookAt(glm::vec<3, T>{0, 0, 1}, glm::vec<3, T>{0, 0, -1},
                     glm::vec<3, T>{0, 1, 0});
}

} // namespace hyperbolic