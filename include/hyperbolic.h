#pragma once

#include <complex>

#include "glm/glm.hpp"
#include "glm/ext/matrix_clip_space.hpp"
#include "glm/ext/matrix_transform.hpp"

#include "vector_types.h"
#include "glm_lorentz.h"

#include "imaginary_extension.h"


namespace hyperbolic {

namespace lorentz = glm::lorentz;

inline scalar_t distance(vec4 const &p, vec4 const &q);

inline vec4 move_d_from_p_to_q(scalar_t d, vec4 const &p, vec4 const &q);

namespace detail {

inline void assert_close(scalar_t a, scalar_t b,
                         scalar_t epsilon = static_cast<scalar_t>(1e-6)) {
  assert(abs(a - b) < epsilon);
}

inline void assert_unit_tangent(vec4 const &v) {
  assert_close(lorentz::dot(v, v), 1);
}

inline void assert_positive(vec4 const &v) {
  assert(lorentz::dot(v, v) > 0);
}

inline void assert_negative(vec4 const &v) {
  assert(lorentz::dot(v, v) < 0);
}

inline void assert_point(vec4 const &v) {
  assert_close(lorentz::dot(v, v), -1);
}

inline void assert_null(vec4 const &v) {
  assert_close(lorentz::dot(v, v), 0);
}

// A class that represents an angle as a pair of cos/sin or cosh/sinh values.
template <int SIGN> class angle {
  scalar_t c, s;
  static_assert(SIGN * SIGN == 1, "SIGN must be +/-1");

public:
  angle() : c{1}, s{0} {}
  explicit angle(scalar_t theta)
      : c{(SIGN == 1) ? std::cosh(theta) : std::cos(theta)}, s{(SIGN == 1) ? std::sinh(theta)
                                                                 : std::sin(theta)} {
  }

  scalar_t get_angle() const { return (SIGN == 1) ? asinh(s) : atan2(s, c); }
  scalar_t get_c() const { return c; }
  scalar_t get_s() const { return s; }

  operator scalar_t() const { return get_angle(); };

  angle operator+(angle const &a) const {
    return angle{c * a.c + SIGN * s * a.s, c * a.s + s * a.c};
  }

  angle& operator+=(angle const &a) {
    scalar_t old_c = c;
    c = c * a.c + SIGN * s * a.s;
    s = old_c * a.s + s * a.c;
    return *this;
  }

  angle operator-(angle const &a) const {
    return angle{c * a.c - SIGN * s * a.s, -c * a.s + s * a.c};
  }

  angle &operator-=(angle const &a) {
    scalar_t old_c = c;
    c = c * a.c - SIGN * s * a.s;
    s = -old_c * a.s + s * a.c;
    return *this;
  }

  angle operator-() const { return angle{c, -s}; }

  static angle unsafe_bypass_checks(scalar_t _c, scalar_t _s) { return angle{_c, _s}; }
  static angle between_points(vec4 const& p, vec4 const& q) {
    scalar_t c = (SIGN == 1) ? -lorentz::dot(p, q) : glm::dot(p, q);
    scalar_t s = glm::sqrt(SIGN * (c * c - 1));
    return angle{c, s};
  }

private:
  angle(scalar_t _c, scalar_t _s) : c{_c}, s{_s} {};
};


} // namespace detail

using spherical_angle = detail::angle<-1>;
using hyperbolic_angle = detail::angle<1>;

inline scalar_t sin(spherical_angle const &angle) { return angle.get_s(); }
inline scalar_t cos(spherical_angle const &angle) { return angle.get_c(); }
inline scalar_t sinh(hyperbolic_angle const &angle) { return angle.get_s(); }
inline scalar_t cosh(hyperbolic_angle const &angle) { return angle.get_c(); }

inline scalar_t distance(vec4 const &p, vec4 const &q) {
  return std::acosh(-lorentz::dot(p, q));
}

inline vec4 move_d_from_p_to_q(hyperbolic_angle d, vec4 const& p, vec4 const& q) {
  hyperbolic_angle const pq = hyperbolic_angle::between_points(p, q);
  return (sinh(pq - d) * p + sinh(d) * q) / sinh(pq);
}

inline vec4 reflect(vec4 const &plane, vec4 const &v) {
  detail::assert_unit_tangent(plane);
  return v - 2 * lorentz::dot(plane, v) * plane;
}

inline vec4 midpoint(vec4 const &p, vec4 const &q) {
  return lorentz::normalize(p + q);
}

inline vec4 transport_unit_tangent_to(vec4 const &u, vec4 const &p) {
  detail::assert_point(p);
  detail::assert_unit_tangent(u);

  scalar_t pu = lorentz::dot(p, u);
  return glm::inversesqrt(1 + pu * pu) * (u + pu * p);
}

inline mat4 transport_frame_to(mat4 const &frame, vec4 const &p) {
  detail::assert_point(p);

  mat4 result{};
  result[0] = transport_unit_tangent_to(frame[0], p);
  result[1] = transport_unit_tangent_to(frame[1], p);
  result[2] = transport_unit_tangent_to(frame[2], p);
  result[3] = p;
  return result;
}

inline vec4 transport_tangent_to(vec4 const &v, vec4 const &p) {
  detail::assert_point(p);
  detail::assert_positive(v);

  scalar_t vv = lorentz::dot(v, v);
  scalar_t pv = lorentz::dot(p, v);
  return sqrt(vv / (vv + pv * pv)) * (v + pv * p);
}

inline vec4 closest_point_to_origin_in_plane(vec4 const &plane) {
  detail::assert_unit_tangent(plane);

  vec4 const v = vec4{0, 0, 0, 1} - plane[3] * plane;
  return glm::inversesqrt(1 + plane[3] * plane[3]) * v;
}

inline vec4 closest_point_in_plane(vec4 const &p, vec4 const &plane) {
  detail::assert_point(p);
  detail::assert_unit_tangent(plane);

  scalar_t const a = lorentz::dot(p, plane);
  vec4 v = p - a * plane;
  return glm::inversesqrt(1 + a * a) * v;
}

inline mat4 horographic(scalar_t fov, scalar_t aspect) {
  return glm::perspective(fov, aspect, 0.00001f, 2.0f) *
         glm::lookAt(glm::vec3{0, 0, 1}, glm::vec3{0, 0, -1},
                     glm::vec3{0, 1, 0});
}

} // namespace hyperbolic