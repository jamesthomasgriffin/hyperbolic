#pragma once

#include <assert.h>

#include <glm/glm.hpp>

#include "hyperbolic.h"

namespace hyperbolic {

// Represents a sphere, horosphere, pseudo-sphere (NB a plane is also a pseudo-sphere)
template <typename T, qualifier Q> struct implicit_sphere {
  enum class Type { Sphere, Horosphere, Pseudosphere, Plane };

  glm::vec<4, T, Q> center; // Should have norm -1 for sphere, 0 for horosphere,
                            // 1 for pseudosphere and plane
  hyperbolic_angle<T> radius{};
  Type type;

  T curvature() const;

  bool contains(glm::vec<4, T, Q> const &p) const;
  implicit_sphere grow_by(hyperbolic_angle<T> const &r) const;

  T distance_from(glm::vec<4, T, Q> const &p) const;

  bool intersects(implicit_sphere const &s) const;
  T time_to_intersection(implicit_sphere const &s,
                         glm::vec<4, T, Q> const &v) const;

  glm::vec<4, T, Q> closest_point(glm::vec<4, T, Q> const &p) const;
  std::array<glm::vec<4, T, Q>, 2>
  closest_points(implicit_sphere const &sphere2) const;
};

template <typename T, qualifier Q>
inline implicit_sphere<T, Q> operator*(glm::mat<4, 4, T, Q> const &m,
                                    implicit_sphere<T, Q> const &s) {
  return implicit_sphere<T, Q>{s.type, m * s.center, s.radius};
}

template <typename T, qualifier Q> inline T implicit_sphere<T, Q>::curvature() const {

  if (type == Type::Sphere)
    return 1 / std::pow(hyperbolic::sinh(radius), 2);
  else if (type == Type::Horosphere)
    return 0;
  else // if (type == Type::Pseudosphere)
    return -1 / std::pow(hyperbolic::cosh(radius), 2);
}

template <typename T, qualifier Q>
inline T
implicit_sphere<T, Q>::distance_from(glm::vec<4, T, Q> const &p) const {

  T const cp = lorentz::dot(c, p);
  if (type == Type::Sphere)
    return glm::acosh(-cp) - static_cast<T>(radius);
  else if (type == Type::Horosphere)
    return glm::log(-cp);
  else // if (type == Type::Pseudosphere)
    return glm::asinh(cp) - static_cast<T>(radius);
}

template <typename T, qualifier Q>
inline bool implicit_sphere<T, Q>::contains(glm::vec<4, T, Q> const &p) const {
  T const cp = glm::lorentz::dot(center, p);

  if (type == Type::Horosphere)
    return cp < 1;
  else if (type == Type::Sphere)
    return -cp < hyperbolic::cosh(radius);
  else // if (type == Type::Pseudosphere)
    return cp < hyperbolic::sinh(radius);
}

template <typename T, qualifier Q>
inline bool implicit_sphere<T, Q>::intersects(implicit_sphere const &s) const {
  implicit_sphere *s1 = this;
  implicit_sphere *s2 = &s;
  if (static_cast<int>(s1->type) > static_cast<int>(s2->type))
    std::swap(s1, s2);

  if (s1->type == Type::Sphere) {
    if (s2->type == Type::Sphere) {
      hyperbolic_angle<T> const sum_radii = s1->radius + s2->radius;
      T const coshd = -lorentz::dot(s1->center, s2->center);
      return coshd < hyperbolic::cosh(sum_radii);
    } else if (s2->type == Type::Horosphere) {

    } else if (s2->type == Type::Sphere) {

    }
  } else if (s1->type == Type::Horosphere) {
    if (s2->type == Type::Horosphere) {

    } else if (s2->type == Type::Sphere) {

    }
  } else if (s1->type == Type::Pseudosphere) {
    if (s2->type == Type::Pseudosphere) {

    } 
  }
  assert(false);
  return false;
}

template <typename T, qualifier Q>
inline glm::vec<4, T, Q>
implicit_sphere<T, Q>::closest_point(glm::vec<4, T, Q> const &p) const {
  return move_d_from_p_to_q(radius, center, p);
}

template <typename T, qualifier Q>
inline std::array<glm::vec<4, T, Q>, 2>
implicit_sphere<T, Q>::closest_points(implicit_sphere const &sphere2) const {
  return {move_d_from_p_to_q(radius, center, sphere2.center),
          move_d_from_p_to_q(sphere2.radius, sphere2.center, center)};
}

template <typename T, qualifier Q>
inline implicit_sphere<T, Q>
implicit_sphere<T, Q>::grow_by(hyperbolic_angle<T> const &r) const {
  return implicit_sphere{center, radius + r, type};
}

} // namespace hyperbolic