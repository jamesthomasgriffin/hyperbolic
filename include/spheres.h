#pragma once

#include <assert.h>

#include <glm/glm.hpp>

#include "hyperbolic.h"

namespace hyperbolic {

namespace spheres {

namespace contains_point {

template <typename T, qualifier Q>
bool sphere(glm::vec<4, T, Q> const &center, hyperbolic_angle<T> const &radius,
            glm::vec<4, T, Q> const &p) {
  T cp = lorentz::dot(center, p);
  return -cp < cosh(radius);
}

template <typename T, qualifier Q>
bool horosphere(glm::vec<4, T, Q> const &center,
                hyperbolic_angle<T> const &radius, glm::vec<4, T, Q> const &p) {
  T cp = lorentz::dot(center, p);
  return -cp < exp(radius);
}

template <typename T, qualifier Q>
bool pseudosphere(glm::vec<4, T, Q> const &center,
                  hyperbolic_angle<T> const &radius,
                  glm::vec<4, T, Q> const &p) {
  T cp = lorentz::dot(center, p);
  return -cp < hyperbolic::sinh(radius);
}
} // namespace contains_point

namespace distance_from {

template <typename T, qualifier Q>
hyperbolic_angle<T> sphere(glm::vec<4, T, Q> const &center,
                           hyperbolic_angle<T> const &radius,
                           glm::vec<4, T, Q> const &p) {
  return hyperbolic_angle<T>::between_points(p, center) - radius;
}

template <typename T, qualifier Q>
hyperbolic_angle<T> horosphere(glm::vec<4, T, Q> const &center,
                               hyperbolic_angle<T> const &radius,
                               glm::vec<4, T, Q> const &p) {
  T exp_a = -lorentz::dot(center, p);
  T exp_minus_a = 1 / exp_a;
  hyperbolic_angle<T> d = hyperbolic_angle<T>::unsafe_bypass_checks(
      (exp_a + exp_minus_a) / 2, (exp_a - exp_minus_a) / 2);
  return d - radius;
}

template <typename T, qualifier Q>
hyperbolic_angle<T> pseudosphere(glm::vec<4, T, Q> const &center,
                                 hyperbolic_angle<T> const &radius,
                                 glm::vec<4, T, Q> const &p) {
  return -hyperbolic_angle<T>::between_point_and_plane(p, center) - radius;
}
} // namespace distance_from

namespace distance_between {

template <typename T, qualifier Q>
hyperbolic_angle<T>
spheres(glm::vec<4, T, Q> const &c1, hyperbolic_angle<T> const &r1,
        glm::vec<4, T, Q> const &c2, hyperbolic_angle<T> const &r2) {
  return hyperbolic_angle<T>::between_points(c1, c2) - r1 - r2;
}

template <typename T, qualifier Q>
hyperbolic_angle<T>
horospheres(glm::vec<4, T, Q> const &c1, hyperbolic_angle<T> const &r1,
            glm::vec<4, T, Q> const &c2, hyperbolic_angle<T> const &r2);

template <typename T, qualifier Q>
hyperbolic_angle<T>
pseudospheres(glm::vec<4, T, Q> const &c1, hyperbolic_angle<T> const &r1,
              glm::vec<4, T, Q> const &c2, hyperbolic_angle<T> const &r2);

template <typename T, qualifier Q>
hyperbolic_angle<T> sphere_and_horosphere(glm::vec<4, T, Q> const &c1,
                                          hyperbolic_angle<T> const &r1,
                                          glm::vec<4, T, Q> const &c2,
                                          hyperbolic_angle<T> const &r2) {
  return distance_from::horosphere(c2, r2, c1) - r1;
}

template <typename T, qualifier Q>
hyperbolic_angle<T> sphere_and_pseudosphere(glm::vec<4, T, Q> const &c1,
                                            hyperbolic_angle<T> const &r1,
                                            glm::vec<4, T, Q> const &c2,
                                            hyperbolic_angle<T> const &r2) {
  return -hyperbolic_angle<T>::between_point_and_plane(c1, c2) - (r1 + r2);
}

template <typename T, qualifier Q>
hyperbolic_angle<T> horosphere_and_pseudosphere(glm::vec<4, T, Q> const &c1,
                                                hyperbolic_angle<T> const &r1,
                                                glm::vec<4, T, Q> const &c2,
                                                hyperbolic_angle<T> const &r2);

} // namespace distance_between

namespace time_until_intersection_between {

template <typename T, qualifier Q>
hyperbolic_angle<T>
spheres(glm::vec<4, T, Q> const &c1, hyperbolic_angle<T> const &r1,
        glm::vec<4, T, Q> const &c2, hyperbolic_angle<T> const &r2,
        glm::vec<4, T, Q> const &v);

template <typename T, qualifier Q>
hyperbolic_angle<T> sphere_and_pseudosphere(glm::vec<4, T, Q> const &c1,
                                            hyperbolic_angle<T> const &r1,
                                            glm::vec<4, T, Q> const &c2,
                                            hyperbolic_angle<T> const &r2,
                                            glm::vec<4, T, Q> const &v);

template <typename T, qualifier Q>
hyperbolic_angle<T> sphere_and_horosphere(glm::vec<4, T, Q> const &c1,
                                          hyperbolic_angle<T> const &r1,
                                          glm::vec<4, T, Q> const &c2,
                                          hyperbolic_angle<T> const &r2,
                                          glm::vec<4, T, Q> const &v);

namespace fast {

template <typename T, qualifier Q>
inline T spheres(hyperbolic_angle<T> const &sum_radii, glm::vec<4, T, Q> const &c2,
          glm::vec<3, T, Q> const &v) {
  T denom = glm::dot(glm::vec<3, T, Q>{c2}, v);
  return (denom > 0) ? -(cosh(sum_radii) - c2[3]) / denom : -1;
}

template <typename T, qualifier Q>
inline T sphere_and_pseudosphere(hyperbolic_angle<T> const &sum_radii,
                          glm::vec<4, T, Q> const &c2,
                          glm::vec<3, T, Q> const &v) {
  T denom = glm::dot(glm::vec<3, T, Q>{c2}, v);
  return (denom > 0) ? -(sinh(sum_radii) - c2[3]) / denom : -1;
}

template <typename T, qualifier Q>
inline T sphere_and_horosphere(hyperbolic_angle<T> const &sum_radii,
                        glm::vec<4, T, Q> const &c2,
                        glm::vec<3, T, Q> const &v) {
  T denom = glm::dot(glm::vec<3, T, Q>{c2}, v);
  return (denom > 0) ? -(exp(sum_radii) - c2[3]) / denom : -1;
}

} // namespace fast

} // namespace time_until_intersection_between

namespace closest_point_to {

template <typename T, qualifier Q>
inline glm::vec<4, T, Q> sphere(glm::vec<4, T, Q> const &point,
                                glm::vec<4, T, Q> const &center,
                                hyperbolic_angle<T> const &radius) {
  return move_d_from_p_to_q(radius, center, point);
}

template <typename T, qualifier Q>
inline glm::vec<4, T, Q> pseudosphere(glm::vec<4, T, Q> const& point,
  glm::vec<4, T, Q> const& plane,
  hyperbolic_angle<T> const& radius) {

  T const a = lorentz::dot(point, plane);
  glm::vec<4, T, Q> const n = -(a * point + plane) / glm::sqrt(1 + a * a);
  hyperbolic_angle<T> dist = distance_from::pseudosphere(plane, radius, point);
  return move_d_from_p_along_n(dist, point, -n);
}

template <typename T, qualifier Q>
inline glm::vec<4, T, Q> horosphere(glm::vec<4, T, Q> const &point,
                             glm::vec<4, T, Q> const &point_at_inf,
                             hyperbolic_angle<T> const &radius) {

  glm::vec<4, T, Q> const n =
      point + (1 / lorentz::dot(point, point_at_inf)) * point_at_inf;
  hyperbolic_angle<T> dist =
      distance_from::horosphere(point_at_inf, radius, point);
  return move_d_from_p_along_n(dist, point, -n);
}

} // namespace closest_point_to

} // namespace spheres

} // namespace hyperbolic