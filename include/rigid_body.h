#pragma once

#include <array>
#include <math.h>

//#include "vector_types.h"
#include <glm/glm.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include "lorentz_lie_algebra.h"

#include "hyperbolic.h"

namespace hyperbolic {

template<typename T, qualifier Q = glm::highp>
struct rigid_body_t {
  using velocity_t = lorentz_lie_algebra_t<T, Q>;
  using scalar_t = T;
  using mat4_t = glm::mat<4, 4, T, Q>;
  using vec3_t = glm::vec<3, T, Q>;
  using vec4_t = glm::vec<4, T, Q>;

  mat4_t frame{1}; // A Lorentz matrix
  velocity_t velocity{};
  vec3_t distributional_inertia{1.0f};
  scalar_t total_mass{1};

  scalar_t rotational_energy() const;
  scalar_t boost_energy() const;
  scalar_t kinetic_energy() const;
  velocity_t local_momentum() const;
  velocity_t momentum() const;
  velocity_t calculate_free_acceleration() const;
  void free_motion(T dt);
  velocity_t apply_inertia(velocity_t const &v) const;
  velocity_t apply_inverse_inertia(velocity_t const &p) const;

  void apply_impulse(vec4_t const &location, vec4_t const &direction);
  void apply_local_impulse(vec4_t const &location, vec4_t const &direction);
  static velocity_t calculate_local_impulse(vec4_t const &location,
                                            vec4_t const &direction);

  static rigid_body_t
  from_distribution_matrix(mat4_t const &M); // Returns a rigid body at rest
  static scalar_t
  impulse_magnitude_from_elastic_collision(rigid_body_t const &b1,
                                           rigid_body_t const &b2,
                                           vec4_t const &p, vec4_t const &n);
  static scalar_t impulse_magnitude_from_inelastic_collision(
      rigid_body_t const &b1, rigid_body_t const &b2, vec4_t const &p,
      vec4_t const &n, T energy_loss_factor);
};

template <typename T, qualifier Q>
inline T rigid_body_t<T, Q>::rotational_energy() const {
  velocity_t const p = local_momentum();
  return glm::dot(p.rotational, velocity.rotational);
}

template <typename T, qualifier Q>
inline T rigid_body_t<T, Q>::boost_energy() const {
  velocity_t const p = local_momentum();
  return glm::dot(p.boost, velocity.boost);
}

template <typename T, qualifier Q>
inline T rigid_body_t<T, Q>::kinetic_energy() const {
  velocity_t const p = local_momentum();
  return glm::dot(p.rotational, velocity.rotational) +
         glm::dot(p.boost, velocity.boost);
}

template <typename T, qualifier Q>
inline typename rigid_body_t<T, Q>::velocity_t rigid_body_t<T, Q>::local_momentum() const {
  return apply_inertia(velocity);
}

template <typename T, qualifier Q>
inline typename rigid_body_t<T, Q>::velocity_t rigid_body_t<T, Q>::momentum() const {
  return local_momentum().conjugated_by(frame);
}

template <typename T, qualifier Q>
inline typename rigid_body_t<T, Q>::velocity_t
rigid_body_t<T, Q>::calculate_free_acceleration() const {

  velocity_t const Ma = bracket(local_momentum(), velocity);

  return apply_inverse_inertia(Ma);
}

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q> mass_distribution_matrix_point(glm::vec<4, T, Q> const &position, T mass) {
  return mass * glm::lorentz::outerProduct(position, position);
}

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q> mass_distribution_matrix_sphere_shell(T mass, T radius) {
  glm::mat<4, 4, T, Q> M{};
  T a = mass * std::pow(std::sinh(radius), 2) / 3;
  M[0][0] = M[1][1] = M[2][2] = a;
  M[3][3] = -mass - 3 * a;
  return M;
}

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q> mass_distribution_matrix_sphere_shell(glm::vec<4, T, Q> const &position,
                                                  T mass, T radius) {

  using mat4_t = glm::mat<4, 4, T, Q>;
  mat4_t P = hyperbolic::transport_frame_to(mat4_t{1}, position);
  return P * mass_distribution_matrix_sphere_shell(mass, radius) *
         glm::transpose(P);
}

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q> mass_distribution_matrix(rigid_body_t<T, Q> const &b) {
  using vec4_t = glm::vec<4, T, Q>;

  glm::mat<4, 4, T, Q> M{0, 1, 1, -1, 1, 0, 1, -1, 1, 1, 0, -1, 0, 0, 0, -1};

  vec4_t const diagonal =
      glm::inverse(M) * vec4_t{b.distributional_inertia, b.total_mass};

  return b.frame * glm::diagonal4x4(diagonal) *
         glm::lorentz::transpose(b.frame);
}

template <typename T, qualifier Q>
inline void rigid_body_t<T, Q>::free_motion(T dt) {

  // Euler method
  frame = frame * exponential(velocity * dt);
  velocity += calculate_free_acceleration() * dt;
}

template <typename T, qualifier Q>
inline typename rigid_body_t<T, Q>::velocity_t
rigid_body_t<T, Q>::apply_inertia(velocity_t const &v) const {
  vec3_t const &m = distributional_inertia;
  return velocity_t{m * v.rotational, (m + total_mass) * v.boost};
}

template <typename T, qualifier Q>
inline typename rigid_body_t<T, Q>::velocity_t
rigid_body_t<T, Q>::apply_inverse_inertia(velocity_t const &p) const {
  vec3_t const &m = distributional_inertia;
  return velocity_t{p.rotational / m, p.boost / (m + total_mass)};
}

template <typename T, qualifier Q>
inline void rigid_body_t<T, Q>::apply_impulse(vec4_t const &location,
                                        vec4_t const &direction) {
  vec4_t local_location = frame * location;
  vec4_t local_direction = frame * direction;

  apply_local_impulse(local_location, local_direction);
}

template <typename T, qualifier Q>
inline void rigid_body_t<T, Q>::apply_local_impulse(vec4_t const &location,
                                              vec4_t const &direction) {
  // Convert impulse into lie algebra element
  velocity_t imp = calculate_local_impulse(location, direction);

  // Convert back to velocity
  velocity += apply_inverse_inertia(imp);
}

template <typename T, qualifier Q>
inline typename rigid_body_t<T, Q>::velocity_t
rigid_body_t<T, Q>::calculate_local_impulse(vec4_t const &location,
                                      vec4_t const &direction) {
  return wedge(location, direction) * static_cast<T>(0.5);
}

template <typename T, qualifier Q>
inline std::array<int, 2> find_pivot(glm::mat<4, 4, T, Q> const &M) {
  using pivot_t = std::array<int, 2>;
  pivot_t pivot{0, 1};

  T error = abs(M[pivot[0]][pivot[1]]);

  constexpr std::array<pivot_t, 5> candidates = {pivot_t{0, 2}, pivot_t{0, 3},
                                                 pivot_t{1, 2}, pivot_t{1, 3},
                                                 pivot_t{2, 3}};
  for (auto const &p : candidates) {
    auto v = abs(M[p[0]][p[1]]);
    if (v > error) {
      pivot = p;
      error = v;
    }
  }
  return pivot;
}

template <typename T, qualifier Q>
inline rigid_body_t<T, Q>
rigid_body_t<T, Q>::from_distribution_matrix(mat4_t const &dist_matrix) {
  using T = mat4::value_type;
  using pivot_t = std::array<int, 2>;
  mat4 M = dist_matrix;
  mat4 P{1.0f};

  // Apply a Lorentzian version of the Jacobi eigenvalue algorithm to transform
  // M to a diagonal matrix
  while (1) {
    // Find the largest pivot
    pivot_t pivot = find_pivot(M);
    if (abs(M[pivot[0]][pivot[1]]) < 1e-6)
      continue;

    if (pivot[1] == 3) {
      // Find the hyperbolic angle
      T alpha = std::atanh(2 * M[pivot[0]][pivot[1]] /
                           (M[pivot[1]][pivot[1]] - M[pivot[0]][pivot[0]]));

      // Apply the Givens boost matrix

    } else {
      // Find the angle
      T alpha = std::atan2(2 * M[pivot[0]][pivot[1]],
                           M[pivot[1]][pivot[1]] - M[pivot[0]][pivot[0]]);

      // Apply the Givens rotation matrix
    }
  }

  rigid_body_t result;
  result.frame = P;
  result.total_mass = -(dist_matrix[0][0] + dist_matrix[1][1]) -
                      (dist_matrix[2][2] + dist_matrix[3][3]);
  result.distributional_inertia =
      vec3{M[1][1] + M[2][2], M[0][0] + M[2][2], M[0][0] + M[1][1]};
  return result;
}

template <typename T, qualifier Q>
inline T
rigid_body_t<T, Q>::impulse_magnitude_from_elastic_collision(rigid_body_t<T, Q> const &b1,
                                                       rigid_body_t<T, Q> const &b2,
                                                       vec4_t const &p,
                                                       vec4_t const &n) {

  velocity_t const delta_v = b1.velocity - b2.velocity;
  T coeff1 = glm::lorentz::dot(n, delta_v * p);

  velocity_t const c = b1.apply_inverse_inertia(wedge(p, n)) +
                       b2.apply_inverse_inertia(wedge(p, n));
  T coeff2 = glm::lorentz::dot(n, c * p) / 2;

  return -coeff1 / coeff2;
}

template <typename T, qualifier Q>
inline T
rigid_body_t<T, Q>::impulse_magnitude_from_inelastic_collision(
    rigid_body_t const &b1, rigid_body_t const &b2, vec4_t const &p,
    vec4_t const &n, T energy_loss_factor) {

  velocity_t const delta_v = b1.velocity - b2.velocity;
  T const coeff1 = glm::lorentz::dot(n, delta_v * p);

  velocity_t const c = b1.apply_inverse_inertia(wedge(p, n)) +
                       b2.apply_inverse_inertia(wedge(p, n));
  T const coeff2 = glm::lorentz::dot(n, c * p) / 2;

  T const coeff0 =
      energy_loss_factor * (b1.kinetic_energy() + b2.kinetic_energy());

  T const disc = coeff1 * coeff1 - 4 * coeff0 * coeff2;

  if (disc >= 0)
    return (-coeff1 + std::sqrt(disc)) / (2 * coeff2);

  // If no collision exists that satisfies the loss then return the minimum
  // energy as the best attempt
  return -coeff1 / (2 * coeff2);
}

} // namespace hyperbolic