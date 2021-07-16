#pragma once

#include <array>
#include <math.h>

#include "vector_types.h"
#include <glm/glm.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include "lorentz_lie_algebra.h"

#include "hyperbolic.h"


struct rigid_body_t {
  using velocity_t = lorentz_lie_algebra_t;
  using scalar_t = vec3::value_type;

  mat4 frame{1.0f};               // A Lorentz matrix
  velocity_t velocity{};
  vec3 distributional_inertia{1.0f};
  scalar_t total_mass{1};

  scalar_t rotational_energy() const;
  scalar_t boost_energy() const;
  scalar_t kinetic_energy() const;
  velocity_t local_momentum() const;
  velocity_t calculate_free_acceleration() const;
  void free_motion(mat4::value_type dt);
  velocity_t apply_inertia(velocity_t const& v) const;
  velocity_t apply_inverse_inertia(velocity_t const &p) const;

  void apply_impulse(vec4 const &location, vec4 const &direction);
  void apply_local_impulse(vec4 const &location, vec4 const &direction);
  static velocity_t calculate_local_impulse(vec4 const &location,
                                     vec4 const &direction);

  static rigid_body_t from_distribution_matrix(mat4 const &M);  // Returns a rigid body at rest
  static scalar_t
  impulse_magnitude_from_elastic_collision(rigid_body_t const &b1,
                                           rigid_body_t const &b2,
                                           vec4 const &p, vec4 const &n);
  static scalar_t impulse_magnitude_from_inelastic_collision(
      rigid_body_t const &b1, rigid_body_t const &b2, vec4 const &p,
      vec4 const &n, scalar_t energy_loss_factor);
};


inline rigid_body_t::scalar_t rigid_body_t::rotational_energy() const {
  velocity_t p = local_momentum();
  return glm::dot(p.rotational, velocity.rotational);
}
inline rigid_body_t::scalar_t rigid_body_t::boost_energy() const {
  velocity_t p = local_momentum();
  return glm::dot(p.boost, velocity.boost);
}
inline rigid_body_t::scalar_t rigid_body_t::kinetic_energy() const {
  velocity_t p = local_momentum();
  return glm::dot(p.rotational, velocity.rotational) +
         glm::dot(p.boost, velocity.boost);
}

inline rigid_body_t::velocity_t rigid_body_t::local_momentum() const {
  return apply_inertia(velocity);
}

inline rigid_body_t::velocity_t
rigid_body_t::calculate_free_acceleration() const {

  velocity_t const Ma = bracket(velocity, local_momentum());

  return apply_inverse_inertia(Ma);
}

inline mat4 mass_distribution_matrix_point(vec4 const &position, float mass) {
  return mass * glm::lorentz::outerProduct(position, position);
}

inline mat4 mass_distribution_matrix_sphere_shell(float mass, float radius) {
  mat4 M = mat4{};
  mat4::value_type a = mass * powf(sinh(radius), 2) / 3;
  M[0][0] = M[1][1] = M[2][2] = a;
  M[3][3] = -mass - 3 * a;
  return M;
}

inline mat4 mass_distribution_matrix_sphere_shell(vec4 const &position,
                                                  float mass, float radius) {
  mat4 P = glm::hyperbolic::transport_frame_to(mat4{1}, {position});
  return P * mass_distribution_matrix_sphere_shell(mass, radius) * glm::transpose(P);
}

inline mat4 mass_distribution_matrix(rigid_body_t const& b) {

  mat4 constexpr M{0, 1, 1, -1, 1, 0, 1, -1, 1, 1, 0, -1, 0, 0, 0, -1};
  vec4 const diagonal = glm::inverse(M) * vec4{b.distributional_inertia, b.total_mass};
  
  return b.frame * glm::diagonal4x4(diagonal) *
         glm::lorentz::transpose(b.frame);
}

inline void rigid_body_t::free_motion(mat4::value_type dt) {

  velocity += calculate_free_acceleration() * dt;

  frame = frame * exponential(velocity * dt);
}

inline rigid_body_t::velocity_t
rigid_body_t::apply_inertia(velocity_t const &v) const {
  vec3 const &m = distributional_inertia;
  return velocity_t{m * v.rotational, (m + total_mass) * v.boost};
}

inline rigid_body_t::velocity_t
rigid_body_t::apply_inverse_inertia(velocity_t const &p) const {
  vec3 const &m = distributional_inertia;
  return velocity_t{p.rotational / m, p.boost / (m + total_mass)};
}

inline void rigid_body_t::apply_impulse(vec4 const &location,
                                        vec4 const &direction) {
  vec4 local_location = frame * location;
  vec4 local_direction = frame * direction;

  apply_local_impulse(local_location, local_direction);
}

inline void rigid_body_t::apply_local_impulse(vec4 const &location,
                                              vec4 const &direction) {
  // Convert impulse into lie algebra element
  lorentz_lie_algebra_t imp = calculate_local_impulse(location, direction);

  // Convert back to velocity
  velocity += apply_inverse_inertia(imp);
}

inline rigid_body_t::velocity_t
rigid_body_t::calculate_local_impulse(vec4 const &location,
                                      vec4 const &direction) {
  return wedge(location, direction) * static_cast<scalar_t>(0.5);
}

inline std::array<int, 2> find_pivot(mat4 const &M) {
  using pivot_t = std::array<int, 2>;
  pivot_t pivot{0, 1};

  mat4::value_type error = abs(M[pivot[0]][pivot[1]]);

  constexpr std::array<pivot_t, 5> candidates = {
      pivot_t{0, 2}, pivot_t{0, 3}, pivot_t{1, 2}, pivot_t{1, 3}, pivot_t{2, 3}};
  for (auto const& p : candidates) {
    auto v = abs(M[p[0]][p[1]]);
    if (v > error) {
      pivot = p;
      error = v;
    }
  }
  return pivot;
}

inline rigid_body_t rigid_body_t::from_distribution_matrix(mat4 const &dist_matrix) {
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
      T alpha = std::atan2(2 * M[pivot[0]][pivot[1]], M[pivot[1]][pivot[1]] - M[pivot[0]][pivot[0]]);

      // Apply the Givens rotation matrix
    }
  }

  rigid_body_t result;
  result.frame = P;
  result.total_mass =
      -(dist_matrix[0][0] + dist_matrix[1][1]) - (dist_matrix[2][2] + dist_matrix[3][3]);
  result.distributional_inertia =
      vec3{M[1][1] + M[2][2], M[0][0] + M[2][2], M[0][0] + M[1][1]};
  return result;
}

inline rigid_body_t::scalar_t
rigid_body_t::impulse_magnitude_from_elastic_collision(rigid_body_t const &b1,
                                                       rigid_body_t const &b2,
                                                       vec4 const &p,
                                                       vec4 const &n) {

  velocity_t const delta_v = b1.velocity - b2.velocity;
  scalar_t coeff1 = glm::lorentz::dot(n, delta_v * p);

  velocity_t const c = b1.apply_inverse_inertia(wedge(p, n)) +
                       b2.apply_inverse_inertia(wedge(p, n));
  scalar_t coeff2 = glm::lorentz::dot(n, c * p) / 2;

  return -coeff1 / coeff2;
}

inline rigid_body_t::scalar_t
rigid_body_t::impulse_magnitude_from_inelastic_collision(
    rigid_body_t const &b1, rigid_body_t const &b2, vec4 const &p,
    vec4 const &n, scalar_t energy_loss_factor) {

  velocity_t const delta_v = b1.velocity - b2.velocity;
  scalar_t const coeff1 = glm::lorentz::dot(n, delta_v * p);

  velocity_t const c = b1.apply_inverse_inertia(wedge(p, n)) +
                       b2.apply_inverse_inertia(wedge(p, n));
  scalar_t const coeff2 = glm::lorentz::dot(n, c * p) / 2;

  scalar_t const coeff0 =
      energy_loss_factor * (b1.kinetic_energy() + b2.kinetic_energy());

  scalar_t const disc = coeff1 * coeff1 - 4 * coeff0 * coeff2;

  if (disc >= 0)
    return (-coeff1 + std::sqrt(disc)) / (2 * coeff2);

  // If no collision exists that satisfies the loss then return the minimum
  // energy as the best attempt
  return -coeff1 / (2 * coeff2);
}