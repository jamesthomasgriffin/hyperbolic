#pragma once

#include <array>

#include <glm/glm.hpp>

#include "hyperbolic.h"
#include "lorentz_lie_algebra.h"

namespace hyperbolic {

template <typename T, qualifier Q = glm::highp> struct rigid_body_t {
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
  velocity_t global_velocity() const;
  velocity_t local_momentum() const;
  velocity_t momentum() const;
  velocity_t calculate_free_acceleration() const;
  void free_motion(T dt);
  velocity_t apply_inertia(velocity_t const &v) const;
  velocity_t apply_inverse_inertia(velocity_t const &p) const;

  void apply_impulse(vec4_t const &location, vec4_t const &direction);
  void apply_impulse(velocity_t const &impulse);
  void apply_local_impulse(vec4_t const &location, vec4_t const &direction);
  void apply_local_impulse(velocity_t const &impulse);
  static velocity_t calculate_local_impulse(vec4_t const &location,
                                            vec4_t const &direction);

  static rigid_body_t
  from_distribution_matrix(mat4_t const &M); // Returns a rigid body at rest
  static rigid_body_t from_diagonal_distribution_matrix(
      vec4_t const &diagonal); // Returns a rigid body at rest

  static scalar_t
  impulse_magnitude_from_elastic_collision(rigid_body_t const &b1,
                                           rigid_body_t const &b2,
                                           vec4_t const &p, vec4_t const &n);
  static scalar_t
  impulse_magnitude_from_elastic_collision(rigid_body_t const &b1,
                                           rigid_body_t const &b2,
                                           velocity_t const &imp_dir);
  static scalar_t impulse_magnitude_from_elastic_collision(
      rigid_body_t const &b1, rigid_body_t const &b2, mat4_t const &rel_frame,
      velocity_t const &imp_dir);

  static scalar_t
  impulse_magnitude_from_elastic_collision(rigid_body_t const &b,
                                           velocity_t const &imp_dir);

  static scalar_t impulse_magnitude_from_inelastic_collision(
      rigid_body_t const &b1, rigid_body_t const &b2, vec4_t const &p,
      vec4_t const &n, T energy_loss_factor);
  static scalar_t impulse_magnitude_from_inelastic_collision(
      rigid_body_t<T, Q> const &b1, rigid_body_t<T, Q> const &b2,
      velocity_t const &imp_dir, T energy_loss_factor);
};

//template <typename T, qualifier Q>
//inline T rigid_body_t<T, Q>::rotational_energy() const {
//  velocity_t const p = local_momentum();
//  return static_cast<T>(-0.5) * glm::dot(p.rotational, velocity.rotational);
//}
//
//template <typename T, qualifier Q>
//inline T rigid_body_t<T, Q>::boost_energy() const {
//  velocity_t const p = local_momentum();
//  return static_cast<T>(-0.5) * glm::dot(p.boost, velocity.boost);
//}

// The kinetic energy of a moving frame is -1/2 * tr(V^t P), where V
// is the velocity, and P is the momentum M(V) (for mass tensor M).
template <typename T, qualifier Q>
inline T rigid_body_t<T, Q>::kinetic_energy() const {
  velocity_t const p = local_momentum();
  return static_cast<T>(-0.5) * hyperbolic::dot(p, velocity);
}

template <typename T, qualifier Q>
inline typename rigid_body_t<T, Q>::velocity_t
rigid_body_t<T, Q>::local_momentum() const {
  return apply_inertia(velocity);
}

template <typename T, qualifier Q>
inline typename rigid_body_t<T, Q>::velocity_t
rigid_body_t<T, Q>::global_velocity() const {
  return velocity.conjugated_by(frame);
}

template <typename T, qualifier Q>
inline typename rigid_body_t<T, Q>::velocity_t
rigid_body_t<T, Q>::momentum() const {
  return local_momentum().conjugated_by(frame);
}

template <typename T, qualifier Q>
inline typename rigid_body_t<T, Q>::velocity_t
rigid_body_t<T, Q>::calculate_free_acceleration() const {

  velocity_t const Ma = bracket(local_momentum(), velocity);

  return apply_inverse_inertia(Ma);
}

namespace mass_distribution {

// The mass distribution M is the integral of
// -density(p)pp^T
// where p in R^4 ranges over the points of the object in the hyperboloid model
// and p^T is the Lorentzian transpose of p. 
//  
// Note that the trace of M is the integral of 
// -density(p)trace(pp^T) = -density(p) * p.p = density(p),
// i.e. the trace of M is the mass of the object.
// 
// M is symmetric for the Lorentzian transpose and furthermore satisfies the
// property that q^T M q < 0 for all q on the hyperboloid.  This is enough to
// guarantee that M is diagonalisable via a Lorentzian orthogonal matrix
// with three negative and one positive eigenvalues.
//
// Hence it is possible to always find a frame such that M is diagonal.
// Assuming that M is diagonal with values (a1, a2, a3, b), we define
// the directional inertia (d.i.) vector to be -(a2 + a3, a1 + a3, a1 + a2),
// note that since a1, a2 and a3 are negative the d.i. vector has positive
// entries.  These values along with the mass m = a1 + a2 + a3 + b, allow
// the velocity of a rigid frame to be transformed to a momentum.
//

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q> of_point(glm::vec<4, T, Q> const &position,
                                     T mass) {
  return (-mass) * glm::lorentz::outerProduct(position, position);
}

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q> of_sphere_shell(T mass, T radius) {
  glm::mat<4, 4, T, Q> M{};
  T a = -mass * std::pow(std::sinh(radius), 2) / 3;
  M[0][0] = M[1][1] = M[2][2] = a;
  M[3][3] = mass - 3 * a;
  return M;
}

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q> of_sphere_shell(glm::vec<4, T, Q> const &position,
                                            T mass, T radius) {
  using mat4_t = glm::mat<4, 4, T, Q>;
  mat4_t P = hyperbolic::transport_frame_to(mat4_t{1}, position);
  return P * mass_distribution_matrix_sphere_shell(mass, radius) *
         glm::transpose(P);
}

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q> of_rigid_body(rigid_body_t<T, Q> const &b) {
  using vec3_t = glm::vec<3, T, Q>;
  using vec4_t = glm::vec<4, T, Q>;
  using mat4_t = glm::mat<4, 4, T, Q>;

  vec3_t const &di = b.distributional_inertia;
  T sum_3_evals = static_cast<T>(-0.5) * (di.x + di.y + di.z);

  vec4_t const diag{vec3_t{sum_3_evals} + di, b.total_mass - sum_3_evals};

  mat4_t const &F = b.frame;
  mat4_t const FD{diag[0] * F[0], diag[1] * F[1], diag[2] * F[2],
                  diag[3] * F[3]};
  mat4_t const Ft = glm::lorentz::transpose(b.frame);
  mat4_t const FDFt = FD * Ft;

  return FDFt;
}

} // namespace mass_distribution

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
inline void rigid_body_t<T, Q>::apply_impulse(velocity_t const& impulse) {
  velocity_t local_impulse =
      impulse.conjugated_by(glm::lorentz::transpose(frame));

  apply_local_impulse(local_impulse);
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
inline void rigid_body_t<T, Q>::apply_local_impulse(velocity_t const& impulse) {
  // Convert back to velocity
  velocity += apply_inverse_inertia(impulse);
}

template <typename T, qualifier Q>
inline typename rigid_body_t<T, Q>::velocity_t
rigid_body_t<T, Q>::calculate_local_impulse(vec4_t const &location,
                                            vec4_t const &direction) {
  return wedge(location, direction) * static_cast<T>(0.5);
}

// UNTESTED
//template <typename T, qualifier Q>
//inline std::array<int, 2> find_pivot(glm::mat<4, 4, T, Q> const &M) {
//  using pivot_t = std::array<int, 2>;
//  pivot_t pivot{0, 1};
//
//  T error = abs(M[pivot[0]][pivot[1]]);
//
//  constexpr std::array<pivot_t, 5> candidates{pivot_t{0, 2}, pivot_t{0, 3},
//                                                 pivot_t{1, 2}, pivot_t{1, 3},
//                                                 pivot_t{2, 3}};
//  for (auto const &p : candidates) {
//    auto v = abs(M[p[0]][p[1]]);
//    if (v > error) {
//      pivot = p;
//      error = v;
//    }
//  }
//  return pivot;
//}

// UNTESTED
template <typename T, qualifier Q>
inline rigid_body_t<T, Q>
rigid_body_t<T, Q>::from_diagonal_distribution_matrix(glm::vec<4, T, Q> const &diag) {
  rigid_body_t<T, Q> result{};
  result.total_mass = -((diag[0] + diag[1]) + (diag[2] + diag[3]));
  result.distributional_inertia =
      glm::vec<3, T, Q>{diag[1] + diag[2], diag[0] + diag[2], diag[0] + diag[1]};
  return result;
}

// INCOMPLETE
template <typename T, qualifier Q>
inline rigid_body_t<T, Q>
rigid_body_t<T, Q>::from_distribution_matrix(mat4_t const &dist_matrix) {
  using pivot_t = std::array<int, 2>;
  glm::mat<4, 4, T, Q> M = dist_matrix;
  glm::mat<4, 4, T, Q> P{1.0f};

  // Apply a Lorentzian version of the Jacobi eigenvalue algorithm to transform
  // M to a diagonal matrix
  while (1) {
    // Find the largest pivot
    pivot_t pivot = find_pivot(M);
    if (abs(M[pivot[0]][pivot[1]]) < static_cast<T>(1e-6))
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
      glm::vec<3, T, Q>{M[1][1] + M[2][2], M[0][0] + M[2][2], M[0][0] + M[1][1]};
  return result;
}

// Given the direction of an impulse applied to two bodies, this returns the non-zero factor
// to multiply that direction by to get an energy preserving impulse.
// The impulse direction is in the frame of b1.
template <typename T, qualifier Q>
inline T rigid_body_t<T, Q>::impulse_magnitude_from_elastic_collision(
    rigid_body_t<T, Q> const &b1, rigid_body_t<T, Q> const &b2, velocity_t const& imp_dir) {

  mat4_t const frame1_from_frame2 = glm::lorentz::transpose(b1.frame) * b2.frame;
  return impulse_magnitude_from_elastic_collision(b1, b2, frame1_from_frame2,
                                                  imp_dir);
}

// Given the direction of an impulse applied to two bodies, this returns the
// non-zero factor to multiply that direction by to get an energy preserving
// impulse. The impulse direction is in the frame of b1.
// For this variant of the function, the relative frame of b2 to b1 is provided.
// This is to allow for carrying out the task even when b1 and b2 are in different frames.
template <typename T, qualifier Q>
inline T rigid_body_t<T, Q>::impulse_magnitude_from_elastic_collision(
    rigid_body_t<T, Q> const &b1, rigid_body_t<T, Q> const &b2, mat4_t const& rel_frame,
    velocity_t const &imp_dir) {

  // Delta energy has form
  // a * imp_dir . delta_v - 0.5 * a^2 imp_dir (M1inv + M2inv) imp_dir
  // = a * coeff1 - a*a * coeff2
  // where a is a constant we multiply the given impulse by, we solve for a
  // non-zero
  mat4_t const& frame1_from_frame2 = rel_frame;
  mat4_t const frame2_from_frame1 = glm::lorentz::transpose(rel_frame);
  velocity_t const delta_v =
      b2.velocity.conjugated_by(frame1_from_frame2) - b1.velocity;
  T const coeff1 = dot(imp_dir, delta_v);

  velocity_t const imp_in_frame2 = imp_dir.conjugated_by(frame2_from_frame1);
  velocity_t const c =
      b1.apply_inverse_inertia(imp_dir) +
      b2.apply_inverse_inertia(imp_in_frame2).conjugated_by(frame1_from_frame2);
  T const coeff2 = dot(c, imp_dir) / 2;

  return coeff1 / coeff2;
}

template <typename T, qualifier Q>
inline T rigid_body_t<T, Q>::impulse_magnitude_from_elastic_collision(
    rigid_body_t<T, Q> const &b1, rigid_body_t<T, Q> const &b2, vec4_t const &p,
    vec4_t const &n) {

  velocity_t impulse = calculate_local_impulse(p, n);
  return impulse_magnitude_from_elastic_collision(b1, b2, impulse);
}

// Collision of a rigid body with a massive object, preserving energy of the rigid body
template <typename T, qualifier Q>
inline T rigid_body_t<T, Q>::impulse_magnitude_from_elastic_collision(
  rigid_body_t const& b, velocity_t const& imp_dir) {
  T const coeff1 = dot(b.velocity, imp_dir);
  T const coeff2 =
      static_cast<T>(0.5) * dot(imp_dir, b.apply_inverse_inertia(imp_dir));
  return -coeff1 / coeff2;
}

// UNTESTED
template <typename T, qualifier Q>
inline T rigid_body_t<T, Q>::impulse_magnitude_from_inelastic_collision(
    rigid_body_t const &b1, rigid_body_t const &b2, vec4_t const &p,
    vec4_t const &n, T energy_loss_factor) {

  velocity_t impulse = calculate_local_impulse(p, n);
  return impulse_magnitude_from_inelastic_collision(b1, b2, impulse, energy_loss_factor);
}

// UNTESTED
//template <typename T, qualifier Q>
//inline T rigid_body_t<T, Q>::impulse_magnitude_from_inelastic_collision(
//    rigid_body_t<T, Q> const &b1, rigid_body_t<T, Q> const &b2,
//    velocity_t const &imp_dir, T energy_loss_factor) {
//
//  velocity_t const delta_v = b1.velocity - b2.velocity;
//  T const coeff1 = glm::lorentz::dot(n, delta_v * p);
//
//  velocity_t const c = b1.apply_inverse_inertia(imp_dir) +
//                       b2.apply_inverse_inertia(imp_dir);
//  T const coeff2 = glm::lorentz::dot(n, c * p) / 2;
//
//  T const coeff0 =
//      energy_loss_factor * (b1.kinetic_energy() + b2.kinetic_energy());
//
//  T const disc = coeff1 * coeff1 - 4 * coeff0 * coeff2;
//
//  if (disc >= 0)
//    return (-coeff1 + glm::sqrt(disc)) / (2 * coeff2);
//
//  // If no collision exists that satisfies the loss then return the minimum
//  // energy as the best attempt
//  return -coeff1 / (2 * coeff2);
//}

} // namespace hyperbolic
