#pragma once

#include <complex>

#include <glm/glm.hpp>

#include "vector_types.h"
#include "glm_lorentz.h"
#include "glm_traits.h"

namespace hyperbolic {

struct lorentz_lie_algebra_t {
  using scalar_t = typename vec3::value_type;
  using vec3_t = vec3;
  using vec4_t = vec4;
  using mat4_t = mat4;
  vec3_t rotational{0.0f}; // Anti-symmetric part of Lie algebra
  vec3_t boost{0.0f};      // Symmetric part of Lie algebra

  bool operator==(lorentz_lie_algebra_t const &b) const;
  bool operator!=(lorentz_lie_algebra_t const &b) const;
  lorentz_lie_algebra_t operator+(lorentz_lie_algebra_t const &b) const;
  lorentz_lie_algebra_t operator-(lorentz_lie_algebra_t const &b) const;
  lorentz_lie_algebra_t operator*(scalar_t c) const;
  vec4 operator*(vec4 const &v) const;

  lorentz_lie_algebra_t &operator+=(lorentz_lie_algebra_t const &b);
  lorentz_lie_algebra_t &operator-=(lorentz_lie_algebra_t const &b);
  lorentz_lie_algebra_t &operator*=(scalar_t c);

  lorentz_lie_algebra_t conjugated_by(mat4 const &P) const;
  lorentz_lie_algebra_t conjugated_by(mat3 const &P) const;
  mat4 to_matrix() const;
  static lorentz_lie_algebra_t from_matrix(mat4 const &A);
};

inline lorentz_lie_algebra_t operator*(lorentz_lie_algebra_t::scalar_t c,
                                       lorentz_lie_algebra_t const &b) {
  return b * c;
}

inline lorentz_lie_algebra_t bracket(lorentz_lie_algebra_t const &a,
                                     lorentz_lie_algebra_t const &b) {
  return {glm::cross(a.rotational, b.rotational) - glm::cross(a.boost, b.boost),
          glm::cross(a.rotational, b.boost) +
              glm::cross(a.boost, b.rotational)};
}

// Result of v'w - w'v
inline lorentz_lie_algebra_t wedge(vec4 const &v, vec4 const &w) {
  return {glm::cross(vec3{w}, vec3{v}), vec3{w[3] * v - v[3] * w}};
}

inline bool
lorentz_lie_algebra_t::operator==(lorentz_lie_algebra_t const &b) const {
  return (rotational == b.rotational) && (boost == b.boost);
}

inline bool
lorentz_lie_algebra_t::operator!=(lorentz_lie_algebra_t const &b) const {
  return (rotational != b.rotational) || (boost != b.boost);
}

inline lorentz_lie_algebra_t
lorentz_lie_algebra_t::operator+(lorentz_lie_algebra_t const &b) const {
  return {rotational + b.rotational, boost + b.boost};
}

inline lorentz_lie_algebra_t
lorentz_lie_algebra_t::operator-(lorentz_lie_algebra_t const &b) const {
  return {rotational - b.rotational, boost - b.boost};
}

inline lorentz_lie_algebra_t
lorentz_lie_algebra_t::operator*(scalar_t c) const {
  return {rotational * c, boost * c};
}

inline vec4 lorentz_lie_algebra_t::operator*(vec4 const &v) const {
  return vec4{glm::cross(rotational, vec3{v}) + v[3] * boost,
              glm::dot(boost, vec3{v})};
}

inline lorentz_lie_algebra_t &
lorentz_lie_algebra_t::operator+=(lorentz_lie_algebra_t const &b) {
  rotational += b.rotational;
  boost += b.boost;
  return *this;
}

inline lorentz_lie_algebra_t &
lorentz_lie_algebra_t::operator-=(lorentz_lie_algebra_t const &b) {
  rotational -= b.rotational;
  boost -= b.boost;
  return *this;
}

inline lorentz_lie_algebra_t &lorentz_lie_algebra_t::operator*=(scalar_t c) {
  rotational *= c;
  boost *= c;
  return *this;
}

inline lorentz_lie_algebra_t
lorentz_lie_algebra_t::conjugated_by(mat4 const &P) const {
  mat4 M = P * to_matrix() * glm::lorentz::transpose(P);
  return from_matrix(M);
}

inline lorentz_lie_algebra_t
lorentz_lie_algebra_t::conjugated_by(mat3 const &P) const {
  return lorentz_lie_algebra_t{P * rotational, P * boost};
}

inline mat4 lorentz_lie_algebra_t::to_matrix() const {
  mat4 result{0};
  result[3] = vec4{boost, 0};
  result[0][3] = boost[0];
  result[1][3] = boost[1];
  result[2][3] = boost[2];
  result[0][1] = rotational[2];
  result[1][0] = -rotational[2];
  result[1][2] = rotational[0];
  result[2][1] = -rotational[0];
  result[0][2] = -rotational[1];
  result[2][0] = rotational[1];
  return result;
}

inline lorentz_lie_algebra_t lorentz_lie_algebra_t::from_matrix(mat4 const &A) {
  return lorentz_lie_algebra_t{vec3{A[1][2], A[2][0], A[0][1]}, vec3{A[3]}};
}

namespace detail {

template<typename M2>
inline std::complex<typename M2::value_type> det(ImaginaryExtension<M2> const &A) {
  using T = typename mat2::value_type;
  std::complex<T> a{A.real[0][0], A.imag[0][0]};
  std::complex<T> b{A.real[1][0], A.imag[1][0]};
  std::complex<T> c{A.real[0][1], A.imag[0][1]};
  std::complex<T> d{A.real[1][1], A.imag[1][1]};
  return a * d - b * c;
}
  
template<typename M2>
struct sl2_lie_algebra {
  ImaginaryExtension<M2> matrix{M2{1}};

  static sl2_lie_algebra from(lorentz_lie_algebra_t const &v) {
    using vec_t = typename glm::traits::to_vec<3, M2>::type;
    vec_t const &w = v.rotational;
    vec_t const &b = v.boost;

    return {ImaginaryExtension<M2>{
        0.5f * M2{b[2], b[0] + w[1], b[0] - w[1], -b[2]},
        0.5f * M2{w[2], w[0] - b[1], w[0] + b[1], -w[2]}}};
  };
};

template <typename M2>
inline ImaginaryExtension<M2> exponential(sl2_lie_algebra<M2> const& A) {
  using T = typename M2::value_type;

  std::complex<T> const detA = det(A.matrix);

  std::complex<T> const a =
      std::sqrt(detA); // Below values are independent of chosen root

  std::complex<T> const sinc_a = (std::norm(a) < 1e-12) ? a : std::sin(a) / a;
  std::complex<T> const cos_a = std::cos(a);

  // exp(A) = cos(a) I + sin(a) / a A
  ImaginaryExtension<M2> expA =
      ImaginaryExtension<M2>{M2(cos_a.real()), M2(cos_a.imag())} +
      ImaginaryExtension<T>{sinc_a.real(), sinc_a.imag()} * A.matrix;

  return expA;
}

template<typename M2>
inline typename glm::traits::to_vec<4, M2>::type hermitian_matrix_to_vec(ImaginaryExtension<M2> const &M) {
  using vec_t = typename glm::traits::to_vec<4, M2>::type;
  vec_t result{};
  result[0] = M.real[0][1];
  result[1] = M.imag[1][0];
  result[2] = (M.real[0][0] - M.real[1][1]) / 2;
  result[3] = (M.real[0][0] + M.real[1][1]) / 2;
  return result;
}

// a simplification of hermitian_matrix_to_vec for when M is real and symmetric
template <typename M2>
inline typename glm::traits::to_vec<4, M2>::type
symmetric_matrix_to_vec(M2 const &M) {
  using vec_t = typename glm::traits::to_vec<4, M2>::type;
  vec_t result{};
  result[0] = M[0][1];
  result[1] = 0;
  result[2] = (M[0][0] - M[1][1]) / 2;
  result[3] = (M[0][0] + M[1][1]) / 2;
  return result;
}

// This converts a 2x2 complex matrix to a 4x4 real matrix representing the
// action of M on v = (x, y, z, w) by M H(v) M*, where H(v) is the Hermitian
// matrix:
//
// | w +  z  x - iy |
// | x + iy  w -  z |
//
// and M* is the conjugate transpose of M
template <typename M2>
inline typename glm::traits::to_mat<4, 4, M2>::type
action_via_hermitian_matrix(ImaginaryExtension<M2> const &M) {

  ImaginaryExtension<M2> const Mct{glm::transpose(M.real),
                                     -glm::transpose(M.imag)};
  using return_t = typename glm::traits::to_mat<4, 4, M2>::type;

  // NB Ew is the identity so isn't needed
  ImaginaryExtension<M2> Ex{M2{0, 1, 1, 0}, M2{0, 0, 0, 0}};
  ImaginaryExtension<M2> Ey{M2{0, 0, 0, 0}, M2{0, -1, 1, 0}};
  ImaginaryExtension<M2> Ez{M2{1, 0, 0, -1}, M2{0, 0, 0, 0}};

  return_t result{};
  result[0] = hermitian_matrix_to_vec(M * Ex * Mct);
  result[1] = hermitian_matrix_to_vec(M * Ey * Mct);
  result[2] = hermitian_matrix_to_vec(M * Ez * Mct);
  result[3] = hermitian_matrix_to_vec(M * Mct);
  return result;
}

// a simplification of action_via_hermitian_matrix for when M is real and
// symmetric
template <typename M2>
inline typename glm::traits::to_mat<4, 4, M2>::type action_via_symmetric_matrix(
    M2 const &M) {
  using return_t = typename glm::traits::to_mat<4, 4, M2>::type;

  M2 const Mct = glm::transpose(M);

  // NB Ew is the identity so isn't needed
  M2 Ex{0, 1, 1, 0};
  M2 imagEy{0, -1, 1, 0};
  M2 Ez{1, 0, 0, -1};

  return_t result{};
  result[0] = symmetric_matrix_to_vec(M * Ex * Mct);
  result[1] = vec4{0, (M * imagEy * Mct)[1][0], 0, 0};
  result[2] = symmetric_matrix_to_vec(M * Ez * Mct);
  result[3] = symmetric_matrix_to_vec(M * Mct);
  return result;
}

} // namespace detail

inline mat4 exponential(lorentz_lie_algebra_t const &v) {
  using sl2 = detail::sl2_lie_algebra<mat2>;
  using T = mat4::value_type;

  vec3 const &w = v.rotational;
  vec3 const &b = v.boost;

  sl2 A = sl2::from(v);

  ImaginaryExtension<mat2> expA = exponential(A);

  return detail::action_via_hermitian_matrix(expA);
}

} // namespace hyperbolic