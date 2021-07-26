#pragma once

#include <complex>

#include <glm/glm.hpp>

#include "glm_lorentz.h"

#include "imaginary_extension.h"

namespace hyperbolic {

template <typename T> using ImaginaryExtension = jtg::ImaginaryExtension<T>;

using qualifier = glm::qualifier;
using length_t = glm::length_t;

template <typename T, qualifier Q = glm::highp> struct lorentz_lie_algebra_t {
  using scalar_t = T;
  using vec3_t = glm::vec<3, T, Q>;
  using vec4_t = glm::vec<4, T, Q>;
  using mat4_t = glm::mat<4, 4, T, Q>;
  using mat3_t = glm::mat<3, 3, T, Q>;

  vec3_t rotational{0, 0, 0}; // Anti-symmetric part of Lie algebra
  vec3_t boost{0, 0, 0};      // Symmetric part of Lie algebra

  bool operator==(lorentz_lie_algebra_t const &b) const;
  bool operator!=(lorentz_lie_algebra_t const &b) const;
  lorentz_lie_algebra_t operator+(lorentz_lie_algebra_t const &b) const;
  lorentz_lie_algebra_t operator-(lorentz_lie_algebra_t const &b) const;
  lorentz_lie_algebra_t operator*(scalar_t c) const;
  lorentz_lie_algebra_t operator-() const;
  vec4_t operator*(vec4_t const &v) const;

  lorentz_lie_algebra_t &operator+=(lorentz_lie_algebra_t const &b);
  lorentz_lie_algebra_t &operator-=(lorentz_lie_algebra_t const &b);
  lorentz_lie_algebra_t &operator*=(scalar_t c);

  lorentz_lie_algebra_t conjugated_by(mat4_t const &P) const;
  lorentz_lie_algebra_t conjugated_by(mat3_t const &P) const;
  mat4_t to_matrix() const;
  static lorentz_lie_algebra_t from_matrix(mat4_t const &A);
};

// This is equivalent to the Frobenius inner product of matrices for the Lorentian
// inner product, i.e. tr(G^T H), where G and H are the respective matrix
// representations of lie algebra elements g and h and G^T is the Lorentzian
// transpose.
template <typename T, qualifier Q>
inline T dot(lorentz_lie_algebra_t<T, Q> const &g,
             lorentz_lie_algebra_t<T, Q> const &h) {
  return 2 * glm::dot(g.rotational, h.rotational) -
         2 * glm::dot(g.boost, h.boost);
}

template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q>
operator*(typename lorentz_lie_algebra_t<T, Q>::scalar_t c,
          lorentz_lie_algebra_t<T, Q> const &b) {
  return b * c;
}

template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q>
lorentz_lie_algebra_t<T, Q>::operator-() const {
  return {-rotational, -boost};
}

template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q>
bracket(lorentz_lie_algebra_t<T, Q> const &a,
        lorentz_lie_algebra_t<T, Q> const &b) {
  return {glm::cross(a.rotational, b.rotational) - glm::cross(a.boost, b.boost),
          glm::cross(a.rotational, b.boost) +
              glm::cross(a.boost, b.rotational)};
}

// Result of vw^t - wv^t
template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q> wedge(glm::vec<4, T, Q> const &v,
                                         glm::vec<4, T, Q> const &w) {
  using vec3_t = glm::vec<3, T, Q>;
  return {glm::cross(vec3_t{w}, vec3_t{v}), vec3_t{v[3] * w - w[3] * v}};
}

template <typename T, qualifier Q>
inline bool lorentz_lie_algebra_t<T, Q>::operator==(
    lorentz_lie_algebra_t<T, Q> const &b) const {
  return (rotational == b.rotational) && (boost == b.boost);
}

template <typename T, qualifier Q>
inline bool lorentz_lie_algebra_t<T, Q>::operator!=(
    lorentz_lie_algebra_t<T, Q> const &b) const {
  return (rotational != b.rotational) || (boost != b.boost);
}

template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q> lorentz_lie_algebra_t<T, Q>::operator+(
    lorentz_lie_algebra_t<T, Q> const &b) const {
  return {rotational + b.rotational, boost + b.boost};
}

template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q> lorentz_lie_algebra_t<T, Q>::operator-(
    lorentz_lie_algebra_t<T, Q> const &b) const {
  return {rotational - b.rotational, boost - b.boost};
}

template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q>
lorentz_lie_algebra_t<T, Q>::operator*(scalar_t c) const {
  return {rotational * c, boost * c};
}

template <typename T, qualifier Q>
inline glm::vec<4, T, Q>
lorentz_lie_algebra_t<T, Q>::operator*(vec4_t const &v) const {
  return vec4_t{glm::cross(rotational, vec3_t{v}) + v[3] * boost,
                glm::dot(boost, vec3_t{v})};
}

template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q> &
lorentz_lie_algebra_t<T, Q>::operator+=(lorentz_lie_algebra_t<T, Q> const &b) {
  rotational += b.rotational;
  boost += b.boost;
  return *this;
}

template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q> &
lorentz_lie_algebra_t<T, Q>::operator-=(lorentz_lie_algebra_t<T, Q> const &b) {
  rotational -= b.rotational;
  boost -= b.boost;
  return *this;
}

template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q> &
lorentz_lie_algebra_t<T, Q>::operator*=(scalar_t c) {
  rotational *= c;
  boost *= c;
  return *this;
}

template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q>
lorentz_lie_algebra_t<T, Q>::conjugated_by(mat4_t const &P) const {
  mat4_t M = P * to_matrix() * glm::lorentz::transpose(P);
  return from_matrix(M);
}

template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q>
lorentz_lie_algebra_t<T, Q>::conjugated_by(mat3_t const &P) const {
  return lorentz_lie_algebra_t{P * rotational, P * boost};
}

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q> lorentz_lie_algebra_t<T, Q>::to_matrix() const {
  mat4_t result{0};
  result[3] = vec4_t{boost, 0};
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

template <typename T, qualifier Q>
inline lorentz_lie_algebra_t<T, Q>
lorentz_lie_algebra_t<T, Q>::from_matrix(mat4_t const &A) {
  return lorentz_lie_algebra_t{vec3_t{A[1][2], A[2][0], A[0][1]}, vec3_t{A[3]}};
}

namespace detail {

template <typename T, qualifier Q>
inline std::complex<T> det(ImaginaryExtension<glm::mat<2, 2, T, Q>> const &A) {
  std::complex<T> a{A.real[0][0], A.imag[0][0]};
  std::complex<T> b{A.real[1][0], A.imag[1][0]};
  std::complex<T> c{A.real[0][1], A.imag[0][1]};
  std::complex<T> d{A.real[1][1], A.imag[1][1]};
  return a * d - b * c;
}

template <typename T, qualifier Q> struct sl2_lie_algebra {
  using vec3_t = glm::vec<3, T, Q>;
  using mat2_t = glm::mat<2, 2, T, Q>;

  ImaginaryExtension<mat2_t> matrix{mat2_t{1}};

  static sl2_lie_algebra from(lorentz_lie_algebra_t<float> const &v) {
    vec3_t const &w = v.rotational;
    vec3_t const &b = v.boost;

    return {ImaginaryExtension<glm::mat<2, 2, T, Q>>{
        0.5f * mat2_t{b[2], b[0] + w[1], b[0] - w[1], -b[2]},
        0.5f * mat2_t{w[2], w[0] - b[1], w[0] + b[1], -w[2]}}};
  };
};

template <typename T, qualifier Q>
inline ImaginaryExtension<glm::mat<2, 2, T, Q>>
exponential(sl2_lie_algebra<T, Q> const &A) {

  using mat2_t = glm::mat<2, 2, T, Q>;
  std::complex<T> const detA = det(A.matrix);

  std::complex<T> const a =
      std::sqrt(detA); // Below values are independent of chosen root

  std::complex<T> const sinc_a = (std::norm(a) < 1e-10) ? a : std::sin(a) / a;
  std::complex<T> const cos_a = std::cos(a);

  // exp(A) = cos(a) I + sin(a) / a A
  ImaginaryExtension<mat2_t> expA =
      ImaginaryExtension<mat2_t>{mat2_t(cos_a.real()), mat2_t(cos_a.imag())} +
      ImaginaryExtension<T>{sinc_a.real(), sinc_a.imag()} * A.matrix;

  return expA;
}

template <typename T, qualifier Q>
inline glm::vec<4, T, Q>
hermitian_matrix_to_vec(ImaginaryExtension<glm::mat<2, 2, T, Q>> const &M) {
  glm::vec<4, T, Q> result{};
  result[0] = M.real[0][1];
  result[1] = M.imag[1][0];
  result[2] = (M.real[0][0] - M.real[1][1]) / 2;
  result[3] = (M.real[0][0] + M.real[1][1]) / 2;
  return result;
}

// a simplification of hermitian_matrix_to_vec for when M is real and symmetric
template <typename T, qualifier Q>
inline glm::vec<4, T, Q>
symmetric_matrix_to_vec(glm::mat<2, 2, T, Q> const &M) {
  glm::vec<4, T, Q> result{};
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
template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q>
action_via_hermitian_matrix(ImaginaryExtension<glm::mat<2, 2, T, Q>> const &M) {

  using mat2_t = glm::mat<2, 2, T, Q>;
  ImaginaryExtension<mat2_t> const Mct{glm::transpose(M.real),
                                       -glm::transpose(M.imag)};

  // NB Ew is the identity so isn't needed
  ImaginaryExtension<mat2_t> Ex{mat2_t{0, 1, 1, 0}, mat2_t{0, 0, 0, 0}};
  ImaginaryExtension<mat2_t> Ey{mat2_t{0, 0, 0, 0}, mat2_t{0, -1, 1, 0}};
  ImaginaryExtension<mat2_t> Ez{mat2_t{1, 0, 0, -1}, mat2_t{0, 0, 0, 0}};

  glm::mat<4, 4, T, Q> result{};
  result[0] = hermitian_matrix_to_vec(M * Ex * Mct);
  result[1] = hermitian_matrix_to_vec(M * Ey * Mct);
  result[2] = hermitian_matrix_to_vec(M * Ez * Mct);
  result[3] = hermitian_matrix_to_vec(M * Mct);
  return result;
}

// a simplification of action_via_hermitian_matrix for when M is real and
// symmetric
template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q>
action_via_symmetric_matrix(glm::mat<2, 2, T, Q> const &M) {
  
  using mat2_t = glm::mat<2, 2, T, Q>;
  using mat4_t = glm::mat<4, 4, T, Q>;
  using vec4_t = glm::vec<4, T, Q>;

  mat2_t const Mct = glm::transpose(M);

  // NB Ew is the identity so isn't needed
  mat2_t Ex{0, 1, 1, 0};
  mat2_t imagEy{0, -1, 1, 0};
  mat2_t Ez{1, 0, 0, -1};

  mat4_t result{};
  result[0] = symmetric_matrix_to_vec(M * Ex * Mct);
  result[1] = vec4_t{0, (M * imagEy * Mct)[1][0], 0, 0};
  result[2] = symmetric_matrix_to_vec(M * Ez * Mct);
  result[3] = symmetric_matrix_to_vec(M * Mct);
  return result;
}

} // namespace detail

template <typename T, qualifier Q>
inline glm::mat<4, 4, T, Q> exponential(lorentz_lie_algebra_t<T, Q> const &v) {
  using mat2_t = glm::mat<2, 2, T, Q>;
  using sl2 = detail::sl2_lie_algebra<T, Q>;
  using vec3_t = glm::vec<3, T, Q>;

  vec3_t const &w = v.rotational;
  vec3_t const &b = v.boost;

  sl2 A = sl2::from(v);

  ImaginaryExtension<mat2_t> expA = exponential(A);

  return detail::action_via_hermitian_matrix(expA);
}

} // namespace hyperbolic