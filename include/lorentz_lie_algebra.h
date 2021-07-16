#pragma once

#include <complex>

#include <glm/glm.hpp>

#include "vector_types.h"
#include "glm_lorentz.h"

struct lorentz_lie_algebra_t {
  using scalar_t = typename vec3::value_type;
  vec3 rotational{0.0f}; // Anti-symmetric part of Lie algebra
  vec3 boost{0.0f};      // Symmetric part of Lie algebra

  lorentz_lie_algebra_t operator+(lorentz_lie_algebra_t const &b) const;
  lorentz_lie_algebra_t operator-(lorentz_lie_algebra_t const &b) const;
  lorentz_lie_algebra_t operator*(scalar_t c) const;
  vec4 operator*(vec4 const &v) const;

  lorentz_lie_algebra_t &operator+=(lorentz_lie_algebra_t const &b);
  lorentz_lie_algebra_t &operator-=(lorentz_lie_algebra_t const &b);
  lorentz_lie_algebra_t &operator*=(scalar_t c);

  lorentz_lie_algebra_t conjugated_by(mat4 const &P) const;
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

inline lorentz_lie_algebra_t
lorentz_lie_algebra_t::operator+(lorentz_lie_algebra_t const &b) const {
  return {rotational + b.rotational, boost + b.boost};
}

inline lorentz_lie_algebra_t
lorentz_lie_algebra_t::operator-(lorentz_lie_algebra_t const &b) const {
  return {rotational - b.rotational, boost - b.boost};
}

inline lorentz_lie_algebra_t lorentz_lie_algebra_t::operator*(scalar_t c) const {
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

inline mat4 to_matrix(lorentz_lie_algebra_t const &v) {
}

inline ImaginaryExtension<mat2> to_sl2(lorentz_lie_algebra_t const &v) {
  vec3 const &e = v.rotational;
  vec3 const &f = v.boost;

  return ImaginaryExtension<mat2>{0.5f * mat2{f[2], f[0] - e[1], f[0] + e[1], -f[2]},
                                  0.5f * mat2{-e[2], e[0] - f[1], e[0] + f[1], e[2]}};
}

inline mat4 exponential(lorentz_lie_algebra_t const &v) {
  vec3 const &e = v.rotational;
  vec3 const &f = v.boost;

  using T = mat4::value_type;

  ImaginaryExtension<mat2> A = to_sl2(v);

  std::complex<T> const detA{glm::dot(e, e) / 4 - glm::dot(f, f) / 4,
                             glm::dot(e, f) / 2};

  std::complex<T> const a =
      std::sqrt(detA); // Below values are independent of chosen root

  std::complex<T> const sinc_a = (std::norm(a) < 1e-12) ? a : sin(a) / a;
  std::complex<T> const cos_a = cos(a);


  // exp(A) = cos(a) I + sin(a) / a A
  ImaginaryExtension<mat2> expA =
      ImaginaryExtension<mat2>{mat2(cos_a.real()), mat2(cos_a.imag())} +
      ImaginaryExtension<T>{sinc_a.real(), sinc_a.imag()} * A;

  return glm::lorentz::detail::action_via_hermitian_matrix(expA);
}
