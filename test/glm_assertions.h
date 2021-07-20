#pragma once

#include <glm/glm.hpp>

constexpr float epsilon = 1e-6f;

namespace glm {

template <length_t L, typename T, qualifier Q, typename OS>
inline OS &operator<<(OS &ostr, vec<L, T, Q> const &v) {
  ostr << '(' << v[0];
  for (int i = 1; i < L; ++i)
    ostr << ',' << v[i];
  ostr << ')';
  return ostr;
}

template <length_t L, typename T, qualifier Q, typename OS>
inline OS &operator<<(OS &ostr, mat<L, L, T, Q> const &m) {
  ostr << '[' << m[0];
  for (int i = 1; i < L; ++i)
    ostr << '\n' << m[i];
  ostr << ']';
  return ostr;
}

} // namespace glm

template <glm::length_t L, typename T, glm::qualifier Q>
void expect_close(glm::vec<L, T, Q> const &v, glm::vec<L, T, Q> const &w,
                  T threshold = static_cast<T>(epsilon)) {

  EXPECT_LT(glm::length(v - w), epsilon)
      << "First vector: " << v << "\nSecond vector: " << w;
}

template <glm::length_t R, glm::length_t C, typename T, glm::qualifier Q>
void expect_close(glm::mat<R, C, T, Q> const &v, glm::mat<R, C, T, Q> const &w,
                  T threshold = static_cast<T>(epsilon)) {
  glm::mat<R, C, T, Q> diff = v - w;
  EXPECT_LT(glm::trace(diff * glm::transpose(diff)), epsilon)
    << "First matrix:\n" << v << "\nSecond matrix:\n" << w;
}