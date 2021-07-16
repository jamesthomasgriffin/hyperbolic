#pragma once
#include <algorithm>  // For std::min, std::max
#include <glm/glm.hpp>

namespace hyperbolic {

using mat4 = glm::mat4;
using mat3 = glm::mat3;
using mat2 = glm::mat2;
using vec4 = glm::vec4;
using vec3 = glm::vec3;
using vec2 = glm::vec2;

// using mat4 = glm::mat<4,4,double>;
// using mat3 = glm::mat<3,3,double>;
// using mat2 = glm::mat<2,2,double>;
// using vec4 = glm::vec<4,double>;
// using vec3 = glm::vec<3,double>;
// using vec2 = glm::vec<2,double>;

using scalar_t = vec4::value_type;

template <typename T> inline T clamp(T x, T a, T b) {
  x = std::min(x, b);
  x = std::max(x, a);
  return x;
}

} // namespace hyperbolic