#include <gtest/gtest.h>

#include "glm_lorentz.h"
#include "hyperbolic.h"
#include "lorentz_lie_algebra.h"

#include "glm_assertions.h"  // For checking whether vectors/matrices are approx. equal

namespace lorentz = glm::lorentz;

TEST(Hyperbolic, Distance) {
  glm::vec3 dir{0.3, 1.1, -0.2};
  glm::mat4 B = lorentz::boost(dir);

  glm::vec4 v{0, 0, 0, 1};

  EXPECT_NEAR(hyperbolic::distance(B * v, v), glm::length(dir), epsilon);
}

TEST(Hyperbolic, Angles) {
  hyperbolic::hyperbolic_angle a{0.1f}, b{0.2f};

  EXPECT_NEAR(static_cast<float>(a + b), 0.3f, epsilon);
}
