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
  EXPECT_NEAR(static_cast<float>(a - b), -0.1f, epsilon);

  a -= b;
  EXPECT_NEAR(static_cast<float>(a), -0.1f, epsilon);
  EXPECT_NEAR(static_cast<float>(-a), 0.1f, epsilon);

  b += a;
  EXPECT_NEAR(static_cast<float>(b), 0.1f, epsilon);

}

TEST(Hyperbolic, Interpolation) {
  float d = 0.1f;
  hyperbolic::hyperbolic_angle a{d};

  glm::vec4 p = lorentz::normalize(glm::vec4{0.1f, 0.2f, 0.3f, 1.0f});
  glm::vec4 q = glm::vec4{1, 1, 1, 2};

  glm::vec4 r = hyperbolic::move_d_from_p_to_q(a, p, q);

  EXPECT_NEAR(lorentz::dot(r, r), -1.0f, epsilon);
  EXPECT_NEAR(hyperbolic::distance(p, r), d, 2 * epsilon);
}
