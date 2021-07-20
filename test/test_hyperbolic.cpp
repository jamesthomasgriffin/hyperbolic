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


TEST(Hyperbolic, Details) { 
  {
    glm::vec3 const rotation{0.25f, 0.5f, 1.5f};
    float const magnitude = glm::length(rotation);
    ASSERT_LT(magnitude, 3.1415f);

    glm::vec3 const axis = glm::normalize(rotation);
    glm::mat4 const R = glm::rotate(glm::mat4{1}, magnitude, axis);

    glm::mat3 const log_R =
        hyperbolic::detail::log_special_orthogonal_matrix(R);

    ASSERT_EQ(log_R, -glm::transpose(log_R));

    glm::vec3 const r{log_R[1][2], log_R[2][0], log_R[0][1]};
    expect_close(r, rotation);
  }

  {
    // A different path is used for rotations close to pi
    glm::mat3 R{1, 0, 0, 0, -1, 0, 0, 0, -1};

    glm::mat3 const log_R =
        hyperbolic::detail::log_special_orthogonal_matrix(R);
    ASSERT_EQ(log_R, -glm::transpose(log_R));

    glm::vec3 const r{log_R[1][2], log_R[2][0], log_R[0][1]};

    expect_close(glm::abs(r), glm::vec3{3.1415926f, 0, 0});
  }
}