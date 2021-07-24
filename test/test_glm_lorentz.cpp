#include <gtest/gtest.h>

#include <math.h>

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

#include "glm_lorentz.h"

#include "glm_assertions.h"  // For checking whether vectors/matrices are approx. equal


TEST(Lorentz, DotProduct) {
  EXPECT_EQ(glm::lorentz::dot(glm::vec4{0, 0, 0, 1}, glm::vec4{0, 0, 0, 1}),
            -1);
  EXPECT_EQ(glm::lorentz::dot(glm::vec4{1, 2, 3, 0}, glm::vec4{3, 2, 1, 0}),
            10);

  EXPECT_EQ(glm::lorentz::dot(glm::vec3{0, 0, 1}, glm::vec3{0, 0, 1}), -1);
  EXPECT_EQ(glm::lorentz::dot(glm::vec3{2, 3, 0}, glm::vec3{2, 1, 0}), 7);

  EXPECT_EQ(glm::lorentz::dot(glm::vec2{0, 1}, glm::vec2{0, 1}), -1);
  EXPECT_EQ(glm::lorentz::dot(glm::vec2{3, 0}, glm::vec2{2, 0}), 6);

  glm::vec4 v{0.1f, 0.2f, 0.3f, 0.1f};
  glm::vec4 w{-0.7f, 0.2f, 0.3f, 0.1f};

  EXPECT_EQ(glm::dot(glm::lorentz::transpose(v), w), glm::lorentz::dot(v, w));
}

TEST(Lorentz, OuterProduct) {
  glm::vec4 const v{0.1f, 0.2f, 0.3f, 0.1f};
  glm::vec4 const w{-0.7f, 0.2f, 0.3f, 0.1f};

  glm::mat4 const op = glm::lorentz::outerProduct(v, w);

  EXPECT_EQ(glm::trace(op), glm::lorentz::dot(v, w));

  expect_close(op, glm::outerProduct(v, glm::lorentz::transpose(w)));
}

TEST(Lorentz, Trace) { EXPECT_EQ(glm::trace(glm::mat4{1.0f}), 4); }

TEST(Lorentz, Transpose) {
  EXPECT_EQ(glm::lorentz::transpose(glm::mat4{1.0f}), glm::mat4{1.0f});

  glm::mat4 M{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

  EXPECT_EQ(glm::trace(glm::lorentz::transpose(M)), glm::trace(M));

  glm::mat4 N{0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1};

  EXPECT_EQ(glm::lorentz::transpose(M * N),
            glm::lorentz::transpose(N) * glm::lorentz::transpose(M));
}

TEST(Lorentz, Normalize) {
  glm::vec4 v{1, 0, 1, 1};
  glm::vec4 w{2, 2, 0, 3};
  glm::vec4 u{2, 2, 1, 3};

  EXPECT_EQ(glm::lorentz::normalize(v), v);
  EXPECT_EQ(glm::lorentz::normalize(w), w);
  EXPECT_EQ(glm::lorentz::normalize(u), u);

  glm::vec4 p{1, 2, 3, 4};
  glm::vec4 normp = glm::lorentz::normalize(p);
  EXPECT_NEAR(glm::lorentz::dot(normp, normp), -1.0f, epsilon);

  glm::vec4 q{1, 2, 3, 1.5f};
  glm::vec4 normq = glm::lorentz::normalize(q);
  EXPECT_NEAR(glm::lorentz::dot(normq, normq), 1.0f, epsilon);
}

TEST(Lorentz, Boost) { 
  
  glm::vec4 v{0, 0, 0, 1};
  float d = 0.1f;
  float coshd = std::cosh(d), sinhd = std::sinh(d);

  // Basic test of correctness
  glm::vec4 bv = glm::lorentz::boost(v, d, 0);
  expect_close(glm::lorentz::boost(v, d, 0), glm::vec4(sinhd, 0, 0, coshd));
  expect_close(glm::lorentz::boost(glm::mat4{1.0}, glm::vec3{d, 0, 0}) * v, glm::vec4(sinhd, 0, 0, coshd));

  // Boosts in same direction produce another boost
  glm::vec3 dir{0.25f, 0.5f, 0.15};
  glm::mat4 B = glm::lorentz::boost(dir);
  expect_close(glm::lorentz::boost(3.0f * dir), B * B * B);

  // Conjugates properly by a rotation matrix
  glm::mat4 R =
      glm::rotate(glm::mat4{1.0f}, 0.3f, glm::vec3{0.1f, 0.7f, -0.7f});
  glm::vec3 Rdir{R * glm::vec4{dir, 0}};
  expect_close(glm::lorentz::boost(Rdir), R * B * glm::transpose(R));
}

TEST(Lorentz, Parabolic) { 
  glm::vec4 v{0, 0, 0, 1};
  float d = 0.25f;
  glm::mat4 A = glm::lorentz::parabolic(glm::mat4{1.0f}, d, 0, 2, 1);
  glm::mat4 B = glm::lorentz::parabolic(glm::mat4{1.0f}, d, 1, 2, 1);

  expect_close(A * B, B * A);
}