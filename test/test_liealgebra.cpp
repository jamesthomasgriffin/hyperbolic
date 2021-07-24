#include <gtest/gtest.h>

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

#include "glm_lorentz.h"
#include "hyperbolic.h"
#include "lorentz_lie_algebra.h"

#include "glm_assertions.h"  // For checking whether vectors/matrices are approx. equal

namespace lorentz = glm::lorentz;
using lie_algebra_t = hyperbolic::lorentz_lie_algebra_t<float>;


TEST(LieAlgebra, Exponential) {

  glm::vec3 dir{0.2, -0.1, -0.2};

  {
    SCOPED_TRACE("Boost");
    lie_algebra_t boost{{0, 0, 0}, dir};
    expect_close(hyperbolic::exponential(boost), lorentz::boost(dir));
  }


  {
    SCOPED_TRACE("Rotation");
    lie_algebra_t rot{dir, {0, 0, 0}};
    expect_close(hyperbolic::exponential(rot),
               glm::rotate(glm::mat4{1}, glm::length(dir), dir));

    expect_close(rot.to_matrix() * glm::vec4{dir, 0}, glm::vec4{0});
  }
}


TEST(LieAlgebra, MatrixRepresentation) {

  {
    // Explicit check
    glm::vec3 w{1, 2, 3};
    glm::vec3 b{4, 5, 6};
    glm::mat4 W = lie_algebra_t{w, b}.to_matrix();
    glm::mat4 Wexpected{0, 3, -2, 4, -3, 0, 1, 5, 2, -1, 0, 6, 4, 5, 6, 0};
    EXPECT_EQ(W, Wexpected);
  }

  glm::vec3 const rotation{0.1f, 0.2f, 0.3f};
  glm::vec3 const boost{1, 2, 3};
  lie_algebra_t const g{rotation, boost};

  glm::mat4 const G = g.to_matrix();
  glm::mat3 const R = glm::mat3{G};

  ASSERT_EQ(glm::vec3{G[3]}, boost);
  ASSERT_EQ(R, -glm::transpose(R));

  glm::vec3 const v{4, 5, 6};
  ASSERT_EQ(glm::cross(rotation, v), R * v);

  ASSERT_EQ(lie_algebra_t::from_matrix(G), g);
}

TEST(LieAlgebra, Bracket) {
  lie_algebra_t g{{0.3, 1.1, -0.2}, {0.2, -0.1, -0.2}};
  lie_algebra_t h{{0.1, 0.1, 0.2}, {0.3, 0.7, -0.1}};

  glm::mat4 Mg = g.to_matrix();
  glm::mat4 Mh = h.to_matrix();

  glm::mat4 Mgh = hyperbolic::bracket(g, h).to_matrix();

  expect_close(Mgh, Mg * Mh - Mh * Mg);
}

TEST(LieAlgebra, Wedge) {

  glm::vec4 const p{0.3, 1.1, -0.2, 1.5};
  glm::vec4 const n{0.3, 0.1, -0.2, 0.5};
  lie_algebra_t const g = hyperbolic::wedge(p, n);

  glm::mat4 const g_mat = g.to_matrix();

  glm::mat4 const wedge =
      glm::lorentz::outerProduct(p, n) - glm::lorentz::outerProduct(n, p);

  expect_close(g_mat, wedge);
}

TEST(LieAlgebra, Action) {

  glm::vec3 dir{0.2, -0.1, -0.2};
  glm::vec4 dir4{dir, 0};

  {
    SCOPED_TRACE("Boost");
    lie_algebra_t boost{{0, 0, 0}, dir};
    expect_close(boost * dir4, glm::vec4{0, 0, 0, glm::dot(dir, dir)});
  }

  {
    SCOPED_TRACE("Rotation");
    lie_algebra_t rot{dir, {0, 0, 0}};
    expect_close(rot * dir4, {});
  }
}