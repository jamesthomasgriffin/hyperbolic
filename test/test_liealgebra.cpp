#include <gtest/gtest.h>

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

#include "glm_lorentz.h"
#include "hyperbolic.h"
#include "lorentz_lie_algebra.h"

#include "glm_assertions.h"  // For checking whether vectors/matrices are approx. equal

namespace lorentz = glm::lorentz;



TEST(LieAlgebra, Exponential) {

  glm::vec3 dir{0.2, -0.1, -0.2};

  {
    SCOPED_TRACE("Boost");
    hyperbolic::lorentz_lie_algebra_t boost{{0, 0, 0}, dir};
    expect_close(hyperbolic::exponential(boost), lorentz::boost(dir));
  }


  {
    SCOPED_TRACE("Rotation");
    hyperbolic::lorentz_lie_algebra_t rot{dir, {0, 0, 0}};
    expect_close(hyperbolic::exponential(rot),
               glm::rotate(glm::mat4{1}, glm::length(dir), dir));

    expect_close(rot.to_matrix() * glm::vec4{dir, 0}, glm::vec4{0});
  }
}


TEST(LieAlgebra, MatrixRepresentation) {

  {
    SCOPED_TRACE("Expected matrix form");
    glm::vec3 w{1, 2, 3};
    glm::vec3 b{4, 5, 6};
    glm::mat4 W = hyperbolic::lorentz_lie_algebra_t{w, b}.to_matrix();
    glm::mat4 Wexpected{0, 3, -2, 4, -3, 0, 1, 5, 2, -1, 0, 6, 4, 5, 6, 0};
    expect_close(W, Wexpected);
  }
}

TEST(LieAlgebra, Bracket) {
  hyperbolic::lorentz_lie_algebra_t g{{0.3, 1.1, -0.2}, {0.2, -0.1, -0.2}};
  hyperbolic::lorentz_lie_algebra_t h{{0.1, 0.1, 0.2}, {0.3, 0.7, -0.1}};

  glm::mat4 Mg = g.to_matrix();
  glm::mat4 Mh = h.to_matrix();

  glm::mat4 Mgh = hyperbolic::bracket(g, h).to_matrix();

  expect_close(Mgh, Mg * Mh - Mh * Mg);
}

TEST(LieAlgebra, Action) {

  glm::vec3 dir{0.2, -0.1, -0.2};
  glm::vec4 dir4{dir, 0};

  {
    SCOPED_TRACE("Boost");
    hyperbolic::lorentz_lie_algebra_t boost{{0, 0, 0}, dir};
    expect_close(boost * dir4, glm::vec4{0, 0, 0, glm::dot(dir, dir)});
  }

  {
    SCOPED_TRACE("Rotation");
    hyperbolic::lorentz_lie_algebra_t rot{dir, {0, 0, 0}};
    expect_close(rot * dir4, {});
  }
}