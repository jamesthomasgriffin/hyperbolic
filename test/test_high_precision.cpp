#include <gtest/gtest.h>

#include <math.h>

#include <glm/glm.hpp>

#include "glm_lorentz.h"

#include "glm_assertions.h"  // For checking whether vectors/matrices are approx. equal

#include "high_precision.h"
#include "honeycombs.h"
#include "imaginary_extension.h"

namespace honeycombs = hyperbolic::honeycombs;
namespace lorentz = glm::lorentz;

using PSL2R = honeycombs::detail::PSL2R<jtg::GaussianInteger<unsigned int>>;

using polyhedron_t = honeycombs::honeycomb_344<uint64_t, float>;
using hp_frame_t = hyperbolic::high_precision_frame<polyhedron_t>;
using mat4_t = glm::mat4;

TEST(HighPrecision, NormalForm) {

  // NB this is about the limit of size of boost for floating point precision
  mat4_t const B = lorentz::boost(glm::mat4{1}, glm::vec3{1.5f, 1.0f, -2.0f});

  hp_frame_t const hpB = hp_frame_t{B, {}}.normal_form();

  EXPECT_EQ(hpB,
            hp_frame_t{hpB}.normal_form()); // normal_form should be idempotent

  expect_close(B, hpB.to_matrix()); // represents same frame

  expect_close(lorentz::transpose(hpB.frame) * hpB.frame,
               mat4_t{1}); // local frame is orthogonal
}

TEST(HighPrecision, Quotients) {

  mat4_t const B = lorentz::boost(glm::mat4{1}, glm::vec3{1.5f, 1.0f, -2.0f});
  mat4_t const C = lorentz::boost(glm::mat4{1}, glm::vec3{1.0f, -0.7f, 2.0f});

  hp_frame_t const hpB = hp_frame_t{B, {}}.normal_form();

  hp_frame_t hpBC = hp_frame_t{hpB.frame * C, hpB.transformation};

  hpBC.normal_form();

  expect_close(hyperbolic::left_quotient(hpB, hpBC), C);
}

TEST(GaussianIntegers, EuclidsAlgorithm) {

  PSL2R A{2, 1, 3, 2};

  jtg::GaussianInteger<int> a{20, 5}, b{3u, 3u};
  auto result =
      jtg::EuclidsAlgorithm<jtg::GaussianInteger<int>>::bezouts_identity(a, b);

  auto gcd = result.m * a + result.n * b;

  EXPECT_EQ(gcd, result.gcd);
}

