#include <gtest/gtest.h>

#include <math.h>

#include <glm/glm.hpp>

#include "glm_lorentz.h"

#include "glm_assertions.h"  // For checking whether vectors/matrices are approx. equal

#include "high_precision.h"
#include "honeycombs.h"
#include "imaginary_extension.h"

namespace honeycombs = hyperbolic::honeycombs;

using PSL2R = honeycombs::detail::PSL2R<jtg::GaussianInteger<unsigned int>>;

TEST(HighPrecision, NormalForm) {}

TEST(GaussianIntegers, EuclidsAlgorithm) {

  PSL2R A{2, 1, 3, 2};

  jtg::GaussianInteger<int> a{20, 5}, b{3u, 3u};
  auto result =
      jtg::EuclidsAlgorithm<jtg::GaussianInteger<int>>::bezouts_identity(a, b);
  std::cout << result.gcd << ' ' << result.m << ' ' << result.n;





}

