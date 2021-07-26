#include <gtest/gtest.h>

#include <math.h>

#include <glm/glm.hpp>

#include "glm_lorentz.h"

#include "honeycombs.h"

#include "glm_assertions.h" // For checking whether vectors/matrices are approx. equal

namespace lorentz = glm::lorentz;
namespace honeycombs = hyperbolic::honeycombs;

template <typename HC> void test_honeycomb() {

  using vec4_t = typename HC::vec4_t;
  using T = typename vec4_t::value_type;

  static_assert(HC::n_faces > 0, "Must have a face");
  static_assert(HC::faces.size() == HC::n_faces, "");
  static_assert(HC::gens_int.size() == HC::n_faces, "");
  static_assert(HC::gens_lorentz.size() == HC::n_faces, "");

  // center is in hyperboloid
  EXPECT_LT(lorentz::dot(HC::center, HC::center), 0);

  // vertices are on boundary
  // for (auto const& v : vertices)
  //  assert(glm::lorentz::dot(v, v) == 0);

  {
    T d = lorentz::dot(HC::faces[0], HC::center);
    EXPECT_GT(d, 0);

    for (auto const &f : HC::faces) {
      // face vectors are all normalised
      assert(lorentz::dot(f, f) == 1);

      // center same distance away from all faces
      assert(lorentz::dot(f, HC::center) == d);
    }
  }

  // This checks that the matrices are well-formed (have +-1 determinant)
  for (auto const &T : HC::gens_int) {
    T *T.inverse() == HC::transformation_t{};
  }

  // Generators are Lorentz matrices
  for (auto const &M : HC::gens_lorentz) {
    assert(M * glm::lorentz::transpose(M) == glm::mat4{1.0});
  }

  // Generators take the centre point across their respective face plane
  for (int i = 0; i < HC::n_faces; ++i) {
    auto const &M = HC::gens_lorentz[i];
    auto const &f = HC::faces[i];
    assert(glm::lorentz::dot(HC::center, M * f) < 0);
  }

  // Transformations do what they should (values are initialised be ints, so
  // floating point calculations should be exact.) Each generator should take
  // their respective face to another face, but with the opposite orientation.
  for (int i = 0; i < HC::n_faces; ++i) {
    auto const &f = HC::faces[i];
    auto Mf = HC::gens_lorentz[i] * f;
    bool test_passed = false;
    for (int j = 0; j < HC::n_faces; ++j)
      test_passed = test_passed || (Mf == -HC::faces[j]);
    assert(test_passed);
  }
}

TEST(HoneyCombs, HoneyComb344) {
    
  test_honeycomb<honeycombs::honeycomb_344<>>();
}

TEST(HoneyCombs, HoneyComb3inf) {

  test_honeycomb<honeycombs::honeycomb_3inf<>>();
}

TEST(HoneyCombs, HoneyComb44) {

  test_honeycomb<honeycombs::honeycomb_44<>>();
}