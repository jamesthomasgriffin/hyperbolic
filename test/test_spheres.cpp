#include <gtest/gtest.h>

#include "spheres.h"

#include "glm_assertions.h"

using sphere_t = hyperbolic::implicit_sphere<float, glm::highp>;
using vec4_t = glm::vec4;
using vec3_t = glm::vec3;
using hangle_t = hyperbolic::hyperbolic_angle<float>;
namespace lorentz = glm::lorentz;

TEST(Spheres, Contains) {

  vec4_t const p = lorentz::boost(vec3_t{0.09f, 0, 0})[3];
  vec4_t const q = lorentz::boost(vec3_t{0.11f, 0, 0})[3];
  {
    sphere_t const sphere{vec4_t{0, 0, 0, 1}, hangle_t{0.1f},
                          sphere_t::Type::Sphere};

    EXPECT_TRUE(sphere.contains(p));
    EXPECT_FALSE(sphere.contains(q));
  }
  {
    sphere_t const sphere{vec4_t{-1, 0, 0, 1}, hangle_t{0.1f},
                          sphere_t::Type::Horosphere};

    EXPECT_TRUE(sphere.contains(p));
    EXPECT_FALSE(sphere.contains(q));
  }
  {
    sphere_t const sphere{vec4_t{-1, 0, 0, 0}, hangle_t{0.1f},
                          sphere_t::Type::Pseudosphere};

    EXPECT_TRUE(sphere.contains(p));
    EXPECT_FALSE(sphere.contains(q));
  }
}

TEST(Spheres, Distance) {
  float radius = 0.1f;
  float diff = 0.05f;
  vec4_t const p = lorentz::boost(vec3_t{radius - diff, 0, 0})[3];
  vec4_t const q = lorentz::boost(vec3_t{radius + diff, 0, 0})[3];
  {
    sphere_t const sphere{vec4_t{0, 0, 0, 1}, hangle_t{radius},
                          sphere_t::Type::Sphere};

    EXPECT_NEAR(sphere.distance_from(p), -diff, epsilon);
    EXPECT_NEAR(sphere.distance_from(q), diff, epsilon);
    EXPECT_NEAR(sphere.distance_as_angle(p), -diff, epsilon);
    EXPECT_NEAR(sphere.distance_as_angle(q), diff, epsilon);
  }
  {
    sphere_t const sphere{vec4_t{-1, 0, 0, 1}, hangle_t{radius},
                          sphere_t::Type::Horosphere};

    EXPECT_NEAR(sphere.distance_from(p), -diff, epsilon);
    EXPECT_NEAR(sphere.distance_from(q), diff, epsilon);
    EXPECT_NEAR(sphere.distance_as_angle(p), -diff, epsilon);
    EXPECT_NEAR(sphere.distance_as_angle(q), diff, epsilon);
  }
  {
    sphere_t const sphere{vec4_t{-1, 0, 0, 0}, hangle_t{radius},
                          sphere_t::Type::Pseudosphere};

    EXPECT_NEAR(sphere.distance_from(p), -diff, epsilon);
    EXPECT_NEAR(sphere.distance_from(q), diff, epsilon);
    EXPECT_NEAR(sphere.distance_as_angle(p), -diff, epsilon);
    EXPECT_NEAR(sphere.distance_as_angle(q), diff, epsilon);


  }
}

TEST(Spheres, ClosestPoint) {
  vec4_t const p = lorentz::boost(vec3_t{0.05f, 0, 0})[3];
  vec4_t const q = lorentz::boost(vec3_t{0.15f, 0, 0})[3];
  vec4_t const r = lorentz::boost(vec3_t{0.1f, 0, 0})[3];
  {
    SCOPED_TRACE("Sphere case");
    sphere_t const sphere{vec4_t{0, 0, 0, 1}, hangle_t{0.1f},
                          sphere_t::Type::Sphere};

    expect_close(sphere.closest_point(p), r);
    expect_close(sphere.closest_point(q), r);
  }
  {
    SCOPED_TRACE("Horosphere case");
    sphere_t const sphere{vec4_t{-1, 0, 0, 1}, hangle_t{0.1f},
                          sphere_t::Type::Horosphere};

    expect_close(sphere.closest_point(p), r);
    expect_close(sphere.closest_point(q), r);
  }
  {
    SCOPED_TRACE("Pseudosphere case");
    sphere_t const sphere{vec4_t{-1, 0, 0, 0}, hangle_t{0.1f},
                          sphere_t::Type::Pseudosphere};

    vec4_t cp_to_p = sphere.closest_point(p);
    vec4_t cp_to_q = sphere.closest_point(q);
    EXPECT_NEAR(lorentz::dot(cp_to_p, cp_to_p), -1, epsilon);
    EXPECT_NEAR(lorentz::dot(cp_to_q, cp_to_q), -1, epsilon);
    expect_close(cp_to_p, r);
    expect_close(cp_to_q, r);
  }
}