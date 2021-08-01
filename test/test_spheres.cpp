#include <gtest/gtest.h>

#include "spheres.h"

#include "glm_assertions.h"

using vec4_t = glm::vec4;
using vec3_t = glm::vec3;
using hangle_t = hyperbolic::hyperbolic_angle<float>;
namespace lorentz = glm::lorentz;
namespace spheres = hyperbolic::spheres;

TEST(Spheres, Contains) {

  vec4_t const p = lorentz::boost(vec3_t{0.09f, 0, 0})[3];
  vec4_t const q = lorentz::boost(vec3_t{0.11f, 0, 0})[3];
  {
    vec4_t sphere_center{0, 0, 0, 1};
    hangle_t radius{0.1f};

    EXPECT_TRUE(spheres::contains_point::sphere(sphere_center, radius, p));
    EXPECT_FALSE(spheres::contains_point::sphere(sphere_center, radius, q));
  }
  {
    vec4_t point{-1, 0, 0, 1};    
    hangle_t radius{0.1f};

    EXPECT_TRUE(spheres::contains_point::horosphere(point, radius, p));
    EXPECT_FALSE(spheres::contains_point::horosphere(point, radius, q));
  }
  {
    vec4_t plane{-1, 0, 0, 0};
    hangle_t radius{0.1f};

    EXPECT_TRUE(spheres::contains_point::pseudosphere(plane, radius, p));
    EXPECT_FALSE(spheres::contains_point::pseudosphere(plane, radius, q));
  }
}

TEST(Spheres, Distance) {
  float r = 0.1f;
  float diff = 0.05f;
  vec4_t const p = lorentz::boost(vec3_t{r - diff, 0, 0})[3];
  vec4_t const q = lorentz::boost(vec3_t{r + diff, 0, 0})[3];
  {
    vec4_t center{0, 0, 0, 1};
    hangle_t radius{r};

    EXPECT_NEAR(spheres::distance_from::sphere(center, radius, p), -diff,
                epsilon);
    EXPECT_NEAR(spheres::distance_from::sphere(center, radius, q), diff,
                epsilon);
  }
  {
    vec4_t plane{-1, 0, 0, 0};
    hangle_t radius{r};

    EXPECT_NEAR(spheres::distance_from::pseudosphere(plane, radius, p), -diff,
                epsilon);
    EXPECT_NEAR(spheres::distance_from::pseudosphere(plane, radius, q), diff,
                epsilon);
  }
  {
    vec4_t point{-1, 0, 0, 1};
    hangle_t radius{r};

    EXPECT_NEAR(spheres::distance_from::horosphere(point, radius, p), -diff,
                epsilon);
    EXPECT_NEAR(spheres::distance_from::horosphere(point, radius, q), diff,
                epsilon);
  }
}

TEST(Spheres, DistanceBetweenSpheres) {
  vec4_t sphere_center{0, 0, 0, 1};
  hangle_t radius{0.1f};

  {
    vec4_t pseudosphere_plane =
        lorentz::boost(vec3_t{0.5f, 0, 0}) * vec4_t{1, 0, 0, 0};
    hangle_t pseudosphere_radius{0.2f};

    EXPECT_NEAR(
        static_cast<float>(
            hyperbolic::spheres::distance_between::sphere_and_pseudosphere(
                sphere_center, radius, pseudosphere_plane,
                pseudosphere_radius)),
        (0.5f - 0.2f - 0.1f), epsilon);
  }
  {
    vec4_t horosphere_point =
        lorentz::boost(vec3_t{0.5f, 0, 0}) * vec4_t{1, 0, 0, 1};
    hangle_t horosphere_radius{0.2f};

    EXPECT_NEAR(
        static_cast<float>(
            hyperbolic::spheres::distance_between::sphere_and_horosphere(
                sphere_center, radius, horosphere_point, horosphere_radius)),
        (0.5f - 0.2f - 0.1f), epsilon);
  }
  {
    vec4_t second_center =
        lorentz::boost(vec3_t{0.5f, 0, 0}) * vec4_t{0, 0, 0, 1};
    hangle_t second_radius{0.2f};

    EXPECT_NEAR(
        static_cast<float>(
            hyperbolic::spheres::distance_between::spheres(
                sphere_center, radius, second_center, second_radius)),
        (0.5f - 0.2f - 0.1f), epsilon);
  }
}

TEST(Spheres, ClosestPoint) {
  namespace closest = spheres::closest_point_to;

  vec4_t const p = lorentz::boost(vec3_t{0.05f, 0, 0})[3];
  vec4_t const q = lorentz::boost(vec3_t{0.15f, 0, 0})[3];
  vec4_t const r = lorentz::boost(vec3_t{0.1f, 0, 0})[3];

  {
    SCOPED_TRACE("Sphere case");
    vec4_t center{0, 0, 0, 1};
    hangle_t radius{0.1f};

    expect_close(closest::sphere(p, center, radius), r);
    expect_close(closest::sphere(q, center, radius), r);
  }
  {
    SCOPED_TRACE("Pseudosphere case");
    vec4_t plane{-1, 0, 0, 0};
    hangle_t radius{0.1f};

    expect_close(closest::pseudosphere(p, plane, radius), r);
    expect_close(closest::pseudosphere(q, plane, radius), r);
  }
  {
    SCOPED_TRACE("Horosphere case");
    vec4_t point_at_inf{-1, 0, 0, 1};
    hangle_t radius{0.1f};

    expect_close(closest::horosphere(p, point_at_inf, radius), r);
    expect_close(closest::horosphere(q, point_at_inf, radius), r);
  }
}