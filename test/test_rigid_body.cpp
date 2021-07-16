#include <gtest/gtest.h>

#include "glm_lorentz.h"
#include "rigid_body.h"

constexpr float epsilon = 1e-6f;

TEST(RigidBody, FreeMotion) {
  
  float delta_t = 0.001f;

  hyperbolic::rigid_body_t rb{};
  rb.distributional_inertia = glm::vec3{0.1f, 0.2f, 0.3f};

  EXPECT_EQ(rb.kinetic_energy(), 0);

  rb.velocity.boost = glm::vec3{0.1f, 0.2f, 0.4f};
  rb.velocity.rotational = glm::vec3{0.2f, -0.1f, 0.3f};

  float const initial_energy = rb.kinetic_energy();
  auto const initial_momentum = rb.momentum();

  int const n_steps = 10;
  for (int i= 0; i < n_steps; ++i)
    rb.free_motion(delta_t);

  EXPECT_NEAR(rb.kinetic_energy(), initial_energy, 20 * n_steps * delta_t * delta_t);

  auto const delta_p = rb.momentum() - initial_momentum;
  float const norm_delta_p = glm::dot(delta_p.boost, delta_p.boost) +
                             glm::dot(delta_p.rotational, delta_p.rotational);
  EXPECT_LT(std::sqrt(norm_delta_p), 20 * n_steps * delta_t * delta_t);
}