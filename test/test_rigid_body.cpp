#include <gtest/gtest.h>

#include "glm_lorentz.h"
#include "rigid_body.h"

#include "glm_assertions.h"


using rigid_body_t = hyperbolic::rigid_body_t<float>;
using velocity_t = rigid_body_t::velocity_t;
namespace lorentz = glm::lorentz;

void expect_close(velocity_t const &v, velocity_t const &w, float epsilon = 1e-6f) {
  expect_close(v.rotational, w.rotational, epsilon);
  expect_close(v.boost, w.boost, epsilon);
}

TEST(RigidBody, KineticEnergy) {
  rigid_body_t rb{};

  // A rigid body with arbitrary values
  rb.frame = glm::lorentz::boost(glm::mat4{1}, glm::vec3{0.1, -0.2, 0.3});
  rb.velocity.boost = glm::vec3{0.1f, -0.1f, 1.5f};
  rb.velocity.rotational = glm::vec3{0, 0.2f, -0.4f};

  rb.distributional_inertia = glm::vec3{0.1f, 0.01f, 0.03f};
  rb.total_mass = 1.5f;

  EXPECT_NEAR(-0.5f * dot(rb.local_momentum(), rb.velocity), rb.kinetic_energy(),
              epsilon);
  EXPECT_NEAR(-0.5f * dot(rb.momentum(), rb.global_velocity()),
              rb.kinetic_energy(), epsilon);
}

TEST(RigidBody, Inertia) {
  rigid_body_t rb{};
  velocity_t v{};

  // arbitrary values
  v.boost = glm::vec3{0.1f, -0.1f, 1.5f};
  v.rotational = glm::vec3{0, 0.2f, -0.4f};

  rb.distributional_inertia = glm::vec3{0.1f, 0.01f, 0.03f};
  rb.total_mass = 1.5f;

  expect_close(rb.apply_inverse_inertia(rb.apply_inertia(v)), v);
  expect_close(rb.apply_inertia(rb.apply_inverse_inertia(v)), v);
}

//TEST(RigidBody, ImpulseInDirection) {
//  // Create particle at p with velocity v
//  glm::vec4 const p = glm::lorentz::normalize(glm::vec4{0.1f, 0.2f, -0.3f, 1.0f});
//  glm::vec4 const w{1.0f, 0.1f, 0.2f, 0.0f};
//  glm::vec4 const v = w + lorentz::dot(w, p) * p;
//  float const mass = 1.0f;
//
//  glm::vec4 const momentum = mass * v;
//
//  velocity_t impulse = rigid_body_t::calculate_local_impulse(p, v);
//}

TEST(RigidBody, ElasticCollision) {
  rigid_body_t rb1{}, rb2{};

  // Two rigid bodies with arbitrary values
  rb1.frame = glm::lorentz::boost(glm::mat4{1}, glm::vec3{0.1, -0.2, 0.3});
  rb1.velocity.boost = glm::vec3{0.1f, -0.1f, 1.5f};
  rb1.velocity.rotational = glm::vec3{0, 0.2f, -0.4f};
  rb1.distributional_inertia = glm::vec3{0.1f, 0.01f, 0.03f};
  rb1.total_mass = 1.5f;

  rb2.frame = glm::lorentz::boost(glm::mat4{1}, glm::vec3{0.2, -0.1, 0.1});
  rb2.velocity.boost = glm::vec3{0.3f, -0.12f, 0.5f};
  rb2.velocity.rotational = glm::vec3{0.1f, 0, 0.6f};
  rb2.distributional_inertia = glm::vec3{0.01f, 0.02f, 0.3f};
  rb2.total_mass = 1.2f;

  // An impulse with arbitrary values
  // The impulse is in the frame of body 1
  velocity_t impulse{};
  impulse.rotational = glm::vec3{0.1f, -0.2f, 0.3f};
  impulse.boost = glm::vec3{0.5f, -0.1f, 0.1f};


  float factor =
      rigid_body_t::impulse_magnitude_from_elastic_collision(rb1, rb2, impulse);

  float const initial_kinetic_energy =
      rb1.kinetic_energy() + rb2.kinetic_energy();

  glm::mat4 const frame1_to_frame2 = glm::lorentz::transpose(rb2.frame) * rb1.frame;
  rb1.apply_local_impulse(factor * impulse);
  rb2.apply_local_impulse(-factor * impulse.conjugated_by(frame1_to_frame2));

  float const final_kinetic_energy =
      rb1.kinetic_energy() + rb2.kinetic_energy();

  EXPECT_NEAR(initial_kinetic_energy, final_kinetic_energy, epsilon);
}

TEST(RigidBody, Impulse) {
  rigid_body_t rb{};

  // A rigid body with arbitrary values
  rb.frame = glm::lorentz::boost(glm::mat4{1}, glm::vec3{0.1, -0.2, 0.3});
  rb.velocity.boost = glm::vec3{0.1f, -0.1f, 1.5f};
  rb.velocity.rotational = glm::vec3{0, 0.2f, -0.4f};

  rb.distributional_inertia = glm::vec3{0.1f, 0.01f, 0.03f};
  rb.total_mass = 1.5f;

  // An impulse with arbitrary values
  velocity_t impulse{};
  impulse.rotational = glm::vec3{0.1f, -0.2f, 0.3f};
  impulse.boost = glm::vec3{0.5f, -0.1f, 0.1f};

  velocity_t const initial_momentum = rb.momentum();
  rb.apply_impulse(impulse);

  velocity_t const final_momentum = rb.momentum();

  expect_close(final_momentum - initial_momentum, impulse);
}


TEST(RigidBody, FreeMotion) {
  
  float delta_t = 0.001f;

  rigid_body_t rb{};
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