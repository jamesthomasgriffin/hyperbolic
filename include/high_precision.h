#pragma once

#include <array>
#include <optional>
#include <stdint.h>

#include "glm/glm.hpp"
#include "glm_lorentz.h"
#include "hyperbolic.h"

#include "imaginary_extension.h"
#include "rigid_body.h"

namespace hyperbolic {

template <class POLYHEDRON> struct high_precision_point;
template <class POLYHEDRON> struct high_precision_frame;

template <typename POLYHEDRON>
std::optional<size_t> closest_face(typename POLYHEDRON::vec4_t const &p);

template <typename POLYHEDRON>
typename POLYHEDRON::transformation_t
transform_to_domain(typename POLYHEDRON::vec4_t const &p);

template <class POLYHEDRON> struct high_precision_point {

  using vec4_t = typename POLYHEDRON::vec4_t;
  using value_type = typename vec4_t::value_type;
  using mat4_t = typename POLYHEDRON::mat4_t;

  using polyhedron_t = POLYHEDRON;
  using transformation_t = typename POLYHEDRON::transformation_t;

  vec4_t point{1.0f};      // low precision floating point part
  transformation_t transformation{}; // high precision transformation

  mat4_t to_vector() const {
    return to_mat4(transformation, value_type{1}) * point;
  }

  void normal_form() {

    transformation_t J = transform_to_domain<polyhedron_t>(point);
    transformation = transformation * J.inverse();
    point = to_mat4(J) * point;
  }
};

template <class Trans, typename T>
inline glm::mat<4, 4, T> left_quotient(glm::mat<4, 4, T> const &quot_frame,
                                      Trans const &quot_trans,
                                      Trans const &r_trans) {
  return lorentz::transpose(quot_frame) *
         to_mat4(quot_trans.inverse() * r_trans, T{});
}

template <class Trans, typename T>
inline glm::mat<4, 4, T>
left_quotient(glm::mat<4, 4, T> const &quot_frame, Trans const &quot_trans,
              glm::mat<4, 4, T> const &r_frame, Trans const &r_trans) {
  return left_quotient(quot_frame, quot_trans, r_trans) * r_frame;
}

template <class Trans, typename T>
inline glm::vec<4, T> left_quotient(glm::mat<4, 4, T> const &quot_frame, Trans const &quot_trans,
                          glm::vec<4, T> const &r_vec, Trans const &r_trans) {
  return left_quotient(quot_frame, quot_trans, r_trans) * r_vec;
}

template <class POLYHEDRON> struct high_precision_frame {
  
  using vec4_t = typename POLYHEDRON::vec4_t;
  using mat4_t = typename POLYHEDRON::mat4_t;
  using value_type = typename mat4_t::value_type;

  using polyhedron_t = POLYHEDRON;
  using position_t = high_precision_point<POLYHEDRON>;
  using transformation_t = typename POLYHEDRON::transformation_t;

  mat4_t frame{1.0f};                  // low precision floating point part
  transformation_t transformation{}; // high precision transformation

  bool operator==(high_precision_frame const& other) const {
    return (frame == other.frame) && (transformation == other.transformation);
  }

  position_t position() const { return position_t{frame[3], transformation}; }

  mat4_t to_matrix() const {
    return to_mat4(transformation, value_type{}) * frame;
  }

  high_precision_frame& normal_form() {
    transformation_t J = transform_to_domain<polyhedron_t>(frame[3]);
    transformation = transformation * J.inverse();
    frame = to_mat4(J, value_type{}) * frame;
    return *this;
  }
};

template <class POLYHEDRON>
inline typename POLYHEDRON::transformation_t normal_form(typename POLYHEDRON::transformation_t &T, typename POLYHEDRON::mat4_t &frame) {
  typename POLYHEDRON::transformation_t J = transform_to_domain<POLYHEDRON>(frame[3]);
  T = T * J.inverse();
  frame = to_mat4(J, typename POLYHEDRON::scalar_t{}) * frame;
  return J;
}

template <class POLYHEDRON>
inline glm::mat<4, 4, typename POLYHEDRON::scalar_t> left_quotient(high_precision_frame<POLYHEDRON> const &quot,
                          high_precision_frame<POLYHEDRON> const &r_fr) {
  //using Trans_t = typename high_precision_frame<POLYHEDRON>::transformation_t;
  return left_quotient(quot.frame, quot.transformation, r_fr.frame,
                                r_fr.transformation);
}

template <class POLYHEDRON>
inline glm::vec<4, typename POLYHEDRON::scalar_t>
left_quotient(high_precision_frame<POLYHEDRON> const &quot,
                          high_precision_point<POLYHEDRON> const &p) {
  //using Trans_t = typename high_precision_frame<POLYHEDRON>::transformation_t;
  return left_quotient(quot.frame, quot.transformation, p.point,
                                p.transformation);
}

template <class POLYHEDRON> struct high_precision_rigid_body {

  using mat4_t = typename POLYHEDRON::mat4_t;
  using value_type = typename mat4_t::value_type;

  using polyhedron_t = POLYHEDRON;
  using position_t = high_precision_point<POLYHEDRON>;
  using frame_t = high_precision_frame<POLYHEDRON>;
  using transformation_t = typename POLYHEDRON::transformation_t;

  rigid_body_t<value_type> body{};   // low precision floating point part
  transformation_t transformation{}; // high precision transformation

  position_t position() const {
    return position_t{body.frame[3], transformation};
  }

  frame_t frame() const { return frame_t{body.frame, transformation}; }

  void normal_form() {
    transformation_t J = transform_to_domain<polyhedron_t>(body.frame[3]);
    transformation = transformation * J.inverse();
    body.frame = to_mat4(J, value_type{}) * body.frame;
  }
};

template <typename POLYHEDRON>
std::optional<size_t> closest_face(typename POLYHEDRON::vec4_t const &p) {

  // In a regular polyhedron the closed line segment joining its center to a
  // point p passes through face F if the orthogonal distance between p and F is
  // non-negative and greatest amongst all other faces.

  using vec4_t = typename POLYHEDRON::vec4_t;
  using scalar_t = typename vec4_t::value_type;

  std::optional<size_t> result{std::nullopt};

  // If all distances are negative then p is in the polyhedron, so
  // initialise with 0.
  scalar_t furthest_dist = 0;

  for (size_t i = 0; i < POLYHEDRON::n_faces; ++i) {
    scalar_t dist = -glm::lorentz::dot(POLYHEDRON::faces[i], p);
    if (dist > furthest_dist) {
      result = i;
      furthest_dist = dist;
    }
  }
  return result;
}

template <typename POLYHEDRON>
typename POLYHEDRON::transformation_t transform_to_domain(typename POLYHEDRON::vec4_t const &p) {

  using vec4_t = typename POLYHEDRON::vec4_t;
  typename POLYHEDRON::transformation_t result{};
  vec4_t current_pos{p};

  bool done = false;
  int count = 0;

  assert(glm::lorentz::dot(p, p) < 0);

  while (!done) {
    std::optional<size_t> cl_face = closest_face<POLYHEDRON>(current_pos);
    if (cl_face.has_value()) {
      size_t face_ix = *cl_face;
      current_pos = POLYHEDRON::gens_lorentz[face_ix] * current_pos;
      result = POLYHEDRON::gens_int[face_ix] * result;
    } else
      done = true;

    assert(++count < 1000); // Fail here if something isn't working, should be
                            // replaced by an exception once code is mature
  }
  return result;
}

// UNTESTED
// This is a higher precision variant intended for comparing between two types of high precision point
template <typename PH1, typename PH2>
typename PH1::transformation_t
transform_to_domain(high_precision_point<PH2> const &p) {

  typename PH1::transformation_t result{};
  high_precision_point<PH2> current_pos{p};

  bool done = false;
  int count = 0;

  while (!done) {
    std::optional<size_t> cl_face = closest_face<PH1>(current_pos.to_vector());
    if (cl_face.has_value()) {
      size_t face_ix = *cl_face;
      current_pos.transformation =
          PH1::gens_int[face_ix] * current_pos.transformation;
    } else
      done = true;

    assert(++count < 1000); // Fail here if something isn't working, might be
                            // replaced by an exception once code is mature
  }
  return current_pos.transformation * p.transformation.inverse();
}

} // namespace hyperbolic