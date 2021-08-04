#pragma once

#include <array>

#include <glm/glm.hpp>

#include "glm_lorentz.h"

#include "hyperbolic.h"

#include "imaginary_extension.h"

namespace hyperbolic {

namespace honeycombs {

namespace detail {

template <typename T> using ImaginaryExtension = jtg::ImaginaryExtension<T>;
template <typename T> using GaussianInteger = jtg::GaussianInteger<T>;

template <typename R> struct PSL2R {

  R a{1}, b{0}, c{0}, d{1};

  constexpr PSL2R operator*(PSL2R const &y) const {
    return PSL2R{a * y.a + b * y.c, a * y.b + b * y.d, c * y.a + d * y.c,
                 c * y.b + d * y.d};
  }

#pragma warning(push)
#pragma warning(disable : 4146)
  constexpr PSL2R inverse() const { return PSL2R{d, -b, -c, a}; }
  constexpr PSL2R operator-() const { return PSL2R{-a, -b, -c, -d}; }

  constexpr bool operator==(PSL2R const &Y) const {
    return ((a == Y.a) && (b == Y.b) && (c == Y.c) && (d == Y.d)) ||
           ((a == -Y.a) && (b == -Y.b) && (c == -Y.c) && (d == -Y.d));
  }
#pragma warning(pop)

  constexpr bool operator!=(PSL2R const &Y) const { return !(*this == Y); }

  template <typename S> constexpr operator PSL2R<S>() const {
    return PSL2R<S>{static_cast<S>(a), static_cast<S>(b), static_cast<S>(c),
                    static_cast<S>(d)};
  }

  static constexpr PSL2R basic_parabolic(R b) {
    return PSL2R{static_cast<R>(1), b, static_cast<R>(0), static_cast<R>(1)};
  }
};

template <typename OS, typename R> OS &operator<<(OS &ostr, PSL2R<R> const &A) {
  ostr << '[' << '[' << A.a << ' ' << A.b << ']' << '[' << A.c << ' ' << A.d
       << ']' << ']';
  return ostr;
}

template <typename INT, typename T>
inline constexpr ImaginaryExtension<glm::mat<2, 2, T>>
to_complex_mat2(PSL2R<GaussianInteger<INT>> const &J, T) {

  // To correctly convert to a floating point matrix we first convert to a
  // signed integer type
  using int_type = typename std::make_signed<INT>::type;
  using mat2_t = glm::mat<2, 2, T>;

  // Note, our PSL2R class is row major, default glm behaviour is colomn major
  mat2_t real_part{
      static_cast<int_type>(J.a.real), static_cast<int_type>(J.c.real),
      static_cast<int_type>(J.b.real), static_cast<int_type>(J.d.real)};
  mat2_t imag_part{
      static_cast<int_type>(J.a.imag), static_cast<int_type>(J.c.imag),
      static_cast<int_type>(J.b.imag), static_cast<int_type>(J.d.imag)};

  return ImaginaryExtension<mat2_t>{real_part, imag_part};
}

template <typename INT, typename T>
inline constexpr glm::mat<4, 4, T> to_mat4(PSL2R<GaussianInteger<INT>> const &J,
                                           T tag) {

  ImaginaryExtension<glm::mat<2, 2, T>> M = to_complex_mat2(J, tag);
  return hyperbolic::detail::action_via_hermitian_matrix(M);
}

template <typename INT, typename T>
inline constexpr glm::mat<4, 4, T> to_mat4(PSL2R<INT> const &J, T) {

  // To correctly convert to a floating point matrix we first convert to a
  // signed integer type
  using int_type = typename std::make_signed<INT>::type;

  // PSL2R is row major, glm is column major
  glm::mat<2, 2, T> M{static_cast<int_type>(J.a), static_cast<int_type>(J.c),
                      static_cast<int_type>(J.b), static_cast<int_type>(J.d)};
  return hyperbolic::detail::action_via_symmetric_matrix(M);
}

} // namespace detail

// cell an ideal octahedron
template <typename INT = uint64_t, typename T = float> struct honeycomb_344 {

  using vec4_t = glm::vec<4, T>;
  using mat4_t = glm::mat<4, 4, T>;
  using scalar_t = T;

  static constexpr int n_vertices = 6;
  static constexpr int n_edges = 12;
  static constexpr int n_faces = 8;

  using transformation_t = detail::PSL2R<detail::GaussianInteger<INT>>;

  static constexpr vec4_t center{1, 1, 0, 2};

  static constexpr std::array<vec4_t, n_faces> faces{
      vec4_t{0, 1, 0, 0},    vec4_t{1, 0, 0, 0},   vec4_t{-1, 0, -1, -1},
      vec4_t{0, -1, -1, -1}, vec4_t{-1, 0, 1, -1}, vec4_t{0, -1, 1, -1},
      vec4_t{-2, -1, 0, -2}, vec4_t{-1, -2, 0, -2}};

  // static const std::array<vec4_t, n_vertices> vertices;

  static const std::array<transformation_t, n_faces> gens_int;
  static const std::array<mat4_t, n_faces> gens_lorentz;

  static void initialise_static_constants();
};

// cell an ideal triangle extruded to infinity
template <typename INT = uint64_t, typename T = float> struct honeycomb_3inf {

  using vec4_t = glm::vec<4, T>;
  using mat4_t = glm::mat<4, 4, T>;
  using scalar_t = T;

  using transformation_t = detail::PSL2R<INT>;

  static constexpr int n_faces = 3;

  static constexpr vec4_t center{static_cast<scalar_t>(0.5), 0, 0, 1};

  static constexpr std::array<vec4_t, n_faces> faces{
      vec4_t{1, 0, 0, 0}, vec4_t{-1, 0, -1, -1}, vec4_t{-1, 0, 1, -1}};

  // static const std::array<vec4_t, n_vertices> vertices;

  static const std::array<transformation_t, n_faces> gens_int;
  static const std::array<mat4_t, n_faces> gens_lorentz;
};

// cell a square pyramid with single vertex at infinity (0, 0, 1, 1)
template <typename INT = uint64_t, typename T = float> struct honeycomb_44 {

  using vec4_t = glm::vec<4, T>;
  using mat4_t = glm::mat<4, 4, T>;

  using transformation_t = detail::PSL2R<detail::GaussianInteger<INT>>;

  static constexpr int n_faces = 4;

  static constexpr vec4_t center{1, 1, 0, 2};

  static constexpr std::array<vec4_t, n_faces> faces{
      vec4_t{0, 1, 0, 0}, vec4_t{1, 0, 0, 0}, vec4_t{-1, 0, -1, -1},
      vec4_t{0, -1, -1, -1}};

  // static const std::array<vec4_t, n_vertices> vertices;

  static const std::array<transformation_t, n_faces> gens_int;
  static const std::array<mat4_t, n_faces> gens_lorentz;
};

#pragma warning(push)
#pragma warning(disable : 4146)
template <typename INT, typename T>
const std::array<typename honeycomb_344<INT, T>::transformation_t, 4>
    parabolicZi{
        honeycomb_344<INT, T>::transformation_t::basic_parabolic({0, 1}),
        honeycomb_344<INT, T>::transformation_t::basic_parabolic({1, 0}),
        honeycomb_344<INT, T>::transformation_t::basic_parabolic({-1ull, 0}),
        honeycomb_344<INT, T>::transformation_t::basic_parabolic({0, -1ull})};

template <typename INT, typename T>
const typename honeycomb_344<INT, T>::transformation_t involution{
    detail::GaussianInteger<uint64_t>{1u, 0u},
    detail::GaussianInteger<uint64_t>{0u, 0u},
    detail::GaussianInteger<uint64_t>{1u, -1ull},
    detail::GaussianInteger<uint64_t>{-1ull, 0u}};

template <typename INT, typename T>
const std::array<typename honeycomb_344<INT, T>::transformation_t,
                 honeycomb_344<INT, T>::n_faces>
    honeycomb_344<INT, T>::gens_int{
        parabolicZi<INT, T>[0],
        parabolicZi<INT, T>[1],
        parabolicZi<INT, T>[2],
        parabolicZi<INT, T>[3],
        involution<INT, T> *parabolicZi<INT, T>[1] * involution<INT, T>,
        involution<INT, T> *parabolicZi<INT, T>[0] * involution<INT, T>,
        involution<INT, T> *parabolicZi<INT, T>[3] * involution<INT, T>,
        involution<INT, T> *parabolicZi<INT, T>[2] * involution<INT, T>};

template <typename INT, typename T>
const std::array<glm::mat<4, 4, T>, honeycomb_344<INT, T>::n_faces>
    honeycomb_344<INT, T>::gens_lorentz{
        to_mat4(honeycomb_344::gens_int[0], T{}),
        to_mat4(honeycomb_344::gens_int[1], T{}),
        to_mat4(honeycomb_344::gens_int[2], T{}),
        to_mat4(honeycomb_344::gens_int[3], T{}),
        to_mat4(honeycomb_344::gens_int[4], T{}),
        to_mat4(honeycomb_344::gens_int[5], T{}),
        to_mat4(honeycomb_344::gens_int[6], T{}),
        to_mat4(honeycomb_344::gens_int[7], T{})};


template <typename INT, typename T>
const std::array<typename honeycomb_44<INT, T>::transformation_t,
                 honeycomb_44<INT, T>::n_faces>
    honeycomb_44<INT, T>::gens_int{
        honeycomb_44::transformation_t::basic_parabolic({0, 1}),
        honeycomb_44::transformation_t::basic_parabolic({1, 0}),
        honeycomb_44::transformation_t::basic_parabolic({-1ull, 0}),
        honeycomb_44::transformation_t::basic_parabolic({0, -1ull})};

template <typename INT, typename T>
const std::array<glm::mat<4, 4, T>, honeycomb_44<INT, T>::n_faces>
    honeycomb_44<INT, T>::gens_lorentz{to_mat4(honeycomb_44::gens_int[0], T{}),
                                       to_mat4(honeycomb_44::gens_int[1], T{}),
                                       to_mat4(honeycomb_44::gens_int[2], T{}),
                                       to_mat4(honeycomb_44::gens_int[3], T{})};


template <typename INT, typename T>
const std::array<typename honeycomb_3inf<INT, T>::transformation_t, honeycomb_3inf<INT, T>::n_faces>
   honeycomb_3inf<INT, T>::gens_int{
       honeycomb_3inf::transformation_t::basic_parabolic(1),
       honeycomb_3inf::transformation_t::basic_parabolic(-1ull),
       {1, 0, -1ull, 1}};

template <typename INT, typename T>
const std::array<glm::mat<4, 4, T>, honeycomb_3inf<INT, T>::n_faces> honeycomb_3inf<INT, T>::gens_lorentz{
   to_mat4(honeycomb_3inf::gens_int[0], T{}),
   to_mat4(honeycomb_3inf::gens_int[1], T{}),
   to_mat4(honeycomb_3inf::gens_int[2], T{})};

#pragma warning(pop)

} // namespace honeycombs
} // namespace hyperbolic
