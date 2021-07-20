#pragma once

#include <glm/glm.hpp>

namespace glm {
namespace traits {

template <length_t L, length_t M, typename TYPE> struct to_mat;
template <length_t L, typename TYPE> struct to_vec;

template <length_t L, length_t M, typename T, qualifier Q>
struct to_mat<L, M, vec<4, T, Q>> { using type = mat<L, M, T, Q>; };
template <length_t L, length_t M, typename T, qualifier Q>
struct to_mat<L, M, vec<3, T, Q>> { using type = mat<L, M, T, Q>; };
template <length_t L, length_t M, typename T, qualifier Q>
struct to_mat<L, M, vec<2, T, Q>> { using type = mat<L, M, T, Q>; };
template <length_t L, length_t M, typename T, qualifier Q>
struct to_mat<L, M, mat<4, 4, T, Q>> { using type = mat<L, M, T, Q>; };
template <length_t L, length_t M, typename T, qualifier Q>
struct to_mat<L, M, mat<3, 4, T, Q>> { using type = mat<L, M, T, Q>; };
template <length_t L, length_t M, typename T, qualifier Q>
struct to_mat<L, M, mat<2, 4, T, Q>> { using type = mat<L, M, T, Q>; };
template <length_t L, length_t M, typename T, qualifier Q>
struct to_mat<L, M, mat<4, 3, T, Q>> { using type = mat<L, M, T, Q>; };
template <length_t L, length_t M, typename T, qualifier Q>
struct to_mat<L, M, mat<3, 3, T, Q>> { using type = mat<L, M, T, Q>; };
template <length_t L, length_t M, typename T, qualifier Q>
struct to_mat<L, M, mat<2, 3, T, Q>> { using type = mat<L, M, T, Q>; };
template <length_t L, length_t M, typename T, qualifier Q>
struct to_mat<L, M, mat<4, 2, T, Q>> { using type = mat<L, M, T, Q>; };
template <length_t L, length_t M, typename T, qualifier Q>
struct to_mat<L, M, mat<3, 2, T, Q>> { using type = mat<L, M, T, Q>; };
template <length_t L, length_t M, typename T, qualifier Q>
struct to_mat<L, M, mat<2, 2, T, Q>> { using type = mat<L, M, T, Q>; };

template <length_t L, typename T, qualifier Q>
struct to_vec<L, vec<4, T, Q>> { using type = vec<L, T, Q>; };
template <length_t L, typename T, qualifier Q>
struct to_vec<L, vec<3, T, Q>> { using type = vec<L, T, Q>; };
template <length_t L, typename T, qualifier Q>
struct to_vec<L, vec<2, T, Q>> { using type = vec<L, T, Q>; };
template <length_t L, typename T, qualifier Q>
struct to_vec<L, mat<4, 4, T, Q>> { using type = vec<L, T, Q>; };
template <length_t L, typename T, qualifier Q>
struct to_vec<L, mat<3, 4, T, Q>> { using type = vec<L, T, Q>; };
template <length_t L, typename T, qualifier Q>
struct to_vec<L, mat<2, 4, T, Q>> { using type = vec<L, T, Q>; };
template <length_t L, typename T, qualifier Q>
struct to_vec<L, mat<4, 3, T, Q>> { using type = vec<L, T, Q>; };
template <length_t L, typename T, qualifier Q>
struct to_vec<L, mat<3, 3, T, Q>> { using type = vec<L, T, Q>; };
template <length_t L, typename T, qualifier Q>
struct to_vec<L, mat<2, 3, T, Q>> { using type = vec<L, T, Q>; };
template <length_t L, typename T, qualifier Q>
struct to_vec<L, mat<4, 2, T, Q>> { using type = vec<L, T, Q>; };
template <length_t L, typename T, qualifier Q>
struct to_vec<L, mat<3, 2, T, Q>> { using type = vec<L, T, Q>; };
template <length_t L, typename T, qualifier Q>
struct to_vec<L, mat<2, 2, T, Q>> { using type = vec<L, T, Q>; };

} // namespace traits
} // namespace glm