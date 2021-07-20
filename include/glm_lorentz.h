#pragma once

#include "glm/glm.hpp"

namespace glm {

template <length_t R, length_t C, typename T, qualifier Q>
GLM_FUNC_DECL GLM_CONSTEXPR T trace(mat<R, C, T, Q> const &M);

namespace lorentz {

// Declarations

template <length_t L, typename T, qualifier Q>
GLM_FUNC_DECL GLM_CONSTEXPR T dot(vec<L, T, Q> const &x, vec<L, T, Q> const &y);

template <length_t L, typename T, qualifier Q>
GLM_FUNC_DECL vec<L, T, Q> normalize(vec<L, T, Q> const &x);

template <length_t L, typename T, qualifier Q>
GLM_FUNC_DECL vec<L, T, Q> normalize(vec<L, T, Q> const &x, int expected_sign);

template <length_t L, typename T, qualifier Q>
GLM_FUNC_DECL GLM_CONSTEXPR vec<L, T, Q>
cross(vec<L, T, Q> const &x, vec<L, T, Q> const &y, vec<L, T, Q> const &z);

template <length_t L, typename T, qualifier Q>
GLM_FUNC_DECL GLM_CONSTEXPR mat<L, L, T, Q>
complete_orthogonal_basis(vec<L, T, Q> const &p, vec<L, T, Q> const &d);

template <length_t L, typename T, qualifier Q>
GLM_FUNC_DECL mat<L, L, T, Q> transpose(mat<L, L, T, Q> const &m);

template <length_t DA, length_t DB, typename T, qualifier Q>
GLM_FUNC_DECL typename glm::detail::outerProduct_trait<DA, DB, T, Q>::type
outerProduct(vec<DA, T, Q> const &c, vec<DB, T, Q> const &r);

template <length_t C, length_t R, typename T, qualifier Q>
GLM_FUNC_DECL GLM_CONSTEXPR mat<C, R, T, Q> &
gram_schmidt_in_place(mat<C, R, T, Q> &M);

template <length_t C, length_t R, typename T, qualifier Q>
GLM_FUNC_DECL GLM_CONSTEXPR mat<C, R, T, Q>
gram_schmidt(mat<C, R, T, Q> const &M);

template <length_t C, length_t R, typename T, qualifier Q>
GLM_FUNC_DECL GLM_CONSTEXPR mat<C, R, T, Q> &
gram_schmidt_normalised_in_place(mat<C, R, T, Q> &M);

template <length_t C, length_t R, typename T, qualifier Q>
GLM_FUNC_DECL GLM_CONSTEXPR mat<C, R, T, Q>
gram_schmidt_normalised(mat<C, R, T, Q> const &M);

template <length_t L, typename T, qualifier Q>
GLM_FUNC_DECL mat<L + 1, L + 1, T, Q> boost(vec<L, T, Q> const &dir);

template <length_t L, length_t K, typename T, qualifier Q>
GLM_FUNC_DECL mat<L, L, T, Q> boost(mat<L, L, T, Q> const &M,
                                    vec<K, T, Q> const &dir);

template <length_t L, typename T, qualifier Q>
GLM_FUNC_DECL vec<L, T, Q> boost(vec<L, T, Q> const &v, T dist, size_t axis);

template <length_t L, typename T, qualifier Q>
GLM_FUNC_DECL mat<L, L, T, Q> parabolic(mat<L, L, T, Q> const &M, T dist,
                                        int axis1, int axis2, int sgn);

template <length_t L, typename T, qualifier Q>
GLM_FUNC_DECL vec<L, T, Q> parabolic(vec<L, T, Q> const &v, T dist, int axis1,
                                     int axis2, int sgn);

namespace detail {

template <typename V, typename T, bool Aligned> struct compute_dot {};

template <typename V, typename T, bool Aligned> struct compute_cross {};

template <typename V, typename T, bool Aligned>
struct compute_orthogonal_basis {};

template <typename M, typename T, bool Aligned> struct compute_transpose {};

template <typename M, typename T, bool Aligned> struct compute_boost {};

template <typename M, typename T, bool Aligned> struct compute_parabolic {};

template <typename M, typename T, bool Aligned> struct perform_gram_schmidt {};

template <typename M, typename T, bool Aligned>
struct perform_gram_schmidt_normalised {};

template <typename T, qualifier Q, bool Aligned>
struct compute_dot<vec<2, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static T call(vec<2, T, Q> const &a,
                                                 vec<2, T, Q> const &b) {
    vec<2, T, Q> tmp(a * b);
    return tmp.x - tmp.y;
  }
};

template <typename T, qualifier Q, bool Aligned>
struct compute_dot<vec<3, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static T call(vec<3, T, Q> const &a,
                                                 vec<3, T, Q> const &b) {
    vec<3, T, Q> tmp(a * b);
    return tmp.x + tmp.y - tmp.z;
  }
};

template <typename T, qualifier Q, bool Aligned>
struct compute_dot<vec<4, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static T call(vec<4, T, Q> const &a,
                                                 vec<4, T, Q> const &b) {
    vec<4, T, Q> tmp(a * b);
    return (tmp.x + tmp.y) + (tmp.z - tmp.w);
  }
};

template <length_t L, typename T, qualifier Q, bool Aligned>
struct compute_normalize {
  GLM_FUNC_QUALIFIER static vec<L, T, Q> call(vec<L, T, Q> const &v) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                      "'normalize' accepts only floating-point inputs");

    T vv = lorentz::dot(v, v);
    return (vv == 0) ? v : v * inversesqrt(abs(vv));
  }
  GLM_FUNC_QUALIFIER static vec<L, T, Q> call(vec<L, T, Q> const &v,
                                              int expected_sign) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                      "'normalize' accepts only floating-point inputs");
    T vv = lorentz::dot(v, v);
    assert(expected_sign == (vv > 0) - (vv < 0));
    return (vv == 0) ? v : v * inversesqrt(abs(vv));
  }
};

template <typename T, qualifier Q, bool Aligned>
struct compute_cross<vec<4, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static vec<4, T, Q>
  call(vec<4, T, Q> const &a, vec<4, T, Q> const &b, vec<4, T, Q> const &c) {
    std::array<vec<4, T, Q>, 9> const tmp{vec<4, T, Q>{a[1], a[0], a[0], a[0]},
                                          vec<4, T, Q>{a[2], a[2], a[1], a[1]},
                                          vec<4, T, Q>{a[3], a[3], a[3], a[2]},
                                          vec<4, T, Q>{b[1], b[0], b[0], b[0]},
                                          vec<4, T, Q>{b[2], b[2], b[1], b[1]},
                                          vec<4, T, Q>{b[3], b[3], b[3], b[2]},
                                          vec<4, T, Q>{c[1], c[0], c[0], c[0]},
                                          vec<4, T, Q>{c[2], c[2], c[1], c[1]},
                                          vec<4, T, Q>{c[3], c[3], c[3], c[2]}};
    vec<4, T, Q> const dets = tmp[0] * (tmp[4] * tmp[8] - tmp[5] * tmp[7]) -
                              tmp[1] * (tmp[3] * tmp[8] - tmp[5] * tmp[6]) +
                              tmp[2] * (tmp[3] * tmp[7] - tmp[4] * tmp[6]);

    vec<4, T, Q> result{dets};
    result[1] *= -1;

    return result;
  }
};

template <typename T, qualifier Q, bool Aligned>
struct compute_orthogonal_basis<vec<4, T, Q>, T, Aligned> {

  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static mat<4, 4, T, Q>
  call(vec<4, T, Q> const &p, vec<4, T, Q> const &p1, vec<4, T, Q> const &p2) {

    vec<4, T, Q> const p3 = lorentz::cross(p, p1, p2);
    return mat<4, 4, T, Q>{p1, p2, p3, p};
  }

  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static mat<4, 4, T, Q>
  call(vec<4, T, Q> const &p, vec<4, T, Q> const &p1) {

    vec<4, T, Q> const v{static_cast<T>(0.1f), static_cast<T>(-0.15f),
                         static_cast<T>(0.2345f), static_cast<T>(0.95525376f)};

    T const vdotp = lorentz::dot(v, p);
    T const vdotp1 = lorentz::dot(v, p1);

    vec<4, T, Q> const dir = v + vdotp * p - vdotp1 * p1;
    T const L2 = lorentz::dot(dir, dir);
    vec<4, T, Q> const p2 = dir * inversesqrt(L2);

    return call(p, p1, p2);
  }
};

template <typename T, qualifier Q, bool Aligned>
struct compute_transpose<mat<4, 4, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static mat<4, 4, T, Q>
  call(mat<4, 4, T, Q> const &m) {
    mat<4, 4, T, Q> result = glm::transpose(m);
    result[3] *= static_cast<float>(-1);
    result[0][3] *= static_cast<float>(-1);
    result[1][3] *= static_cast<float>(-1);
    result[2][3] *= static_cast<float>(-1);
    result[3][3] *= static_cast<float>(-1);
    return result;
  }
};

template <typename T, qualifier Q, bool Aligned>
struct compute_transpose<mat<3, 3, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static mat<3, 3, T, Q>
  call(mat<3, 3, T, Q> const &m) {
    mat<3, 3, T, Q> result = glm::transpose(m);
    result[2] *= static_cast<float>(-1);
    result[0][2] *= static_cast<float>(-1);
    result[1][2] *= static_cast<float>(-1);
    result[2][2] *= static_cast<float>(-1);
    return result;
  }
};

template <typename T, qualifier Q, bool Aligned>
struct compute_transpose<mat<2, 2, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static mat<2, 2, T, Q>
  call(mat<2, 2, T, Q> const &m) {
    mat<2, 2, T, Q> result = glm::transpose(m);
    result[0][1] *= static_cast<float>(-1);
    result[1][0] *= static_cast<float>(-1);
    return result;
  }
};

template <length_t DA, length_t DB, typename T, qualifier Q, bool Aligned>
struct compute_outer_product {
  GLM_FUNC_QUALIFIER static
      typename glm::detail::outerProduct_trait<DA, DB, T, Q>::type
      call(vec<DA, T, Q> const &col, vec<DB, T, Q> const &row) {
    typename glm::detail::outerProduct_trait<DA, DB, T, Q>::type result =
        glm::outerProduct(col, row);
    result[DB - 1] *= static_cast<float>(-1);
    return result;
  }
};

template <length_t R, typename T, qualifier Q, bool Aligned>
struct perform_gram_schmidt<mat<4, R, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static mat<4, R, T, Q> &
  call(mat<4, R, T, Q> &m) {

    vec<3, T, Q> invnorms{};
    invnorms[2] = 1 / lorentz::dot(m[3], m[3]);

    m[2] -= invnorms[2] * lorentz::dot(m[3], m[2]) * m[3];

    invnorms[1] = 1 / lorentz::dot(m[2], m[2]);

    m[1] -= invnorms[2] * lorentz::dot(m[3], m[1]) * m[3] +
            invnorms[1] * lorentz::dot(m[2], m[1]) * m[2];

    invnorms[0] = 1 / lorentz::dot(m[1], m[1]);

    m[0] -= invnorms[2] * lorentz::dot(m[3], m[0]) * m[3] +
            invnorms[1] * lorentz::dot(m[2], m[0]) * m[2] +
            invnorms[0] * lorentz::dot(m[1], m[0]) * m[1];

    return m;
  }
};

template <length_t R, typename T, qualifier Q, bool Aligned>
struct perform_gram_schmidt_normalised<mat<4, R, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static mat<4, R, T, Q> &
  call(mat<4, R, T, Q> &m) {

    m[3] = lorentz::normalize(m[3], -1);

    m[2] += lorentz::dot(m[3], m[2]) * m[3];
    m[2] = lorentz::normalize(m[2], 1);

    m[1] += lorentz::dot(m[3], m[1]) * m[3] - lorentz::dot(m[2], m[1]) * m[2];
    m[1] = lorentz::normalize(m[1], 1);

    m[0] += lorentz::dot(m[3], m[0]) * m[3] - lorentz::dot(m[2], m[0]) * m[2] -
            lorentz::dot(m[1], m[0]) * m[1];
    m[0] = lorentz::normalize(m[0], 1);

    return m;
  }
};

template <typename T, qualifier Q, bool Aligned>
struct compute_boost<mat<4, 4, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER static mat<4, 4, T, Q> call(vec<3, T, Q> const &dir) {

    T d2 = glm::dot(dir, dir);
    T d = sqrt(d2);
    T sinchd = (d > 1e-6) ? sinh(d) / d : 1 + d2 / 3;
    T cosh_coeff = (d > 1e-5) ? (cosh(d) - 1) / d2 : (12 + d2) / 24;

    mat<4, 4, T, Q> S{0};
    S[3] = vec<4, T, Q>{dir, 0};
    S[0][3] = dir[0];
    S[1][3] = dir[1];
    S[2][3] = dir[2];

    mat<4, 4, T, Q> S2{glm::outerProduct(dir, dir)};
    S2[3][3] = d2;

    return mat<4, 4, T, Q>{1} + sinchd * S + cosh_coeff * S2;
  }

  GLM_FUNC_QUALIFIER static mat<4, 4, T, Q> call(mat<4, 4, T, Q> const &M,
                                                 vec<3, T, Q> const &dir) {
    return M * call(dir);
  }
};

template <typename T, qualifier Q, bool Aligned>
struct compute_boost<vec<4, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER static vec<4, T, Q> call(vec<4, T, Q> const &v, T dist,
                                              int axis) {

    T coshd = cosh(dist);
    T sinhd = sinh(dist);

    vec<4, T, Q> result{v};

    result[axis] = coshd * v[axis] + sinhd * v[3];
    result[3] = coshd * v[3] + sinhd * v[axis];

    return result;
  }
};

template <typename T, qualifier Q, bool Aligned>
struct compute_parabolic<mat<4, 4, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER static mat<4, 4, T, Q> call(T dist, int axis1, int axis2,
                                                 int sgn) {
    T hdd = dist * dist / 2;

    mat<4, 4, T, Q> result(1.0f);
    result[axis2][axis1] = -sgn * dist;
    result[axis1][axis2] = sgn * dist;
    result[axis1][3] = dist;
    result[3][axis1] = dist;
    result[axis2][axis2] = 1 - hdd;
    result[3][3] = 1 + hdd;
    result[axis2][3] = -sgn * hdd;
    result[3][axis2] = sgn * hdd;

    return result;
  }

  GLM_FUNC_QUALIFIER static mat<4, 4, T, Q>
  call(mat<4, 4, T, Q> const &M, T dist, int axis1, int axis2, int sgn) {
    return M * call(dist, axis1, axis2, sgn);
  }

  GLM_FUNC_QUALIFIER static vec<4, T, Q> call(vec<4, T, Q> const &v, T dist,
                                              int axis1, int axis2, int sgn) {
    return call(dist, axis1, axis2, sgn) * v;
  }
};

} // namespace detail

template <length_t L, typename T, qualifier Q>
GLM_FUNC_QUALIFIER GLM_CONSTEXPR T dot(vec<L, T, Q> const &x,
                                       vec<L, T, Q> const &y) {
  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'dot' accepts only floating-point inputs");
  return lorentz::detail::compute_dot<
      vec<L, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(x, y);
}

template <length_t L, typename T, qualifier Q>
GLM_FUNC_QUALIFIER vec<L, T, Q> normalize(vec<L, T, Q> const &x) {
  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'normalize' accepts only floating-point inputs");

  return detail::compute_normalize<L, T, Q,
                                   glm::detail::is_aligned<Q>::value>::call(x);
}

template <length_t L, typename T, qualifier Q>
GLM_FUNC_QUALIFIER vec<L, T, Q> normalize(vec<L, T, Q> const &x,
                                          int expected_sign) {
  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'normalize' accepts only floating-point inputs");

  return detail::compute_normalize<
      L, T, Q, glm::detail::is_aligned<Q>::value>::call(x, expected_sign);
}

template <length_t L, typename T, qualifier Q>
GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<L, T, Q>
cross(vec<L, T, Q> const &x, vec<L, T, Q> const &y, vec<L, T, Q> const &z) {
  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'cross' accepts only floating-point inputs");
  return lorentz::detail::compute_cross<
      vec<L, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(x, y, z);
}

template <length_t L, typename T, qualifier Q>
GLM_FUNC_QUALIFIER GLM_CONSTEXPR mat<L, L, T, Q>
complete_orthogonal_basis(vec<L, T, Q> const &p, vec<L, T, Q> const &d) {
  GLM_STATIC_ASSERT(
      std::numeric_limits<T>::is_iec559,
      "'complete_orthogonal_basis' accepts only floating-point inputs");

  return lorentz::detail::compute_orthogonal_basis<
      vec<L, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(p, d);
}

template <length_t L, typename T, qualifier Q>
GLM_FUNC_QUALIFIER mat<L, L, T, Q> transpose(mat<L, L, T, Q> const &m) {
  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'transpose' accepts only floating-point inputs");

  return lorentz::detail::compute_transpose<
      mat<L, L, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(m);
}

template <length_t L, typename T, qualifier Q>
GLM_FUNC_QUALIFIER vec<L, T, Q> transpose(vec<L, T, Q> const &v) {
  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'transpose' accepts only floating-point inputs");

  vec<L, T, Q> result{v};
  result[L - 1] = -result[L - 1];
  return result;
}

template <length_t DA, length_t DB, typename T, qualifier Q>
GLM_FUNC_QUALIFIER typename glm::detail::outerProduct_trait<DA, DB, T, Q>::type
outerProduct(vec<DA, T, Q> const &c, vec<DB, T, Q> const &r) {
  GLM_STATIC_ASSERT(
      std::numeric_limits<T>::is_iec559,
      "'lorentzian outer product' accepts only floating-point inputs");
  return detail::compute_outer_product<
      DA, DB, T, Q, glm::detail::is_aligned<Q>::value>::call(c, r);
}

template <length_t C, length_t R, typename T, qualifier Q>
GLM_FUNC_QUALIFIER GLM_CONSTEXPR mat<C, R, T, Q> &
gram_schmidt_in_place(mat<C, R, T, Q> &M) {

  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'boost' accepts only floating-point inputs");

  return lorentz::detail::perform_gram_schmidt<
      mat<C, R, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(M);
}

template <length_t C, length_t R, typename T, qualifier Q>
GLM_FUNC_QUALIFIER GLM_CONSTEXPR mat<C, R, T, Q>
gram_schmidt(mat<C, R, T, Q> const &M) {

  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'boost' accepts only floating-point inputs");
  auto result{M};
  return lorentz::detail::perform_gram_schmidt<
      mat<C, R, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(result);
}

template <length_t C, length_t R, typename T, qualifier Q>
GLM_FUNC_QUALIFIER GLM_CONSTEXPR mat<C, R, T, Q> &
gram_schmidt_normalised_in_place(mat<C, R, T, Q> &M) {

  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'boost' accepts only floating-point inputs");

  return lorentz::detail::perform_gram_schmidt_normalised<
      mat<C, R, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(M);
}

template <length_t C, length_t R, typename T, qualifier Q>
GLM_FUNC_QUALIFIER GLM_CONSTEXPR mat<C, R, T, Q>
gram_schmidt_normalised(mat<C, R, T, Q> const &M) {

  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'boost' accepts only floating-point inputs");
  auto result{M};
  return lorentz::detail::perform_gram_schmidt_normalised<
      mat<C, R, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(result);
}

template <length_t L, typename T, qualifier Q>
GLM_FUNC_QUALIFIER mat<L + 1, L + 1, T, Q> boost(vec<L, T, Q> const &dir) {
  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'boost' accepts only floating-point inputs");

  return lorentz::detail::compute_boost<
      mat<L + 1, L + 1, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(dir);
}

template <length_t L, length_t K, typename T, qualifier Q>
GLM_FUNC_QUALIFIER mat<L, L, T, Q> boost(mat<L, L, T, Q> const &M,
                                         vec<K, T, Q> const &dir) {
  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'boost' accepts only floating-point inputs");

  GLM_STATIC_ASSERT(
      L == K + 1,
      "'boost' vector should have dimension one less than the matrix");

  return lorentz::detail::compute_boost<
      mat<L, L, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(M, dir);
}

template <length_t L, typename T, qualifier Q>
GLM_FUNC_QUALIFIER vec<L, T, Q> boost(vec<L, T, Q> const &v, T dist, int axis) {
  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'boost' accepts only floating-point inputs");

  return lorentz::detail::compute_boost<
      vec<L, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(v, dist, axis);
}

template <length_t L, typename T, qualifier Q>
GLM_FUNC_QUALIFIER mat<L, L, T, Q> parabolic(mat<L, L, T, Q> const &M, T dist,
                                             int axis1, int axis2, int sgn) {
  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'parabolic' accepts only floating-point inputs");

  return lorentz::detail::compute_parabolic<
      mat<L, L, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(M, dist,
                                                                   axis1, axis2,
                                                                   sgn);
}

template <length_t L, typename T, qualifier Q>
GLM_FUNC_QUALIFIER vec<L, T, Q> parabolic(vec<L, T, Q> const &v, T dist,
                                          int axis1, int axis2, int sgn) {
  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'parabolic' accepts only floating-point inputs");

  return lorentz::detail::compute_parabolic<
      mat<L, L, T, Q>, T, glm::detail::is_aligned<Q>::value>::call(v, dist,
                                                                   axis1, axis2,
                                                                   sgn);
}

} // namespace lorentz

namespace detail {

template <typename M, typename T, bool Aligned> struct compute_trace {};

template <typename T, qualifier Q, bool Aligned>
struct compute_trace<mat<4, 4, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static T call(mat<4, 4, T, Q> const &m) {
    return (m[0][0] + m[1][1]) + (m[2][2] + m[3][3]);
  }
};

template <typename T, qualifier Q, bool Aligned>
struct compute_trace<mat<3, 3, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static T call(mat<3, 3, T, Q> const &m) {
    return m[0][0] + m[1][1] + m[2][2];
  }
};

template <typename T, qualifier Q, bool Aligned>
struct compute_trace<mat<2, 2, T, Q>, T, Aligned> {
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR static T call(mat<2, 2, T, Q> const &m) {
    return (m[0][0] + m[1][1]);
  }
};

} // namespace detail

template <length_t R, length_t C, typename T, qualifier Q>
GLM_FUNC_QUALIFIER GLM_CONSTEXPR T trace(mat<R, C, T, Q> const &m) {
  GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559,
                    "'trace' accepts only floating-point inputs");
  return detail::compute_trace<mat<R, C, T, Q>, T,
                               glm::detail::is_aligned<Q>::value>::call(m);
}

} // namespace glm