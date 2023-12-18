// Copyright (C) 2023 Alessandro Fornasier, Pieter van Goor.
// Control of Networked Systems, University of Klagenfurt, Austria.
// System Theory and Robotics Lab, Australian Centre for Robotic
// Vision, Australian national University, Australia.
//
// All rights reserved.
//
// This software is licensed under the terms of the BSD-2-Clause-License with
// no commercial use allowed, the full terms of which are made available
// in the LICENSE file. No license in patents is granted.
//
// You can contact the authors at <alessandro.fornasier@ieee.org>,
// <pieter.vangoor@anu.edu.au>.

#ifndef SEn3_HPP
#define SEn3_HPP

#include <array>

#include "SO3.hpp"

namespace group
{
/**
 * @brief the Special Eucldean group with n multiple isometries (SEn3)
 *
 * @tparam FPType. Floating point type (float, double, long double)
 *
 * @note Group Formulation for Consistent Non-Linear Estiamtion
 * @note State Estimation for Robotics [http://asrl.utias.utoronto.ca/~tdb/bib/barfoot_ser17.pdf]
 * @note The Geometry and Kinematics of the MatrixType Lie Group SEK(3) [https://arxiv.org/abs/2012.00950]
 */
template <typename FPType, int n>
class SEn3
{
  static_assert(n > 0, "SEn3: n needs to be grater than 0");

 public:
  using Scalar = FPType;        //!< The underlying scalar type
  using SO3Type = SO3<FPType>;  //!< The underlying SO3 type
  using VectorType =
      Eigen::Matrix<FPType, 3 + 3 * n, 1>;  //!< R6 Vectorspace element type (isomorphic to Lie Algebra se3)
  using MatrixType = Eigen::Matrix<FPType, 3 + n, 3 + n>;              //!< Lie Algebra / Lie Group matrix type
  using TMatrixType = Eigen::Matrix<FPType, 3 + 3 * n, 3 + 3 * n>;     //!< Transformation matrix type
  using IsometriesType = std::array<typename SO3Type::VectorType, n>;  //!< Vector of translations (Isometries)

  /**
   * @brief Construct an identity SEn3 object
   */
  SEn3() : C_(SO3Type()), t_() { t_.fill(SO3Type::VectorType::Zero()); }

  /**
   * @brief Construct a SEn3 object from a given normalized quaternion, and an array of vectors.
   *
   * @param q Quaternion
   * @param t Array of R3 vectors
   */
  SEn3(const typename SO3Type::QuaternionType& q, const IsometriesType& t) : C_(q), t_(t) {}

  /**
   * @brief Construct a SEn3 object from a given rotation matrix, and an array of vectors.
   *
   * @param R Rotation matrix
   * @param t Array of R3 vectors
   */
  SEn3(const typename SO3Type::MatrixType& R, const IsometriesType& t) : C_(R), t_(t) {}

  /**
   * @brief Construct a SEn3 object from a given SO3 object, and an array of vectors.
   *
   * @param C SO3 group element
   * @param t Array of R3 vectors
   */
  SEn3(const SO3Type& C, const IsometriesType& t) : C_(C), t_(t) {}

  /**
   * @brief Construct a SEn3 object from a given matrix
   *
   * @param T SEn3 group element in matrix form
   */
  SEn3(const MatrixType& T) : C_(T.template block<3, 3>(0, 0)), t_()
  {
    for (int i = 0; i < n; ++i)
    {
      t_[i] = T.template block<3, 1>(0, 3 + i);
    }
  }

  /**
   * @brief wedge operator, transform a vector in R(3+3n) to a matrix in sen3
   *
   * @param u R(3+3n) vector
   *
   * @return SEn3 Lie algebra element in matrix form
   */
  [[nodiscard]] static const MatrixType wedge(const VectorType& u)
  {
    MatrixType U = MatrixType::Zero();
    U.template block<3, 3>(0, 0) = SO3Type::wedge(u.template block<3, 1>(0, 0));
    for (int i = 0; i < n; ++i)
    {
      U.template block<3, 1>(0, 3 + i) = u.template block<3, 1>(3 + 3 * i, 0);
    }
    return U;
  }

  /**
   * @brief Transform a matrix in sen3 to a vector in R(3+3n)
   *
   * @param U SEn3 Lie algebra element in matrix form
   *
   * @return R(3+3n) vector
   */
  [[nodiscard]] static const VectorType vee(const MatrixType& U)
  {
    VectorType u = VectorType::Zero();
    u.template block<3, 1>(0, 0) = SO3Type::vee(U.template block<3, 3>(0, 0));
    for (int i = 0; i < n; ++i)
    {
      u.template block<3, 1>(3 + 3 * i, 0) = U.template block<3, 1>(0, 3 + i);
    }
    return u;
  }

  /**
   * @brief sen3 adjoint matrix
   *
   * @param u R(3+3n) vector
   *
   * @return SEn3 Lie algebra adjoint matrix
   */
  [[nodiscard]] static const TMatrixType adjoint(const VectorType& u)
  {
    TMatrixType ad = TMatrixType::Zero();
    typename SO3Type::MatrixType W = SO3Type::wedge(u.template block<3, 1>(0, 0));
    ad.template block<3, 3>(0, 0) = W;
    for (int i = 0; i < n; ++i)
    {
      ad.template block<3, 3>(3 + 3 * i, 3 + 3 * i) = W;
      ad.template block<3, 3>(3 + 3 * i, 0) = SO3Type::wedge(u.template block<3, 1>(3 + 3 * i, 0));
    }
    return ad;
  }

  /**
   * @brief SEn3 left Jacobian matrix
   *
   * @param u R(3+3n) vector
   *
   * @return SEn3 left Jacobian matrix
   */
  [[nodiscard]] static const TMatrixType leftJacobian(const VectorType& u)
  {
    typename SO3Type::VectorType w = u.template block<3, 1>(0, 0);
    FPType ang = w.norm();
    if (ang < eps_)
    {
      return TMatrixType::Identity() + 0.5 * adjoint(u);
    }
    typename SO3Type::MatrixType SO3JL = SO3Type::leftJacobian(w);
    TMatrixType J = TMatrixType::Identity();
    J.template block<3, 3>(0, 0) = SO3JL;
    for (int i = 0; i < n; ++i)
    {
      typename SO3Type::VectorType x = u.template block<3, 1>(3 + 3 * i, 0);
      J.template block<3, 3>(3 + 3 * i, 0) = SE3leftJacobianQ(w, x);
      J.template block<3, 3>(3 + 3 * i, 3 + 3 * i) = SO3JL;
    }
    return J;
  }

  /**
   * @brief SEn3 right Jacobian matrix
   *
   * @param u R(3+3n) vector
   *
   * @return SEn3 right Jacobian matrix
   */
  [[nodiscard]] static const TMatrixType rightJacobian(const VectorType& u) { return leftJacobian(-u); }

  /**
   * @brief The exponential map for SEn3.
   * Returns a SEn3 object given a vector u in R(3+3n) (equivalent to exp(wedge(u)))
   *
   * @param u R(3+3n) vector
   *
   * @return SEn3 group element
   */
  [[nodiscard]] static const SEn3 exp(const VectorType& u)
  {
    typename SO3Type::VectorType w = u.template block<3, 1>(0, 0);
    typename SO3Type::MatrixType SO3JL = SO3Type::leftJacobian(w);
    IsometriesType t;
    for (int i = 0; i < n; ++i)
    {
      t[i] = SO3JL * u.template block<3, 1>(3 + 3 * i, 0);
    }
    return SEn3(SO3Type::exp(w), t);
  }

  // /**
  //  * @brief The exponential map for SEn3.
  //  * Returns a SEn3 object given a matrix U in sen3
  //  *
  //  * @param U SEn3 Lie algebra element in matrix form
  //  *
  //  * @return SEn3 group element
  //  */
  // [[nodiscard]] static SEn3 exp(const MatrixType& U) { return exp(vee(U)); }

  /**
   * @brief The logarithmic map for SEn3.
   * Return a vector given a SEn3 object (equivalent to vee(log(X)))
   *
   * @param X SEn3 group element
   *
   * @return R(3+3n) vector
   */
  [[nodiscard]] static const VectorType log(const SEn3& X)
  {
    VectorType u = VectorType::Zero();
    u.template block<3, 1>(0, 0) = SO3Type::log(X.C_);
    typename SO3Type::MatrixType invSO3JL = SO3Type::leftJacobian(u.template block<3, 1>(0, 0)).inverse();
    for (int i = 0; i < n; ++i)
    {
      u.template block<3, 1>(3 + 3 * i, 0) = invSO3JL * X.t_[i];
    }
    return u;
  }

  // /**
  //  * @brief The logarithmic map for SEn3.
  //  * Return a se3 matrix given a SEn3 object
  //  *
  //  * @param X SEn3 group element
  //  *
  //  * @return SEn3 Lie algebra element in matrix form
  //  */
  // [[nodiscard]] static MatrixType log(const SEn3& X) { return wedge(log(X)); }

  /**
   * @brief Operator * overloading.
   * Implements the SEn3 composition this * other
   *
   * @param other SEn3 group element
   *
   * @return SEn3 group element
   *
   * @note usage: z = x * y
   */
  [[nodiscard]] const SEn3 operator*(const SEn3& other) const
  {
    IsometriesType t;
    for (int i = 0; i < n; ++i)
    {
      t[i] = C_.R() * other.t_[i] + t_[i];
    }
    return SEn3(C_ * other.C_, t);
  }

  /**
   * @brief Operator * overloading.
   * Implements the SEn3 composition this * other with a SEn3 group or alegra element in matrix form
   *
   * @param other SE3 group or algebra element in matrix form
   *
   * @return SEn3 group or algebra element in matrix form
   *
   * @note usage: z = x * y
   */
  [[nodiscard]] const MatrixType operator*(const MatrixType& other) const
  {
    bool is_algebra = true;
    bool is_group = true;

    for (int i = 0; i < n; ++i)
    {
      is_algebra &= (other(3 + i, 3 + i) == 0);
      is_group &= (other(3 + i, 3 + i) == 1);
    }

    if (!(is_algebra || is_group))
    {
      throw std::runtime_error(
          "SEn3: operator* is defined only for composition with matrix form of SEn3 group or algebra elements");
    }

    MatrixType res = other;
    res.template block<3, 3>(0, 0) = C_.R() * other.template block<3, 3>(0, 0);
    for (int i = 0; i < n; ++i)
    {
      if (is_algebra)
      {
        res.template block<3, 1>(0, 3 + i) = C_.R() * other.template block<3, 1>(0, 3 + i);
      }
      else
      {
        res.template block<3, 1>(0, 3 + i) = C_.R() * other.template block<3, 1>(0, 3 + i) + t_[i];
      }
    }
    return res;
  }

  /**
   * @brief Operator * overloading.
   * Implements the SE3 action on a R3 vector
   *
   * @param other R3 vector
   *
   * @return R3 vector
   */
  [[nodiscard]] const typename SO3Type::VectorType operator*(const typename SO3Type::VectorType& other) const
  {
    static_assert(n == 1, "SEn3: SE3 action on R3 (* operator overloading) available only for n = 1");
    return C_.R() * other + t_[0];
  }

  /**
   * @brief Implements the SEn3 composition this = this * other
   *
   * @param other SEn3 group element
   *
   * @return SEn3 group element
   */
  const SEn3& multiplyRight(const SEn3& other)
  {
    for (int i = 0; i < n; ++i)
    {
      t_[i] = (C_.R() * other.t_[i] + t_[i]).eval();
    }
    C_.multiplyRight(other.C_);  // C_ * other.C_
    return *this;
  }

  /**
   * @brief Implements the SEn3 composition this = other * this
   *
   * @param other SEn3 group element
   *
   * @return SEn3 group element
   */
  const SEn3& multiplyLeft(const SEn3& other)
  {
    for (int i = 0; i < n; ++i)
    {
      t_[i] = (other.C_.R() * t_[i] + other.t_[i]).eval();
    }
    C_.multiplyLeft(other.C_);  // other.C_ * th.C_
    return *this;
  }

  /**
   * @brief Get a constant copy of the inverse of the SEn3 object
   *
   * @return SEn3 group element
   */
  [[nodiscard]] const SEn3 inv() const
  {
    IsometriesType t;
    for (int i = 0; i < n; ++i)
    {
      t[i] = -C_.R().transpose() * t_[i];
    }
    return SEn3(C_.R().transpose(), t);
  }

  /**
   * @brief Get a constant copy of the SEn3 object as a matrix
   *
   * @return SEn3 group element in matrix form
   */
  [[nodiscard]] const MatrixType T() const
  {
    MatrixType T = MatrixType::Identity();
    T.template block<3, 3>(0, 0) = C_.R();
    for (int i = 0; i < n; ++i)
    {
      T.template block<3, 1>(0, 3 + i) = t_[i];
    }
    return T;
  }

  /**
   * @brief Get a constant reference to the SEn3 rotation matrix
   *
   * @return Rotation matrix
   */
  [[nodiscard]] const typename SO3Type::MatrixType& R() const { return C_.R(); }

  /**
   * @brief Get a constant reference to the SEn3 normalized quaternion
   *
   * @return Quaternion
   */
  [[nodiscard]] const typename SO3Type::QuaternionType& q() const { return C_.q(); }

  /**
   * @brief Get a constant reference to the SEn3 translation vectors (isometries)
   *
   * @return Array of R3 vectors
   */
  [[nodiscard]] const IsometriesType& t() const { return t_; }

  /**
   * @brief Get a constant referece to the first isometry of SE3
   *
   * @note Method available only for n = 1, thus SE3
   *
   * @return R3 vector
   */
  [[nodiscard]] const typename SO3Type::VectorType& x() const
  {
    static_assert(n == 1, "SEn3: x() method available only for n = 1");
    return t_[0];
  }

  /**
   * @brief Get a constant referece to the first isometry (velocity) of SEn3 with n > 1
   *
   * @note Method available only for n > 1
   *
   * @return R3 vector
   */
  [[nodiscard]] const typename SO3Type::VectorType& v() const
  {
    static_assert(n > 1, "SEn3: v() method available only for n > 1");
    return t_[0];
  }

  /**
   * @brief Get a constant referece to the second isometry (position) of SEn3 with n > 1
   *
   * @note Method available only for n > 1
   *
   * @return R3 vector
   */
  [[nodiscard]] const typename SO3Type::VectorType& p() const
  {
    static_assert(n > 1, "SEn3: v() method available only for n > 1");
    return t_[1];
  }

  /**
   * @brief Get a constant copy of the SEn3 object as a matrix
   *
   * @return SEn3 group element in matrix form
   */
  [[nodiscard]] const MatrixType asMatrix() const { return T(); }

  /**
   * @brief Set SEn3 object value from given matrix
   *
   * @param T SEn3 group element in matrix form
   */
  void fromT(const MatrixType& T)
  {
    C_.fromR(T.template block<3, 3>(0, 0));
    for (int i = 0; i < n; ++i)
    {
      t_[i] = T.template block<3, 1>(0, 3 + i);
    }
  }

  /**
   * @brief SEn3 Adjoint matrix
   *
   * @return SEn3 group Adjoint matrix
   */
  [[nodiscard]] const TMatrixType Adjoint() const
  {
    TMatrixType Ad = TMatrixType::Identity();
    Ad.template block<3, 3>(0, 0) = C_.R();
    for (int i = 0; i < n; ++i)
    {
      Ad.template block<3, 3>(3 + 3 * i, 3 + 3 * i) = C_.R();
      Ad.template block<3, 3>(3 + 3 * i, 0) = SO3Type::wedge(t_[i]) * C_.R();
    }
    return Ad;
  }

  /**
   * @brief SEn3 Inverse Adjoint matrix
   *
   * @return SEn3 group inverse Adjoint matrix
   */
  [[nodiscard]] const TMatrixType invAdjoint() const
  {
    TMatrixType Ad = TMatrixType::Identity();
    Ad.template block<3, 3>(0, 0) = C_.R().transpose();
    for (int i = 0; i < n; ++i)
    {
      Ad.template block<3, 3>(3 + 3 * i, 3 + 3 * i) = C_.R().transpose();
      Ad.template block<3, 3>(3 + 3 * i, 0) = -C_.R().transpose() * SO3Type::wedge(t_[i]);
    }
    return Ad;
  }

 protected:
  /**
   * @brief SE3 left Jacobian Q matrix
   *
   * @param w R3 vector
   * @param v R3 vector
   *
   * @return SE3 left Jacobian Q matrix
   */
  [[nodiscard]] static const typename SO3Type::MatrixType SE3leftJacobianQ(const typename SO3Type::VectorType& w,
                                                                           const typename SO3Type::VectorType& v)
  {
    typename SO3Type::MatrixType p = SO3Type::wedge(w);
    typename SO3Type::MatrixType r = SO3Type::wedge(v);

    FPType ang = w.norm();
    FPType s = sin(ang);
    FPType c = cos(ang);

    FPType ang_p2 = pow(ang, 2);
    FPType ang_p3 = ang_p2 * ang;
    FPType ang_p4 = ang_p3 * ang;
    FPType ang_p5 = ang_p4 * ang;

    FPType c1 = (ang - s) / ang_p3;
    FPType c2 = (0.5 * ang_p2 + c - 1.0) / ang_p4;
    FPType c3 = (ang * (1.0 + 0.5 * c) - 1.5 * s) / ang_p5;

    typename SO3Type::MatrixType m1 = p * r + r * p + p * r * p;
    typename SO3Type::MatrixType m2 = p * p * r + r * p * p - 3.0 * p * r * p;
    typename SO3Type::MatrixType m3 = p * r * p * p + p * p * r * p;

    return 0.5 * r + c1 * m1 + c2 * m2 + c3 * m3;
  }

  SO3Type C_;         //!< SO3 object the rotational component of SEn3
  IsometriesType t_;  //!< R6 vector representing the n translational components of SEn3

  static constexpr FPType eps_ = std::is_same_v<FPType, float> ? 1.0e-6 : 1.0e-9;  //!< Epsilon
};

using SE3d = SEn3<double, 1>;   //!< The SE3 group with double precision floating point
using SE3f = SEn3<float, 1>;    //!< The SE3 group with single precision floating point
using SE23d = SEn3<double, 2>;  //!< The SE23 group with double precision floating point
using SE23f = SEn3<float, 2>;   //!< The SE23 group with single precision floating point

}  // namespace group

#endif  // SEn3_HPP