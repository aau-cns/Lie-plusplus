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

#ifndef SO3_HPP
#define SO3_HPP

#include <Eigen/Dense>

namespace group
{
/**
 * @brief The Special orthogonal group of dimension 3 (SO3) representing 3D rotations
 *
 * @tparam FPType. Floating point type (float, FPType, long FPType)
 *
 * @note Group Formulation for Consistent Non-Linear Estiamtion
 * @note State Estimation for Robotics [http://asrl.utias.utoronto.ca/~tdb/bib/barfoot_ser17.pdf]
 * @note The Geometry and Kinematics of the Matrix Lie Group SEK(3) [https://arxiv.org/abs/2012.00950]
 */
template <typename FPType>
class SO3
{
 public:
  using Scalar = FPType;                             //!< The underlying scalar type
  using MatrixType = Eigen::Matrix<FPType, 3, 3>;    //!< Lie Algebra / Lie Group matrix type
  using TMatrixType = Eigen::Matrix<FPType, 3, 3>;   //!< Transformation matrix type (Linear operator on R3)
  using VectorType = Eigen::Matrix<FPType, 3, 1>;    //!< R3 Vectorspace element type (isomorphic to Lie Algebra so3)
  using QuaternionType = Eigen::Quaternion<FPType>;  //!< Quaternion type
  using AngleAxisType = Eigen::AngleAxis<FPType>;    //!< Angle-axis type

  /**
   * @brief Construct an identity SO3 object
   */
  SO3() : R_(MatrixType::Identity()), q_(QuaternionType::Identity()) { checkq(); }

  /**
   * @brief Construct a SO3 object from a given normalized quaternion.
   *
   * @param q Quaternion
   */
  SO3(const QuaternionType& q) : R_(q.toRotationMatrix()), q_(q) { checkq(); }

  /**
   * @brief Construct a SO3 object from a given rotation matrix
   *
   * @param R Rotation matrix
   */
  SO3(const MatrixType& R) : R_(R), q_(R) { checkq(); }

  /**
   * @brief Construct a SO3 object from two vector such that Ru = v
   *
   * @param u Vector in R3
   * @param v vector in R3
   */
  SO3(const VectorType& u, const VectorType& v) : R_(), q_()
  {
    VectorType un = u.normalized();
    VectorType vn = v.normalized();

    VectorType ax = un.cross(vn);
    FPType s = ax.norm();
    FPType c = un.dot(vn);

    if (s == 0)
    {
      q_ = QuaternionType::Identity();
      R_ = MatrixType::Identity();
    }
    else if (std::abs(1 + c) < eps_)
    {
      ax = un.cross(VectorType::Random());
      ax.normalize();
      q_ = QuaternionType(0.0, ax(0), ax(1), ax(2));
      R_ = q_.toRotationMatrix();
    }
    else
    {
      R_ = MatrixType::Identity() + wedge(ax) + ((1.0 - c) / (s * s)) * (wedge(ax) * wedge(ax));
      q_ = QuaternionType(R_);
    }

    checkq();
  }

  /**
   * @brief wedge operator, transform a vector in R3 to a matrix in so3
   *
   * @param u R3 vector
   * @return SO3 Lie algebra element in matrix form
   */
  [[nodiscard]] static const MatrixType wedge(const VectorType& u)
  {
    MatrixType U;
    U << 0.0, -u(2), u(1), u(2), 0.0, -u(0), -u(1), u(0), 0.0;
    return (U);
  }

  /**
   * @brief Transform a matrix in so3 to a vector in R3
   *
   * @param U SO3 Lie algebra element in matrix form
   *
   * @return R3 vector
   */
  [[nodiscard]] static const VectorType vee(const MatrixType& U)
  {
    VectorType u;
    u << U(2, 1), U(0, 2), U(1, 0);
    return (u);
  }

  /**
   * @brief so3 adjoint matrix
   *
   * @param u R3 vector
   *
   * @return SO3 Lie algebra adjoint matrix
   */
  [[nodiscard]] static const TMatrixType adjoint(const VectorType& u) { return wedge(u); }

  /**
   * @brief SO3 Gamma2 matrix
   *
   * @param u R3 vector
   *
   * @return SO3 Gamma2 matrix
   */
  [[nodiscard]] static const TMatrixType Gamma2(const VectorType& u)
  {
    FPType ang = u.norm();
    if (ang < eps_)
    {
      return 0.5 * TMatrixType::Identity() - 1 / 6 * wedge(u);
    }
    VectorType ax = u / ang;
    FPType ang_p2 = pow(ang, 2);
    FPType s = (ang - sin(ang)) / ang_p2;
    FPType c = (1 - cos(ang)) / ang_p2;
    return c * TMatrixType::Identity() + s * wedge(ax) + (0.5 - c) * ax * ax.transpose();
  }

  /**
   * @brief SO3 left Jacobian matrix
   *
   * @param u R3 vector
   *
   * @return SO3 left Jacobian matrix
   */
  [[nodiscard]] static const TMatrixType leftJacobian(const VectorType& u)
  {
    FPType ang = u.norm();
    if (ang < eps_)
    {
      return TMatrixType::Identity() + 0.5 * wedge(u);
    }
    VectorType ax = u / ang;
    FPType s = sin(ang) / ang;
    FPType c = cos(ang);
    return s * TMatrixType::Identity() + ((1.0 - c) / ang) * wedge(ax) + (1.0 - s) * ax * ax.transpose();
  }

  /**
   * @brief SO3 right Jacobian matrix
   *
   * @param u R3 vector
   *
   * @return SO3 right Jacobian matrix
   */
  [[nodiscard]] static const TMatrixType rightJacobian(const VectorType& u) { return leftJacobian(-u); }

  /**
   * @brief The exponential map for SO3.
   * Returns a SO3 object given a vector u in R3 (equivalent to exp(wedge(u)))
   *
   * @param u R3 vector
   *
   * @return SO3 group element
   */
  [[nodiscard]] static const SO3 exp(const VectorType& u)
  {
    FPType ang = 0.5 * u.norm();
    QuaternionType q;
    if (ang < eps_)
    {
      q.w() = 1.0;
      q.vec() = 0.5 * u;
    }
    else
    {
      q.w() = cos(ang);
      q.vec() = sin(ang) * u.normalized();
    }
    q.normalize();
    return SO3(q);
  }

  // /**
  //  * @brief The exponential map for SO3.
  //  * Returns a SO3 object given a matrix U in so3
  //  *
  //  * @param U SO3 Lie algebra element in matrix form
  //  *
  //  * @return SO3 group element
  //  */
  // [[nodiscard]] static const SO3 exp(const MatrixType& U) { return exp(vee(U)); }

  /**
   * @brief The logarithmic map for SO3.
   * Return a vector given a SO3 object (equivalent to vee(log(X)))
   *
   * @param X SO3 group element
   *
   * @return R3 vector
   */
  [[nodiscard]] static const VectorType log(const SO3& X)
  {
    FPType qw = X.q_.w();
    VectorType qv = X.q_.vec();
    FPType ang = qv.norm();
    VectorType u;
    if (qw < 0)
    {
      qw = -qw;
      qv = -qv;
    }
    if (ang < eps_)
    {
      u = (2.0 / qw) * qv * (1.0 - (pow((ang / qw), 2) / 3));
    }
    else
    {
      u = 2 * atan2(ang, qw) * (qv / ang);
    }
    return u;
  }

  // /**
  //  * @brief The logarithmic map for SO3.
  //  * Return a so3 matrix given a SO3 object
  //  *
  //  * @param X SO3 group element
  //  *
  //  * @return SO3 Lie algebra element in matrix form
  //  */
  // [[nodiscard]] static const MatrixType log(const SO3& X) { return wedge(log(X)); }

  /**
   * @brief Operator * overloading.
   * Implements the SO3 composition this * other
   *
   * @param other SO3 group element
   *
   * @return SO3 group element
   *
   * @note usage: z = x * y
   */
  [[nodiscard]] const SO3 operator*(const SO3& other) const { return SO3(q_ * other.q_); }

  /**
   * @brief Operator * overloading.
   * Implements the SO3 composition this * other with a SO3 group or alegra element in matrix form
   *
   * @param other SO3 group or algebra element in matrix form
   *
   * @return SO3 group or algebra element in matrix form
   *
   * @note usage: z = x * y
   */
  [[nodiscard]] const MatrixType operator*(const MatrixType& other) const { return R_ * other; }

  /**
   * @brief Operator * overloading.
   * Implements the SO3 action on a R3 vector
   *
   * @param other R3 vector
   *
   * @return R3 vector
   */
  [[nodiscard]] const VectorType operator*(const VectorType& other) const
  {
    return q_ * other;
    // return R_ * other;
  }

  /**
   * @brief Implements the SO3 composition this = this * other
   *
   * @param other SO3 group element
   *
   * @return SO3 group element
   */
  const SO3& multiplyRight(const SO3& other)
  {
    q_ = q_ * other.q_;
    R_ = R_ * other.R_;
    return *this;
  }

  /**
   * @brief Implements the SO3 composition this = other * this
   *
   * @param other SO3 group element
   *
   * @return SO3 group element
   */
  const SO3& multiplyLeft(const SO3& other)
  {
    q_ = other.q_ * q_;
    R_ = other.R_ * R_;
    return *this;
  }

  /**
   * @brief Get a constant copy of the inverse of the SO3 object
   *
   * @return SO3 group element
   */
  [[nodiscard]] const SO3 inv() const { return SO3(q_.inverse()); }

  /**
   * @brief Get a constant reference to the SO3 rotation matrix
   *
   * @return Rotation matrix
   */
  [[nodiscard]] const MatrixType& R() const { return R_; }

  /**
   * @brief Get a constant reference to the SO3 normalized quaternion
   *
   * @return Quaternion
   */
  [[nodiscard]] const QuaternionType& q() const { return q_; }

  /**
   * @brief Get a constant reference to the SO3 object as a matrix
   *
   * @return SO3 group element in matrix form
   */
  [[nodiscard]] const MatrixType& asMatrix() const { return R_; }

  /**
   * @brief Set SO3 object value from a given normalized quaternion
   *
   * @param q quaternion
   */
  void fromq(const QuaternionType& q)
  {
    q_ = q;
    R_ = q.toRotationMatrix();

    checkq();
  }

  /**
   * @brief Set SO3 object value from given rotation matrix
   *
   * @param R rotation matrix
   */
  void fromR(const MatrixType& R)
  {
    q_ = R;
    R_ = R;

    checkq();
  }

  /**
   * @brief SO3 Adjoint matrix
   *
   * @return SO3 group Adjoint matrix
   */
  [[nodiscard]] const TMatrixType Adjoint() const { return R_; }

  /**
   * @brief SO3 Inverse Adjoint matrix
   *
   * @return SO3 group inverse Adjoint matrix
   */
  [[nodiscard]] const TMatrixType invAdjoint() const { return R_.transpose(); }

 protected:
  /**
   * @brief Check that q_ is a normalized quaternion
   */
  void checkq()
  {
    if (std::abs(q_.norm() - 1.0) > eps_)
    {
      throw std::invalid_argument("SO3: QuaternionType has to be normalized!");
    }
  }

  MatrixType R_;      //!< Rotation matrix
  QuaternionType q_;  //!< Normalized quaternion

  static constexpr FPType eps_ = std::is_same_v<FPType, float> ? 1.0e-6 : 1.0e-9;  //!< Epsilon
};

using SO3d = SO3<double>;  //!< The SO3 group with double precision floating point
using SO3f = SO3<float>;   //!< The SO3 group with single precision floating point

}  // namespace group

#endif  // SO3_HPP