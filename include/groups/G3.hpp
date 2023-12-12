// Copyright (C) 2023 Alessandro Fornasier, Giulio Delama.
// Control of Networked Systems, University of Klagenfurt, Austria.
//
// All rights reserved.
//
// This software is licensed under the terms of the BSD-2-Clause-License with
// no commercial use allowed, the full terms of which are made available
// in the LICENSE file. No license in patents is granted.
//
// You can contact the authors at <alessandro.fornasier@ieee.org>,
// <giulio.delama@ieee.org>.

#ifndef G3_HPP
#define G3_HPP

#include <array>

#include "SO3.hpp"
#include "SEn3.hpp"

namespace group
{
/**
 * @brief the Inhomogeneous Galileian group (G3)
 *
 * @tparam FPType. Floating point type (float, double, long double)
 *
 * @note Group Formulation for Consistent Non-Linear Estiamtion
 * @note /todo: add reference
 * @note /todo: add reference
 */
template <typename FPType>
class G3
{
 public:
  using Scalar = FPType;        //!< The underlying scalar type
  using SO3Type = SO3<FPType>;  //!< The underlying SO3 type
  using SE23Type = SEn3<FPType, 2>;                   //!< The underlying SE23 type
  using VectorType = Eigen::Matrix<FPType, 10, 1>;    //!< R10 Vectorspace element type (isomorphic to Lie Algebra g3)
  using MatrixType = Eigen::Matrix<FPType, 5, 5>;     //!< Lie Algebra / Lie Group matrix type
  using TMatrixType = Eigen::Matrix<FPType, 10, 10>;  //!< Transformation matrix type

  /**
   * @brief Construct an identity G3 object
   */
  G3() : D_(SE23Type()), s_(0.0) {}

  /**
   * @brief Construct a G3 object from a given normalized quaternion, an array of vectors, and a scalar factor.
   *
   * @param q Quaternion
   * @param t Array of R3 vectors
   * @param s Scalar factor
   */
  G3(const typename SO3Type::QuaternionType& q, const SE23Type::IsometriesType& t, const Scalar& s) : D_(q, t), s_(s) {}

  /**
   * @brief Construct a G3 object from a given rotation matrix, an array of vectors, and a scalar factor.
   *
   * @param R Rotation matrix
   * @param t Array of R3 vectors
   * @param s Scalar factor
   */
  G3(const typename SO3Type::MatrixType& R, const SE23Type::IsometriesType& t, const Scalar& s) : D_(R, t), s_(s) {}

  /**
   * @brief Construct a G3 object from a given SO3 object, an array of vectors, and a scalar factor.
   *
   * @param C SO3 group element
   * @param t Array of R3 vectors
   * @param s Scalar factor
   */
  G3(const SO3Type& C, const SE23Type::IsometriesType& t, const Scalar& s) : D_(C, t), s_(s) {}

  /**
   * @brief Construct a G3 object from a given SE23 object, and a scalar factor.
   *
   * @param D SE23 group element
   * @param s Scalar factor
   */
  G3(const SE23Type& D,const Scalar& s) : D_(D), s_(s) {}

  /**
   * @brief Construct a G3 object from a given matrix
   *
   * @param T G3 group element in matrix form
   */
  G3(const MatrixType& T) : D_(T), s_(T(3, 4)){}

  /**
   * @brief wedge operator, transform a vector in R10 to a matrix in g3
   *
   * @param u R10 vector
   *
   * @return g3 Lie algebra element in matrix form
   */
  [[nodiscard]] static const MatrixType wedge(const VectorType& u)
  {
    MatrixType U = MatrixType::Zero();
    U = SE23Type::wedge(u.template block<9, 1>(0, 0));
    U(3, 4) = u.template block<1, 1>(9, 0);
    return U;
  }

  /**
   * @brief Transform a matrix in g3 to a vector in R10
   *
   * @param U G3 Lie algebra element in matrix form
   *
   * @return R10 vector
   */
  [[nodiscard]] static const VectorType vee(const MatrixType& U)
  {
    VectorType u = VectorType::Zero();
    u.template block<9, 1>(0, 0) = SE23Type::vee(U);
    u.template block<1, 1>(9, 0) = U(3, 4);
    return u;
  }

  /**
   * @brief g3 adjoint matrix
   *
   * @param u R10 vector
   *
   * @return G3 Lie algebra adjoint matrix
   */
  [[nodiscard]] static const TMatrixType adjoint(const VectorType& u)
  {
    TMatrixType ad = TMatrixType::Zero();
    ad.template block<9, 9>(0, 0) = SE23Type::adjoint(u.template block<9, 1>(0, 0));
    ad.template block<3, 3>(6, 3) = -u.template block<1, 1>(9, 0) * SO3Type::MatrixType::Identity();
    ad.template block<3, 1>(6, 9) = u.template block<3, 1>(3, 0);
    return ad;
  }

  /**
   * @brief G3 left Jacobian matrix
   *
   * @param u R10 vector
   *
   * @return G3 left Jacobian matrix
   */
  [[nodiscard]] static const TMatrixType leftJacobian(const VectorType& u)
  {
    TMatrixType J = TMatrixType::Identity();
    J.template block<9, 9>(0, 0) = SE23Type::leftJacobian(u.template block<9, 1>(0, 0));
    J.template block<3, 3>(6, 0) -= u.template block<1, 1>(9, 0) * G3leftJacobianQ2(u.template block<3, 1>(0, 0), u.template block<3, 1>(3, 0));
    J.template block<3, 3>(6, 3) = -u.template block<1, 1>(9, 0) * SO3Type::MatrixType::Identity();
    J.template block<3, 1>(6, 9) = SO3Type::J2(u.template block<3, 1>(0, 0)) * u.template block<3, 1>(3, 0);
    return J;
  }

  /**
   * @brief G3 right Jacobian matrix
   *
   * @param u R10 vector
   *
   * @return G3 right Jacobian matrix
   */
  [[nodiscard]] static const TMatrixType rightJacobian(const VectorType& u) { return leftJacobian(-u); }

  /**
   * @brief The exponential map for G3.
   * Returns a G3 object given a vector u in R10 (equivalent to exp(wedge(u)))
   *
   * @param u R10 vector
   *
   * @return G3 group element
   */
  [[nodiscard]] static const G3 exp(const VectorType& u)
  {
    MatrixType expU = SE23Type.exp(u.template block<9, 1>(0, 0)).asMatrix();
    expU.template block<3, 1>(0, 4) += u.template block<1, 1>(9, 0) * SO3Type::J2(u.template block<3, 1>(0, 0)) * u.template block<3, 1>(3, 0);
    expU(3, 4) = u.template block<1, 1>(9, 0);
    return G3(expU);
  }

  // /**
  //  * @brief The exponential map for G3.
  //  * Returns a G3 object given a matrix U in g3
  //  *
  //  * @param U G3 Lie algebra element in matrix form
  //  *
  //  * @return G3 group element
  //  */
  // [[nodiscard]] static G3 exp(const MatrixType& U) { return exp(vee(U)); }

  /**
   * @brief The logarithmic map for G3.
   * Return a vector given a G3 object (equivalent to vee(log(X)))
   *
   * @param X G3 group element
   *
   * @return R10 vector
   */
  [[nodiscard]] static const VectorType log(const G3& X)
  {
    VectorType u = VectorType::Zero();
    u.template block<3, 1>(0, 0) = SO3Type::log(X.D_.C_);
    typename SO3Type::MatrixType invSO3JL = SO3Type::leftJacobian(u.template block<3, 1>(0, 0)).inverse();
    u.template block<3, 1>(3, 0) = invSO3JL * X.D_.v();
    u.template block<3, 1>(6, 0) = invSO3JL * (X.D_.p() - X.s_ * SO3Type::J2(u.template block<3, 1>(0, 0)) * u.template block<3, 1>(3, 0));
    u.template block<1, 1>(9, 0) = X.s_;
    return u;
  }

  // /**
  //  * @brief The logarithmic map for G3.
  //  * Return a g3 matrix given a G3 object
  //  *
  //  * @param X G3 group element
  //  *
  //  * @return G3 Lie algebra element in matrix form
  //  */
  // [[nodiscard]] static MatrixType log(const G3& X) { return wedge(log(X)); }

  /**
   * @brief Operator * overloading.
   * Implements the G3 composition this * other
   *
   * @param other G3 group element
   *
   * @return G3 group element
   *
   * @note usage: z = x * y
   */
  [[nodiscard]] const G3 operator*(const G3& other) const
  {
    MatrixType D = (D_ * other.D_).asMatrix();
    D.template block<3, 1>(0, 4) += D_.v() * other.s_;
    D.template block<1, 1>(3, 4) = other.s_ + s_;
    return G3(D);
  }

  /**
   * @brief Operator * overloading.
   * Implements the G3 composition with a g3 element this * other
   *
   * @param other G3 Lie algebra element in matrix form
   *
   * @return G3 Lie algebra element in matrix form
   *
   * @note usage: z = x * y
   */
  [[nodiscard]] const MatrixType operator*(const MatrixType& other) const
  {
    MatrixType D = D_ * other;
    D.template block<3, 1>(0, 4) += D_.v() * other.template block<1, 1>(3, 4);
    D.template block<1, 1>(3, 4) = other.template block<1, 1>(3, 4) + s_;
    return D;
  }

  /**
   * @brief Implements the G3 composition this = this * other
   *
   * @param other G3 group element
   *
   * @return G3 group element
   */
  const G3& multiplyRight(const G3& other)
  {
    D_.multiplyRight(other.D_);
    D_.v() += D_.v() * other.s_;
    s_ += other.s_;
    return *this;
  }

  /**
   * @brief Implements the G3 composition this = other * this
   *
   * @param other G3 group element
   *
   * @return G3 group element
   */
  const G3& multiplyLeft(const G3& other)
  {
    D_.multiplyLeft(other.D_);
    D_.v() += other.D_.v() * s_;
    s_ += other.s_;
    return *this;
  }

  /**
   * @brief Get a constant copy of the inverse of the G3 object
   *
   * @return G3 group element
   */
  [[nodiscard]] const G3 inv() const
  {
    MatrixType invG = D_.inv().asMatrix();
    invG.template block<3, 1>(0, 4) += D_.R().transpose() * D_.v() * s_;
    invG.template block<1, 1>(3, 4) = -s_;
    return G3(invG);
  }

  /**
   * @brief Get a constant copy of the SE23 object as a matrix
   *
   * @return G3 group element in matrix form
   */
  [[nodiscard]] const MatrixType T() const { return D_.T(); }

  /**
   * @brief Get a constant reference to the SE23 rotation matrix
   *
   * @return Rotation matrix
   */
  [[nodiscard]] const typename SO3Type::MatrixType& R() const { return D_.R(); }

  /**
   * @brief Get a constant reference to the SE23 normalized quaternion
   *
   * @return Quaternion
   */
  [[nodiscard]] const typename SO3Type::QuaternionType& q() const { return D_.C_.q(); }

  /**
   * @brief Get a constant reference to the SE23 translation vectors (isometries)
   *
   * @return Array of R3 vectors
   */
  [[nodiscard]] const typename SE23Type::IsometriesType& t() const { return D_.t_; }

  /**
   * @brief Get a constant referece to the first isometry (velocity) of SE23
   *
   * @return R3 vector
   */
  [[nodiscard]] const typename SO3Type::VectorType& v() const { return D_.v(); }

  /**
   * @brief Get a constant referece to the second isometry (position) of SEn3 with n > 1
   *
   * @return R3 vector
   */
  [[nodiscard]] const typename SO3Type::VectorType& p() const { return D_.p(); }

  /**
   * @brief Get a constant copy of the G3 object as a matrix
   *
   * @return G3 group element in matrix form
   */
  [[nodiscard]] const MatrixType asMatrix() const {
    MatrixType X = D_.asMatrix();
    X.template block<1, 1>(3, 4) = s_;
    return X;
  }

  /**
   * @brief Set SE23 object value from given matrix
   *
   * @param T SE23 group element in matrix form
   */
  void fromT(const MatrixType& T)
  {
    D_.fromT(T);
  }

  /**
   * @brief G3 Adjoint matrix
   *
   * @return G3 group Adjoint matrix
   */
  [[nodiscard]] const TMatrixType Adjoint() const
  {
    TMatrixType Ad = TMatrixType::Identity();
    Ad.template block<9, 9>(0, 0) = D_.Adjoint();
    Ad.template block<3, 3>(6, 0) = - s_ * SO3Type::wedge(D_.v()) * D_.R();
    Ad.template block<3, 3>(6, 3) = - s_ * D_.R();
    Ad.template block<3, 1>(6, 9) = D_.v();
    return Ad;
  }

  /**
   * @brief G3 Inverse Adjoint matrix
   *
   * @return G3 group inverse Adjoint matrix
   */
  [[nodiscard]] const TMatrixType invAdjoint() const
  {
    TMatrixType Ad = TMatrixType::Identity();
    Ad.template block<9, 9>(0, 0) = D_.invAdjoint();
    Ad.template block<3, 3>(6, 3) = D_.R().transpose() * s_;
    Ad.template block<3, 1>(6, 9) = - D_.R().transpose() * D_.v();
    return Ad;
  }

 protected:
  /**
   * @brief G3 left Jacobian Q2 matrix
   *
   * @param w R3 vector
   * @param v R3 vector
   *
   * @return G3 left Jacobian Q2 matrix
   */
  [[nodiscard]] static const typename SO3Type::MatrixType G3leftJacobianQ2(const typename SO3Type::VectorType& w,
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
    FPType ang_p6 = ang_p5 * ang;
    FPType ang_p7 = ang_p6 * ang;

    FPType c1 = (0.5 * ang_p2 + c - 1.0) / ang_p4;
    FPType c2 = (c1 * ang_p3 - ang + s) / ang_p5;
    FPType c3 = -(2.0 * c + ang * s - 2.0) / ang_p4;
    FPType c4 = (c1 * ang_p3 + ang * c + ang - 2.0 * s) / ang_p5;
    FPType c5 = -(0.75 * ang * c + (0.25 * ang_p2 - 0.75) * s) / ang_p5;
    FPType c6 = (0.25 * ang * c + 0.5 * ang - 0.75 * s) / ang_p5;
    FPType c7 = ((0.25 * ang_p2 - 2.0) * c - 1.25 * ang * s + 2.0) / ang_p6;
    FPType c8 = (c1 * ang_p3 + 1.25 * ang * c + (0.25 * ang_p2 - 1.25) * s) / ang_p7;

    typename SO3Type::MatrixType m1 = r * p;
    typename SO3Type::MatrixType m2 = r * p * p;
    typename SO3Type::MatrixType m3 = p * r;
    typename SO3Type::MatrixType m4 = p * p * r;
    typename SO3Type::MatrixType m5 = p * r * p;
    typename SO3Type::MatrixType m6 = p * p * r * p;
    typename SO3Type::MatrixType m7 = p * r * p * p;
    typename SO3Type::MatrixType m8 = p * p * r * p * p;

    return 1 / 6 * r + c1 * m1 + c2 * m2 + c3 * m3 + c4 * m4 + c5 * m5 + c6 * m6 + c7 * m7 + c8 * m8;
  }

  /**
   * @brief G3 left Jacobian U1 matrix
   *
   * @param w R3 vector
   *
   * @return G3 left Jacobian U1 matrix
   */
  [[nodiscard]] static const typename SO3Type::MatrixType G3leftJacobianU1(const typename SO3Type::VectorType& u)
  {
    FPType ang = u.norm();
    if (ang < eps_)
    {
      return 0.5 * SO3Type::MatrixType::Identity() + 1 / 3 * wedge(u);
    }
    SO3Type::VectorType ax = u / ang;
    FPType ang_p2 = pow(ang, 2);
    FPType s = sin(ang);
    FPType c = cos(ang);
    FPType c1 = (ang * s + c - 1) / ang_p2;
    FPType c2 = (s - ang * c) / ang_p2;
    return c1 * SO3Type::MatrixType::Identity() + c2 * wedge(ax) + (0.5 - c1) * ax * ax.transpose();
  }


  SE23Type D_;    //!< The SE23 element of the symmetry group ancting on the extended pose
  Scalar s_;      //!< Scalar factor of the G3 group

  static constexpr FPType eps_ = 1e-6;  //!< Epsilon
};

using G3d = G3<double>;   //!< The G3 group with double precision floating point
using G3f = G3<float>;    //!< The G3 group with single precision floating point

}  // namespace group

#endif  // G3_HPP