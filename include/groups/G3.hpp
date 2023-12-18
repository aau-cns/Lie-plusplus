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

namespace group
{
/**
 * @brief the Inhomogeneous Galileian group (G3). This group is the group of 3D rotations,
 * translations in space and time, and transformations between frames of reference
 * that differ only by constant relative motion.
 *
 * @tparam FPType. Floating point type (float, double, long double)
 *
 * @note Galilei invariant theories [https://arxiv.org/abs/math-ph/0604002]
 * @note Constructive Equivariant Observer Design for Inertial Velocity-Aided
 * Attitude [https://arxiv.org/pdf/2209.03564.pdf]
 */
template <typename FPType>
class G3
{
 public:
  using Scalar = FPType;                              //!< The underlying scalar type
  using SO3Type = SO3<FPType>;                        //!< The underlying SO3 type
  using VectorType = Eigen::Matrix<FPType, 10, 1>;    //!< R10 Vectorspace element type (isomorphic to Lie Algebra g3)
  using MatrixType = Eigen::Matrix<FPType, 5, 5>;     //!< Lie Algebra / Lie Group matrix type
  using TMatrixType = Eigen::Matrix<FPType, 10, 10>;  //!< Transformation matrix type
  using IsometriesType = std::array<typename SO3Type::VectorType, 2>;  //!< Vector of translations (Isometries)

  /**
   * @brief Construct an identity G3 object
   */
  G3() : C_(SO3Type()), t_(), s_(0.0) { t_.fill(SO3Type::VectorType::Zero()); }

  /**
   * @brief Construct a G3 object from a given normalized quaternion, an array of vectors, and a scalar factor.
   *
   * @param q Quaternion
   * @param t Array of R3 vectors
   * @param s Scalar factor
   */
  G3(const typename SO3Type::QuaternionType& q, const IsometriesType& t, const Scalar& s) : C_(q), t_(t), s_(s) {}

  /**
   * @brief Construct a G3 object from a given rotation matrix, an array of vectors, and a scalar factor.
   *
   * @param R Rotation matrix
   * @param t Array of R3 vectors
   * @param s Scalar factor
   */
  G3(const typename SO3Type::MatrixType& R, const IsometriesType& t, const Scalar& s) : C_(R), t_(t), s_(s) {}

  /**
   * @brief Construct a G3 object from a given SO3 object, an array of vectors, and a scalar factor.
   *
   * @param C SO3 group element
   * @param t Array of R3 vectors
   * @param s Scalar factor
   */
  G3(const SO3Type& C, const IsometriesType& t, const Scalar& s) : C_(C), t_(t), s_(s) {}

  /**
   * @brief Construct a G3 object from a given matrix
   *
   * @param X G3 group element in matrix form
   */
  G3(const MatrixType& X) : C_(X.template block<3, 3>(0, 0)), t_(), s_(X(3, 4))
  {
    t_[0] = X.template block<3, 1>(0, 3);
    t_[1] = X.template block<3, 1>(0, 4);
  }

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
    U.template block<3, 3>(0, 0) = SO3Type::wedge(u.template block<3, 1>(0, 0));
    U.template block<3, 1>(0, 3) = u.template block<3, 1>(3, 0);
    U.template block<3, 1>(0, 4) = u.template block<3, 1>(6, 0);
    U(3, 4) = u(9);
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
    u.template block<3, 1>(0, 0) = SO3Type::vee(U.template block<3, 3>(0, 0));
    u.template block<3, 1>(3, 0) = U.template block<3, 1>(0, 3);
    u.template block<3, 1>(6, 0) = U.template block<3, 1>(0, 4);
    u(9) = U(3, 4);
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
    typename SO3Type::MatrixType W = SO3Type::wedge(u.template block<3, 1>(0, 0));
    ad.template block<3, 3>(0, 0) = W;
    ad.template block<3, 3>(3, 0) = SO3Type::wedge(u.template block<3, 1>(3, 0));
    ad.template block<3, 3>(3, 3) = W;
    ad.template block<3, 3>(6, 0) = SO3Type::wedge(u.template block<3, 1>(6, 0));
    ad.template block<3, 3>(6, 3) = -u(9) * SO3Type::MatrixType::Identity();
    ad.template block<3, 3>(6, 6) = W;
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
    typename SO3Type::VectorType w = u.template block<3, 1>(0, 0);
    typename SO3Type::VectorType v = u.template block<3, 1>(3, 0);
    typename SO3Type::VectorType p = u.template block<3, 1>(6, 0);
    Scalar s = u(9);
    FPType ang = w.norm();
    if (ang < eps_)
    {
      return TMatrixType::Identity() + 0.5 * adjoint(u);
    }
    typename SO3Type::MatrixType SO3JL = SO3Type::leftJacobian(w);
    TMatrixType J = TMatrixType::Identity();
    J.template block<3, 3>(0, 0) = SO3JL;
    J.template block<3, 3>(3, 0) = G3leftJacobianQ1(w, v);
    J.template block<3, 3>(3, 3) = SO3JL;
    J.template block<3, 3>(6, 0) = G3leftJacobianQ1(w, p) - s * G3leftJacobianQ2(w, v);
    J.template block<3, 3>(6, 3) = -s * G3leftJacobianU1(w);
    J.template block<3, 3>(6, 6) = SO3JL;
    J.template block<3, 1>(6, 9) = SO3Type::Gamma2(w) * v;
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
    typename SO3Type::VectorType w = u.template block<3, 1>(0, 0);
    typename SO3Type::VectorType v = u.template block<3, 1>(3, 0);
    typename SO3Type::VectorType p = u.template block<3, 1>(6, 0);
    Scalar s = u(9);
    typename SO3Type::MatrixType SO3JL = SO3Type::leftJacobian(w);
    IsometriesType t;
    t[0] = SO3JL * v;
    t[1] = SO3JL * p + s * SO3Type::Gamma2(w) * v;
    return G3(SO3Type::exp(w), t, s);
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
    u.template block<3, 1>(0, 0) = SO3Type::log(X.C_);
    typename SO3Type::MatrixType invSO3JL = SO3Type::leftJacobian(u.template block<3, 1>(0, 0)).inverse();
    u.template block<3, 1>(3, 0) = invSO3JL * X.t_[0];
    u.template block<3, 1>(6, 0) =
        invSO3JL * (X.t_[1] - X.s_ * SO3Type::Gamma2(u.template block<3, 1>(0, 0)) * u.template block<3, 1>(3, 0));
    u(9) = X.s_;

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
    IsometriesType t;
    t[0] = C_.R() * other.t_[0] + t_[0];
    t[1] = C_.R() * other.t_[1] + t_[0] * other.s_ + t_[1];
    return G3(C_ * other.C_, t, s_ + other.s_);
  }

  /**
   * @brief Operator * overloading.
   * Implements the G3 composition this * other with a G3 group or alegra element in matrix form
   *
   * @param other G3 group or algebra element in matrix form
   *
   * @return G3 group or algebra element in matrix form
   *
   * @note usage: z = x * y
   */
  [[nodiscard]] const MatrixType operator*(const MatrixType& other) const
  {
    bool is_algebra = ((other(3, 3) == 0) && (other(4, 4) == 0));
    bool is_group = ((other(3, 3) == 1) && (other(4, 4) == 1));

    if (!(is_algebra || is_group))
    {
      throw std::runtime_error(
          "G3: operator* is defined only for composition with matrix form of G3 group or algebra elements");
    }

    MatrixType res = other;
    res.template block<3, 3>(0, 0) = C_.R() * other.template block<3, 3>(0, 0);
    if (is_algebra)
    {
      res.template block<3, 1>(0, 3) = C_.R() * other.template block<3, 1>(0, 3);
      res.template block<3, 1>(0, 4) = C_.R() * other.template block<3, 1>(0, 4) + t_[0] * other(3, 4);
    }
    else
    {
      res.template block<3, 1>(0, 3) = C_.R() * other.template block<3, 1>(0, 3) + t_[0];
      res.template block<3, 1>(0, 4) = C_.R() * other.template block<3, 1>(0, 4) + t_[0] * other(3, 4) + t_[1];
      res(3, 4) += s_;
    }
    return res;
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
    t_[1] = (C_.R() * other.t_[1] + t_[0] * other.s_ + t_[1]).eval();
    t_[0] = (C_.R() * other.t_[0] + t_[0]).eval();
    s_ += other.s_;
    C_.multiplyRight(other.C_);
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
    t_[1] = (other.C_.R() * t_[1] + other.t_[0] * s_ + other.t_[1]).eval();
    t_[0] = (other.C_.R() * t_[0] + other.t_[0]).eval();
    s_ += other.s_;
    C_.multiplyLeft(other.C_);
    return *this;
  }

  /**
   * @brief Get a constant copy of the inverse of the G3 object
   *
   * @return G3 group element
   */
  [[nodiscard]] const G3 inv() const
  {
    IsometriesType t;
    t[0] = -C_.R().transpose() * t_[0];
    t[1] = -C_.R().transpose() * (t_[1] - s_ * t_[0]);
    return G3(C_.R().transpose(), t, -s_);
  }

  /**
   * @brief Get a constant reference to the SE23 rotation matrix
   *
   * @return Rotation matrix
   */
  [[nodiscard]] const typename SO3Type::MatrixType& R() const { return C_.R(); }

  /**
   * @brief Get a constant reference to the SE23 normalized quaternion
   *
   * @return Quaternion
   */
  [[nodiscard]] const typename SO3Type::QuaternionType& q() const { return C_.q(); }

  /**
   * @brief Get a constant reference to the SE23 translation vectors (isometries)
   *
   * @return Array of R3 vectors
   */
  [[nodiscard]] const IsometriesType& t() const { return t_; }

  /**
   * @brief Get a constant referece to the first isometry (velocity) of SE23
   *
   * @return R3 vector
   */
  [[nodiscard]] const typename SO3Type::VectorType& v() const { return t_[0]; }

  /**
   * @brief Get a constant referece to the second isometry (position) of SEn3 with n > 1
   *
   * @return R3 vector
   */
  [[nodiscard]] const typename SO3Type::VectorType& p() const { return t_[1]; }

  /**
   * @brief Get a constant reference to the scalar factor
   *
   * @return Scalar factor
   */
  [[nodiscard]] const Scalar& s() const { return s_; }

  /**
   * @brief Get a constant copy of the G3 object as a matrix
   *
   * @return G3 group element in matrix form
   */
  [[nodiscard]] const MatrixType asMatrix() const
  {
    MatrixType X = MatrixType::Identity();
    X.template block<3, 3>(0, 0) = C_.R();
    X.template block<3, 1>(0, 3) = t_[0];
    X.template block<3, 1>(0, 4) = t_[1];
    X(3, 4) = s_;
    return X;
  }

  /**
   * @brief Set G3 object value from given matrix
   *
   * @param X G3 group element in matrix form
   */
  void fromMatrix(const MatrixType& X)
  {
    C_.fromR(X.template block<3, 3>(0, 0));
    t_[0] = X.template block<3, 1>(0, 3);
    t_[1] = X.template block<3, 1>(0, 4);
    s_ = X(3, 4);
  }

  /**
   * @brief G3 Adjoint matrix
   *
   * @return G3 group Adjoint matrix
   */
  [[nodiscard]] const TMatrixType Adjoint() const
  {
    TMatrixType Ad = TMatrixType::Identity();
    Ad.template block<3, 3>(0, 0) = C_.R();
    Ad.template block<3, 3>(3, 0) = SO3Type::wedge(t_[0]) * C_.R();
    Ad.template block<3, 3>(3, 3) = C_.R();
    Ad.template block<3, 3>(6, 0) = SO3Type::wedge(t_[1] - s_ * t_[0]) * C_.R();
    Ad.template block<3, 3>(6, 3) = -s_ * C_.R();
    Ad.template block<3, 3>(6, 6) = C_.R();
    Ad.template block<3, 1>(6, 9) = t_[0];
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
    Ad.template block<3, 3>(0, 0) = C_.R().transpose();
    Ad.template block<3, 3>(3, 0) = -C_.R().transpose() * SO3Type::wedge(t_[0]);
    Ad.template block<3, 3>(3, 3) = C_.R().transpose();
    Ad.template block<3, 3>(6, 0) = -C_.R().transpose() * SO3Type::wedge(t_[1]);
    Ad.template block<3, 3>(6, 3) = C_.R().transpose() * s_;
    Ad.template block<3, 3>(6, 6) = C_.R().transpose();
    Ad.template block<3, 1>(6, 9) = -C_.R().transpose() * t_[0];
    return Ad;
  }

 protected:
  /**
   * @brief G3 left Jacobian Q1 matrix
   *
   * @param w R3 vector
   * @param v R3 vector
   *
   * @return G3 left Jacobian Q1 matrix
   */
  [[nodiscard]] static const typename SO3Type::MatrixType G3leftJacobianQ1(const typename SO3Type::VectorType& w,
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

    FPType c0 = 1 / 6;
    FPType c1 = (0.5 * ang_p2 + c - 1.0) / ang_p4;
    FPType c2 = (c0 * ang_p3 - ang + s) / ang_p5;
    FPType c3 = -(2.0 * c + ang * s - 2.0) / ang_p4;
    FPType c4 = (c0 * ang_p3 + ang * c + ang - 2.0 * s) / ang_p5;
    FPType c5 = -(0.75 * ang * c + (0.25 * ang_p2 - 0.75) * s) / ang_p5;
    FPType c6 = (0.25 * ang * c + 0.5 * ang - 0.75 * s) / ang_p5;
    FPType c7 = ((0.25 * ang_p2 - 2.0) * c - 1.25 * ang * s + 2.0) / ang_p6;
    FPType c8 = (c0 * ang_p3 + 1.25 * ang * c + (0.25 * ang_p2 - 1.25) * s) / ang_p7;

    typename SO3Type::MatrixType m1 = r * p;
    typename SO3Type::MatrixType m2 = r * p * p;
    typename SO3Type::MatrixType m3 = p * r;
    typename SO3Type::MatrixType m4 = p * p * r;
    typename SO3Type::MatrixType m5 = p * r * p;
    typename SO3Type::MatrixType m6 = p * p * r * p;
    typename SO3Type::MatrixType m7 = p * r * p * p;
    typename SO3Type::MatrixType m8 = p * p * r * p * p;

    return c0 * r + c1 * m1 + c2 * m2 + c3 * m3 + c4 * m4 + c5 * m5 + c6 * m6 + c7 * m7 + c8 * m8;
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
      return 0.5 * SO3Type::MatrixType::Identity() + 1 / 3 * SO3Type::wedge(u);
    }
    typename SO3Type::VectorType ax = u / ang;
    FPType ang_p2 = pow(ang, 2);
    FPType s = sin(ang);
    FPType c = cos(ang);
    FPType c1 = (ang * s + c - 1) / ang_p2;
    FPType c2 = (s - ang * c) / ang_p2;
    return c1 * SO3Type::MatrixType::Identity() + c2 * SO3Type::wedge(ax) + (0.5 - c1) * ax * ax.transpose();
  }

  SO3Type C_;         //!< The SE23 element of the symmetry group ancting on the extended pose
  IsometriesType t_;  //!< The translation vectors (isometries) of the G3 element
  Scalar s_;          //!< Scalar factor of the G3 group

  static constexpr FPType eps_ = std::is_same_v<FPType, float> ? 1.0e-6 : 1.0e-9;  //!< Epsilon
};

using G3d = G3<double>;  //!< The G3 group with double precision floating point
using G3f = G3<float>;   //!< The G3 group with single precision floating point

}  // namespace group

#endif  // G3_HPP
