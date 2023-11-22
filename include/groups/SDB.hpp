// Copyright (C) 2023 Alessandro Fornasier, Pieter van Goor.
// Control of Networked Systems, University of Klagenfurt, Austria.
// System Theory and Robotics Lab, Australian Centre for Robotic
// Vision, Australian National University.
//
// All rights reserved.
//
// This software is licensed under the terms of the BSD-2-Clause-License with
// no commercial use allowed, the full terms of which are made available
// in the LICENSE file. No license in patents is granted.
//
// You can contact the authors at <alessandro.fornasier@ieee.org>

#ifndef SDB_HPP
#define SDB_HPP

#include "SEn3.hpp"

namespace group
{
/**
 * @brief The Semi Direct Bias group. This represent the core components of the symmetry group for Inertial Navigation
 * Systems (INS), including the SE23 |x se3 semi-direct product acting on the extended pose and IMU biases
 *
 * @tparam FPType. Floating point type (float, double, long double)
 *
 * @note Equivariant Symmetries for Inertial Navigation Systems [https://arxiv.org/abs/2309.03765]
 * @note MSCEqF: A Multi State Constraint Equivariant Filter for Vision-aided Inertial Navigation
 * [https://arxiv.org/abs/2311.11649]
 */
template <typename FPType>
class SemiDirectBias
{
 public:
  using Scalar = FPType;                            //!< The underlying scalar type
  using SE3Type = group::SEn3<FPType, 1>;           //!< The underlying SE3 type
  using SE23Type = group::SEn3<FPType, 2>;          //!< The underlying SE23 type
  using VectorType = Eigen::Matrix<FPType, 15, 1>;  //!< The underlying R15 vector type
  using Vector6Type = Eigen::Matrix<FPType, 6, 1>;  //!< The underlying R6 vector type

  /**
   * @brief Construct an identity Semi Direct Bias Group object
   */
  SemiDirectBias() : D_(), delta_(Vector6Type::Zero()){};

  /**
   * @brief Construct a Semi Direct Bias Group object from a given SE23 object, and a R6 vector
   *
   * @param D SE23 group element
   * @param delta R6 vector
   */
  SemiDirectBias(const SE23Type& D, const Vector6Type& delta) : D_(D), delta_(delta){};

  /**
   * @brief The exponential map for the Semi Direct Bias Group.
   * Returns a SemiDirectBias object given a vector u in R15 (equivalent to exp(wedge(u)))
   *
   * @param u R15 vector
   *
   * @return SemiDirectBias group element
   */
  [[nodiscard]] static const SemiDirectBias exp(const VectorType& u)
  {
    return SemiDirectBias(SE23Type::exp(u.template block<9, 1>(0, 0)),
                          SE3Type::leftJacobian(u.template block<6, 1>(0, 0)) * u.template block<6, 1>(9, 0));
  }

  /**
   * @brief The logarithmic map for the Semi Direct Bias Group.
   * Return a vector given a SemiDirectBias object (equivalent to vee(log(X)))
   *
   * @param X SemiDirectBias group element
   *
   * @return R15 vector
   */
  [[nodiscard]] static const VectorType log(const SemiDirectBias& X)
  {
    VectorType u = VectorType::Zero();
    u.template block<9, 1>(0, 0) = SE23Type::log(X.D_);
    u.template block<6, 1>(9, 0) = (SE3Type::leftJacobian(u.template block<6, 1>(0, 0))).inverse() * X.delta_;
    return u;
  }

  /**
   * @brief Get a constant copy of B (the SE3 subgroup of SE23 composed by the rotational component (R) and the first
   * isometry (v))
   *
   * @return SE3 group element
   */
  [[nodiscard]] const SE3Type B() const { return SE3Type(D_.q(), {D_.v()}); }

  /**
   * @brief Get a constant copy of C (the SE3 subgroup of SE23 composed by the rotational component (R) and the
   * second isometry (p))
   *
   * @return SE3 group element
   */
  [[nodiscard]] const SE3Type C() const { return SE3Type(D_.q(), {D_.p()}); }

  /**
   * @brief Get a constant reference to D (the SE23 element)
   *
   * @return SE23 group element
   */
  [[nodiscard]] const SE23Type& D() const { return D_; }

  /**
   * @brief Get a constant reference to delta (the R6 element)
   *
   * @return R6 vector
   */
  [[nodiscard]] const Vector6Type& delta() const { return delta_; }

  /**
   * @brief Get a constant copy of the inverse of the SemiDirectBias object
   *
   * @return SemiDirectBias group element
   */
  [[nodiscard]] const SemiDirectBias inv() const { return SemiDirectBias(D_.inv(), -B().invAdjoint() * delta_); }

  /**
   * @brief Operator * overloading.
   * Implements the SemiDirectBias composition this * other
   *
   * @param other SemiDirectBias group element
   *
   * @return SemiDirectBias group element
   *
   * @note usage: z = x * y
   */
  [[nodiscard]] const SemiDirectBias operator*(const SemiDirectBias& other) const
  {
    return SemiDirectBias(D_ * other.D_, delta_ + B().Adjoint() * other.delta_);
  }

  /**
   * @brief Implements the SemiDirectBias composition this = this * other
   *
   * @param other SemiDirectBias group element
   *
   * @return SemiDirectBias group element
   */
  const SemiDirectBias& multiplyRight(const SemiDirectBias& other)
  {
    delta_ = (delta_ + B().Adjoint() * other.delta_).eval();
    D_.multiplyRight(other.D_);  // D_ * other.D_
    return *this;
  }

  /**
   * @brief Implements the SemiDirectBias composition this = other * this
   *
   * @param other SemiDirectBias group element
   *
   * @return SemiDirectBias group element
   */
  const SemiDirectBias& multiplyLeft(const SemiDirectBias& other)
  {
    delta_ = (other.delta_ + other.B().Adjoint() * delta_).eval();
    D_.multiplyLeft(other.D_);  // other.D_ * th.D_
    return *this;
  }

 private:
  SE23Type D_;         //!< The SE23 element of the symmetry group ancting on the extended pose
  Vector6Type delta_;  //!< The R6 element of the symmetry group acting on the IMU biases
};

using SDBf = SemiDirectBias<float>;   //!< The Semi Direct Bias group with single precision floating point
using SDBd = SemiDirectBias<double>;  //!< The Semi Direct Bias group with double precision floating point

}  // namespace group

#endif  // SDB_HPP