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

#ifndef TG_HPP
#define TG_HPP

#include "SEn3.hpp"

namespace group
{
/**
 * @brief The Tangent group. This represent the core components of the symmetry group for Inertial Navigation
 * Systems (INS), including the SE23 |x se23 semi-direct product acting on the extended pose and IMU biases, including a
 * virtual velocity bias
 *
 * @tparam FPType. Floating point type (float, double, long double)
 *
 * @note Equivariant Filter Design for Inertial Navigation Systems with Input Measurement Biases
 * [https://arxiv.org/abs/2202.02058] []
 * @note Equivariant Symmetries for Inertial Navigation Systems [https://arxiv.org/abs/2309.03765]
 */
template <typename FPType>
class Tangent
{
 public:
  using SE3Type = group::SEn3<FPType, 1>;           //!< The underlying SE3 type
  using SE23Type = group::SEn3<FPType, 2>;          //!< The underlying SE23 type
  using VectorType = Eigen::Matrix<FPType, 18, 1>;  //!< The underlying R18 vector type
  using Vector9Type = Eigen::Matrix<FPType, 9, 1>;  //!< The underlying R9 vector type

  /**
   * @brief Construct an identity Tangent Group object
   */
  Tangent() : D_(), delta_(Vector9Type::Zero()){};

  /**
   * @brief Construct a Tangent Group object from a given SE23 object, and a vector
   *
   * @param D SE23 group element
   * @param delta R9 vector
   */
  Tangent(const SE23Type& D, const Vector9Type& delta) : D_(D), delta_(delta){};

  /**
   * @brief The exponential map for the Tangent Group TG.
   * Returns a Tangent object given a vector u in R18 (equivalent to exp(wedge(u)))
   *
   * @param u R18 vector
   *
   * @return TG group element
   *
   */
  [[nodiscard]] static const Tangent exp(const VectorType& u)
  {
    return Tangent(SE23Type::exp(u.template block<9, 1>(0, 0)),
                   SE23Type::leftJacobian(u.template block<9, 1>(0, 0)) * u.template block<9, 1>(9, 0));
  }

  /**
   * @brief The logarithmic map for the Tangent Group.
   * Return a vector given a Tangent object (equivalent to vee(log(X)))
   *
   * @param X TG Group element
   *
   * @return R18 vector
   */
  [[nodiscard]] static const VectorType log(const Tangent& X)
  {
    VectorType u = VectorType::Zero();
    u.template block<9, 1>(0, 0) = SE23Type::log(X.D_);
    u.template block<9, 1>(9, 0) = (SE23Type::leftJacobian(u.template block<9, 1>(0, 0))).inverse() * X.delta_;
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
   * @brief Get a constant reference to delta (the R9 element)
   *
   * @return R9 vector
   */
  [[nodiscard]] const Vector9Type& delta() const { return delta_; }

  /**
   * @brief Get a constant copy of the inverse of the Tangent object
   *
   * @return TG group element
   */
  [[nodiscard]] const Tangent inv() const { return Tangent(D_.inv(), -D_.invAdjoint() * delta_); }

  /**
   * @brief Operator * overloading.
   * Implements the Tangent composition this * other
   *
   * @param other TG group element
   *
   * @return TG group element
   *
   * @note usage: z = x * y
   */
  [[nodiscard]] const Tangent operator*(const Tangent& other) const
  {
    return Tangent(D_ * other.D_, delta_ + D_.Adjoint() * other.delta_);
  }

  /**
   * @brief Implements the Tangent composition this = this * other
   *
   * @param other TG group element
   *
   * @return TG group element
   */
  const Tangent& multiplyRight(const Tangent& other)
  {
    delta_ = (delta_ + D_.Adjoint() * other.delta_).eval();
    D_.multiplyRight(other.D_);  // D_ * other.D_
    return *this;
  }

  /**
   * @brief Implements the Tangent composition this = other * this
   *
   * @param other TG group element
   *
   * @return TG group element
   */
  const Tangent& multiplyLeft(const Tangent& other)
  {
    delta_ = (other.delta_ + other.D_.Adjoint() * delta_).eval();
    D_.multiplyLeft(other.D_);  // other.D_ * D_
    return *this;
  }

 private:
  SE23Type D_;         //!< The SE23 element of the symmetry group ancting on the extended pose
  Vector9Type delta_;  //!< The R9 element of the symmetry group acting on the IMU biases
};

using TGf = Tangent<float>;   //!< The Tangent group with single precision floating point
using TGd = Tangent<double>;  //!< The Tangent group with double precision floating point

}  // namespace group

#endif  // TG_HPP