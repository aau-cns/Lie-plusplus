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
#include "G3.hpp"

namespace group
{
/**
 * @brief The Tangent group. The tangent group of a Lie group G is a semi-direct product group between G and
 * its Lie algebra g, denoted G ⋉ g.
 *
 * @tparam GroupType. The underlying Lie group type
 *
 * @note Equivariant Filter Design for Inertial Navigation Systems with Input Measurement Biases
 * [https://arxiv.org/abs/2202.02058] [https://ieeexplore.ieee.org/document/9811778/]
 * @note Equivariant Symmetries for Inertial Navigation Systems [https://arxiv.org/abs/2309.03765]
 */
template <typename GroupType>
class Tangent
{
 public:
  using Group = GroupType;                                        //!< The underlying Lie group type
  using Scalar = typename Group::Scalar;                          //!< The underlying scalar type
  static constexpr int N = Group::VectorType::RowsAtCompileTime;  //!< The dimension of the group Lie algebra
  using VectorType = Eigen::Matrix<Scalar, 2 * N, 1>;             //!< The underlying vector type

  /**
   * @brief Construct an identity Tangent Group object
   */
  Tangent() : G_(), g_(Group::VectorType::Zero()){};

  /**
   * @brief Construct a Tangent Group object from a given group object, and a vector
   *
   * @param G Lie group element
   * @param g Lie algebra vector
   */
  Tangent(const Group& G, const typename Group::VectorType& g) : G_(G), g_(g){};

  /**
   * @brief The exponential map for the Tangent Group TG.
   * Returns a Tangent object given a VetorType u (equivalent to exp(wedge(u)))
   *
   * @param u VetorType
   *
   * @return TG group element
   *
   */
  [[nodiscard]] static const Tangent exp(const VectorType& u)
  {
    return Tangent(Group::exp(u.template block<N, 1>(0, 0)),
                   Group::leftJacobian(u.template block<N, 1>(0, 0)) * u.template block<N, 1>(N, 0));
  }

  /**
   * @brief The logarithmic map for the Tangent Group.
   * Return a vector given a Tangent object (equivalent to vee(log(X)))
   *
   * @param X TG element
   *
   * @return VectorType vector
   */
  [[nodiscard]] static const VectorType log(const Tangent& X)
  {
    VectorType u = VectorType::Zero();
    u.template block<N, 1>(0, 0) = Group::log(X.G_);
    u.template block<N, 1>(N, 0) = (Group::leftJacobian(u.template block<N, 1>(0, 0))).inverse() * X.g_;
    return u;
  }

  /**
   * @brief Get a constant reference to G (the Lie group element)
   *
   * @return Group element
   */
  [[nodiscard]] const Group& G() const { return G_; }

  /**
   * @brief Get a constant reference to g (the Lie algebra element)
   *
   * @return Lie algebra vector
   */
  [[nodiscard]] const typename Group::VectorType& g() const { return g_; }

  /**
   * @brief Get a constant copy of the inverse of the Tangent object
   *
   * @return TG group element
   */
  [[nodiscard]] const Tangent inv() const { return Tangent(G_.inv(), -G_.invAdjoint() * g_); }

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
    return Tangent(G_ * other.G_, g_ + G_.Adjoint() * other.g_);
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
    g_ = (g_ + G_.Adjoint() * other.g_).eval();
    G_.multiplyRight(other.G_);  // G_ * other.G_
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
    g_ = (other.g_ + other.G_.Adjoint() * g_).eval();
    G_.multiplyLeft(other.G_);  // other.G_ * G_
    return *this;
  }

 private:
  Group G_;                       //!< The Lie group element of the symmetry group
  typename Group::VectorType g_;  //!< The Lie algebra vector of the symmetry group
};

/**
 * @brief The SEn3 Tangent group. This derived class represents the core components of the symmetry group for Inertial
 * Navigation Systems (INS), including the SE23 |x se23 semi-direct product acting on the extended pose and IMU biases,
 * including a virtual velocity bias
 *
 * @tparam FPType. Floating point type (float, double, long double)
 * @tparam n. The dimension of the underlying vector space
 *
 */

template <typename FPType, int n>
class SEn3TG : public Tangent<SEn3<FPType, n>>
{
 public:
  using BaseType = Tangent<SEn3<FPType, n>>;
  using SE3Type = SEn3<FPType, 1>;  //!< The underlying SE3 type

  /**
   * @brief Default constructor
   */
  SEn3TG() : BaseType() {}

  /**
   * @brief Construct a SEn3TG object from a given group object and vector
   *
   * @param G Lie group element
   * @param g Lie algebra vector
   */
  SEn3TG(const typename BaseType::Group& G, const typename BaseType::Group::VectorType& g) : BaseType(G, g) {}

  /**
   * @brief Get a constant copy of B (the SE3 subgroup of SE23 composed by the rotational component (R) and the first
   * isometry (v))
   *
   * @return SE3 group element
   */
  [[nodiscard]] const SE3Type B() const { return SE3Type(this->G().q(), {this->G().v()}); }

  /**
   * @brief Get a constant copy of C (the SE3 subgroup of SE23 composed by the rotational component (R) and the
   * second isometry (p))
   *
   * @return SE3 group element
   */
  [[nodiscard]] const SE3Type C() const { return SE3Type(this->G().q(), {this->G().p()}); }
};

/**
 * @brief The G3 Tangent group. This derived class represents the core components of the symmetry group
 * for equivariant IMU preintegration
 *
 * @tparam FPType. Floating point type (float, double, long double)
 *
 */
template <typename FPType>
class G3TG : public Tangent<G3<FPType>>
{
 public:
  using BaseType = Tangent<G3<FPType>>;
  using SE3Type = SEn3<FPType, 1>;  //!< The underlying SE3 type

  /**
   * @brief Default constructor
   */
  G3TG() : BaseType() {}

  /**
   * @brief Construct a SEn3TG object from a given group object and vector
   *
   * @param G Lie group element
   * @param g Lie algebra vector
   */
  G3TG(const typename BaseType::Group& G, const typename BaseType::Group::VectorType& g) : BaseType(G, g) {}

  /**
   * @brief Get a constant copy of B (the SE3 subgroup of SE23 composed by the rotational component (R) and the first
   * isometry (v))
   *
   * @return SE3 group element
   */
  [[nodiscard]] const SE3Type B() const { return SE3Type(this->G().q(), {this->G().v()}); }

  /**
   * @brief Get a constant copy of C (the SE3 subgroup of SE23 composed by the rotational component (R) and the
   * second isometry (p))
   *
   * @return SE3 group element
   */
  [[nodiscard]] const SE3Type C() const { return SE3Type(this->G().q(), {this->G().p()}); }
};

using SE23TGd = SEn3TG<double, 2>;  //!< The SE23 tangent group with double precision floating point
using SE23TGf = SEn3TG<float, 2>;   //!< The SE23 tangent group with single precision floating point
using G3TGd = G3TG<double>;         //!< The G3 tangent group with double precision floating point
using G3TGf = G3TG<float>;          //!< The G3 tangent group with single precision floating point

}  // namespace group

#endif  // TG_HPP