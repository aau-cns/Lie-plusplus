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

#ifndef TEST_GROUPS_HPP
#define TEST_GROUPS_HPP

#include <unsupported/Eigen/MatrixFunctions>

#include "groups/In.hpp"
#include "groups/SDB.hpp"
#include "groups/SOT3.hpp"
#include "groups/TG.hpp"

#include "utils/tools.hpp"

namespace test
{
using namespace group;

typedef testing::Types<SO3f, SO3d> SO3Groups;
typedef testing::Types<SOT3f, SOT3d> SOT3Groups;
typedef testing::Types<SE3d, SE3f> SE3Groups;
typedef testing::Types<SE3d, SE3f, SE23f, SE23d> SEn3Groups;
typedef testing::Types<SDBf, SDBd> SDBGroups;
typedef testing::Types<TGf, TGd> TGGroups;
typedef testing::Types<Inf, Ind> INGroups;

/**
 * @brief Tangent group specific tests
 */
template <typename T>
class TGGroupsTest : public testing::Test
{
};
TYPED_TEST_SUITE(TGGroupsTest, TGGroups);

TYPED_TEST(TGGroupsTest, TGGroupsConstructors)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    {
      auto X = TypeParam();
      MatrixEquality(X.D().asMatrix(), TypeParam::SE23Type::MatrixType::Identity());
      MatrixEquality(X.delta(), TypeParam::Vector9Type::Zero());
    }
    auto q = TypeParam::SE23Type::SO3Type::QuaternionType::UnitRandom();
    typename TypeParam::SE23Type::IsometriesType t = {TypeParam::SE23Type::SO3Type::VectorType::Random(),
                                                      TypeParam::SE23Type::SO3Type::VectorType::Random()};
    typename TypeParam::Vector9Type delta = TypeParam::Vector9Type::Random();
    {
      auto X = TypeParam(typename TypeParam::SE23Type(q, t), delta);
      MatrixEquality(X.D().asMatrix(), typename TypeParam::SE23Type(q, t).asMatrix());
      MatrixEquality(X.delta(), delta);
      auto Y = X;
      MatrixEquality(X.D().asMatrix(), Y.D().asMatrix());
      MatrixEquality(X.delta(), Y.delta());
    }
  }
}

TYPED_TEST(TGGroupsTest, TestExpLog)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    typename TypeParam::VectorType x = TypeParam::VectorType::Zero();
    typename TypeParam::VectorType y = TypeParam::log(TypeParam::exp(x));
    MatrixEquality(x, y);

    x = 1e-12 * TypeParam::VectorType::Random();
    auto X = TypeParam::exp(x);
    y = TypeParam::log(X);
    MatrixEquality(x, y);

    x = TypeParam::VectorType::Random();
    y = TypeParam::log(TypeParam::exp(x));
    MatrixEquality(x, y);

    X = TypeParam::exp(x);

    Eigen::Matrix<typename TypeParam::Scalar, 10, 10> W = Eigen::Matrix<typename TypeParam::Scalar, 10, 10>::Zero();
    W.block(0, 0, 9, 9) = SEn3<typename TypeParam::Scalar, 2>::adjoint(x.segment(0, 9));
    W.block(0, 9, 9, 1) = x.segment(9, 9);

    auto E = W.exp();

    MatrixEquality(X.D().Adjoint(), E.block(0, 0, 9, 9));
    MatrixEquality(X.delta(), E.block(0, 9, 9, 1));
  }
}

TYPED_TEST(TGGroupsTest, TestAssociativity)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto X1 = TypeParam::exp(TypeParam::VectorType::Random());
    auto X2 = TypeParam::exp(TypeParam::VectorType::Random());
    auto X3 = TypeParam::exp(TypeParam::VectorType::Random());

    auto Z1 = (X1 * X2) * X3;
    auto Z2 = X1 * (X2 * X3);

    MatrixEquality(Z1.D().asMatrix(), Z2.D().asMatrix());
    MatrixEquality(Z1.delta(), Z2.delta());
  }
}

TYPED_TEST(TGGroupsTest, TestIdentity)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto X = TypeParam::exp(TypeParam::VectorType::Random());
    auto I = TypeParam();
    typename TypeParam::SE23Type::MatrixType Imat = TypeParam::SE23Type::MatrixType::Identity();

    MatrixEquality(I.D().asMatrix(), Imat);
    MatrixEquality(I.delta(), TypeParam::Vector9Type::Zero());

    auto X1 = X * I;
    auto X2 = I * X;

    MatrixEquality(X.D().asMatrix(), X1.D().asMatrix());
    MatrixEquality(X.delta(), X1.delta());
    MatrixEquality(X.D().asMatrix(), X2.D().asMatrix());
    MatrixEquality(X.delta(), X2.delta());
  }
}

TYPED_TEST(TGGroupsTest, TestInverse)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto X = TypeParam::exp(TypeParam::VectorType::Random());
    auto X_inv = X.inv();
    auto I = TypeParam();

    auto I1 = X * X_inv;
    auto I2 = X_inv * X;

    MatrixEquality(I.D().asMatrix(), I1.D().asMatrix());
    MatrixEquality(I.delta(), I1.delta());
    MatrixEquality(I.D().asMatrix(), I2.D().asMatrix());
    MatrixEquality(I.delta(), I2.delta());
  }
}

TYPED_TEST(TGGroupsTest, TestGroupProduct)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto X = TypeParam::exp(TypeParam::VectorType::Random());
    auto Y = TypeParam::exp(TypeParam::VectorType::Random());

    auto Z = X * Y;

    MatrixEquality(Z.D().asMatrix(), X.D().asMatrix() * Y.D().asMatrix());
    MatrixEquality(Z.delta(), X.delta() + X.D().Adjoint() * Y.delta());
  }
  {
    auto X = TypeParam::exp(TypeParam::VectorType::Random());
    auto Y = TypeParam::exp(TypeParam::VectorType::Random());

    auto Z = X * Y;
    auto W = Y * X;

    auto X1 = X;
    auto X2 = X;

    X1.multiplyRight(Y);

    MatrixEquality(Z.D().asMatrix(), X1.D().asMatrix());
    MatrixEquality(Z.delta(), X1.delta());

    X2.multiplyLeft(Y);

    MatrixEquality(W.D().asMatrix(), X2.D().asMatrix());
    MatrixEquality(W.delta(), X2.delta());
  }
}

/**
 * @brief Semi Direct Bias group specific tests
 */
template <typename T>
class SDBGroupsTest : public testing::Test
{
};
TYPED_TEST_SUITE(SDBGroupsTest, SDBGroups);

TYPED_TEST(SDBGroupsTest, SDBGroupsConstructors)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    {
      auto X = TypeParam();
      MatrixEquality(X.D().asMatrix(), TypeParam::SE23Type::MatrixType::Identity());
      MatrixEquality(X.delta(), TypeParam::Vector6Type::Zero());
    }
    auto q = TypeParam::SE23Type::SO3Type::QuaternionType::UnitRandom();
    typename TypeParam::SE23Type::IsometriesType t = {TypeParam::SE23Type::SO3Type::VectorType::Random(),
                                                      TypeParam::SE23Type::SO3Type::VectorType::Random()};
    typename TypeParam::Vector6Type delta = TypeParam::Vector6Type::Random();
    {
      auto X = TypeParam(typename TypeParam::SE23Type(q, t), delta);
      MatrixEquality(X.D().asMatrix(), typename TypeParam::SE23Type(q, t).asMatrix());
      MatrixEquality(X.delta(), delta);
      auto Y = X;
      MatrixEquality(X.D().asMatrix(), Y.D().asMatrix());
      MatrixEquality(X.delta(), Y.delta());
    }
  }
}

TYPED_TEST(SDBGroupsTest, TestExpLog)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    typename TypeParam::VectorType x = TypeParam::VectorType::Zero();
    typename TypeParam::VectorType y = TypeParam::log(TypeParam::exp(x));
    MatrixEquality(x, y);

    x = 1e-12 * TypeParam::VectorType::Random();
    auto X = TypeParam::exp(x);
    y = TypeParam::log(X);
    MatrixEquality(x, y);

    x = TypeParam::VectorType::Random();
    y = TypeParam::log(TypeParam::exp(x));
    MatrixEquality(x, y);
  }
}

TYPED_TEST(SDBGroupsTest, TestAssociativity)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto X1 = TypeParam::exp(TypeParam::VectorType::Random());
    auto X2 = TypeParam::exp(TypeParam::VectorType::Random());
    auto X3 = TypeParam::exp(TypeParam::VectorType::Random());

    auto Z1 = (X1 * X2) * X3;
    auto Z2 = X1 * (X2 * X3);

    MatrixEquality(Z1.D().asMatrix(), Z2.D().asMatrix());
    MatrixEquality(Z1.delta(), Z2.delta());
  }
}

TYPED_TEST(SDBGroupsTest, TestIdentity)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto X = TypeParam::exp(TypeParam::VectorType::Random());
    auto I = TypeParam();
    typename TypeParam::SE23Type::MatrixType Imat = TypeParam::SE23Type::MatrixType::Identity();

    MatrixEquality(I.D().asMatrix(), Imat);
    MatrixEquality(I.delta(), TypeParam::Vector6Type::Zero());

    auto X1 = X * I;
    auto X2 = I * X;

    MatrixEquality(X.D().asMatrix(), X1.D().asMatrix());
    MatrixEquality(X.delta(), X1.delta());
    MatrixEquality(X.D().asMatrix(), X2.D().asMatrix());
    MatrixEquality(X.delta(), X2.delta());
  }
}

TYPED_TEST(SDBGroupsTest, TestInverse)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto X = TypeParam::exp(TypeParam::VectorType::Random());
    auto X_inv = X.inv();
    auto I = TypeParam();

    auto I1 = X * X_inv;
    auto I2 = X_inv * X;

    MatrixEquality(I.D().asMatrix(), I1.D().asMatrix());
    MatrixEquality(I.delta(), I1.delta());
    MatrixEquality(I.D().asMatrix(), I2.D().asMatrix());
    MatrixEquality(I.delta(), I2.delta());
  }
}

TYPED_TEST(SDBGroupsTest, TestGroupProduct)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto X = TypeParam::exp(TypeParam::VectorType::Random());
    auto Y = TypeParam::exp(TypeParam::VectorType::Random());

    auto Z = X * Y;

    MatrixEquality(Z.D().asMatrix(), X.D().asMatrix() * Y.D().asMatrix());
    MatrixEquality(Z.delta(), X.delta() + X.B().Adjoint() * Y.delta());
  }
  {
    auto X = TypeParam::exp(TypeParam::VectorType::Random());
    auto Y = TypeParam::exp(TypeParam::VectorType::Random());

    auto Z = X * Y;
    auto W = Y * X;

    auto X1 = X;
    auto X2 = X;

    X1.multiplyRight(Y);

    MatrixEquality(Z.D().asMatrix(), X1.D().asMatrix());
    MatrixEquality(Z.delta(), X1.delta());

    X2.multiplyLeft(Y);

    MatrixEquality(W.D().asMatrix(), X2.D().asMatrix());
    MatrixEquality(W.delta(), X2.delta());
  }
}

/**
 * @brief Intrinsic group specific tests
 */
template <typename T>
class InGroupsTest : public testing::Test
{
};
TYPED_TEST_SUITE(InGroupsTest, INGroups);

TYPED_TEST(InGroupsTest, InGroupsConstructors)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    {
      auto X = TypeParam();
      MatrixEquality(X.K(), TypeParam::MatrixType::Identity());
    }
    fp fx = utils::random<fp>(400, 800);
    fp fy = utils::random<fp>(400, 800);
    fp cx = utils::random<fp>(200, 400);
    fp cy = utils::random<fp>(200, 400);
    typename TypeParam::MatrixType L;
    L << fx, 0.0, cx, 0.0, fy, cy, 0.0, 0.0, 1.0;
    {
      auto X = TypeParam(L);
      MatrixEquality(X.K(), L);
    }
    {
      typename TypeParam::VectorType l;
      l << fx, fy, cx, cy;
      auto X = TypeParam(l);
      MatrixEquality(X.K(), L);
    }
    {
      auto X = TypeParam(fx, fy, cx, cy);
      MatrixEquality(X.K(), L);
      auto Y = X;
      MatrixEquality(X.K(), Y.K());
    }
  }
}

TYPED_TEST(InGroupsTest, InAction)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    fp fx = utils::random<fp>(400, 800);
    fp fy = utils::random<fp>(400, 800);
    fp cx = utils::random<fp>(200, 400);
    fp cy = utils::random<fp>(200, 400);
    typename TypeParam::Vector3Type x = TypeParam::Vector3Type::Random();
    auto X = TypeParam(fx, fy, cx, cy);
    MatrixEquality(X * x, X.K() * x);
  }
}

/**
 * @brief SO3 specific tests
 */
template <typename T>
class SO3GroupsTest : public testing::Test
{
};
TYPED_TEST_SUITE(SO3GroupsTest, SO3Groups);

TYPED_TEST(SO3GroupsTest, SO3Constructors)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    {
      auto X = TypeParam();
      MatrixEquality(X.R(), TypeParam::MatrixType::Identity());
      QuaternionEquality(X.q(), TypeParam::QuaternionType::Identity());
    }
    auto q = TypeParam::QuaternionType::UnitRandom();
    {
      auto X = TypeParam(q);
      MatrixEquality(X.R(), q.toRotationMatrix());
      QuaternionEquality(X.q(), q);
    }
    {
      auto X = TypeParam(q.toRotationMatrix());
      MatrixEquality(X.R(), q.toRotationMatrix());
      QuaternionEquality(X.q(), q);
      auto Y = X;
      MatrixEquality(X.R(), Y.R());
      QuaternionEquality(X.q(), Y.q());
    }
  }
}

TYPED_TEST(SO3GroupsTest, SO3Setters)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto X = TypeParam();
    auto q = TypeParam::QuaternionType::UnitRandom();
    X.fromq(q);
    MatrixEquality(X.R(), q.toRotationMatrix());
    QuaternionEquality(X.q(), q);
    q = TypeParam::QuaternionType::UnitRandom();
    X.fromR(q.toRotationMatrix());
    MatrixEquality(X.R(), q.toRotationMatrix());
    QuaternionEquality(X.q(), q);
  }
}

TYPED_TEST(SO3GroupsTest, SO3Action)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto q = TypeParam::QuaternionType::UnitRandom();
    auto X = TypeParam(q);
    typename TypeParam::VectorType x = TypeParam::VectorType::Random();
    MatrixEquality(X * x, q.toRotationMatrix() * x);
  }
}

/**
 * @brief SOT3 specific tests
 */
template <typename T>
class SOT3GroupsTest : public testing::Test
{
};
TYPED_TEST_SUITE(SOT3GroupsTest, SOT3Groups);

TYPED_TEST(SOT3GroupsTest, SOT3Constructors)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    {
      auto X = TypeParam();
      MatrixEquality(X.Q(), TypeParam::MatrixType::Identity());
      MatrixEquality(X.R(), TypeParam::SO3Type::MatrixType::Identity());
      QuaternionEquality(X.q(), TypeParam::SO3Type::QuaternionType::Identity());
      ScalarEquality(X.s(), 1.0);
    }
    auto q = TypeParam::SO3Type::QuaternionType::UnitRandom();
    fp s = utils::random<fp>(0.1, 10);
    typename TypeParam::MatrixType Q = TypeParam::MatrixType::Identity();
    Q.template block<3, 3>(0, 0) = q.toRotationMatrix();
    Q(3, 3) = s;
    {
      auto X = TypeParam(q, s);
      MatrixEquality(X.Q(), Q);
      MatrixEquality(X.R(), q.toRotationMatrix());
      QuaternionEquality(X.q(), q);
      ScalarEquality(X.s(), s);
    }
    {
      auto X = TypeParam(q.toRotationMatrix(), s);
      MatrixEquality(X.Q(), Q);
      MatrixEquality(X.R(), q.toRotationMatrix());
      QuaternionEquality(X.q(), q);
      ScalarEquality(X.s(), s);
    }
    {
      auto X = TypeParam(typename TypeParam::SO3Type(q), s);
      MatrixEquality(X.Q(), Q);
      MatrixEquality(X.R(), q.toRotationMatrix());
      QuaternionEquality(X.q(), q);
      ScalarEquality(X.s(), s);
    }
    {
      auto X = TypeParam(Q);
      MatrixEquality(X.Q(), Q);
      MatrixEquality(X.R(), q.toRotationMatrix());
      QuaternionEquality(X.q(), q);
      ScalarEquality(X.s(), s);
      auto Y = X;
      MatrixEquality(X.Q(), Y.Q());
      MatrixEquality(X.R(), Y.R());
      QuaternionEquality(X.q(), Y.q());
      ScalarEquality(X.s(), X.s());
    }
  }
}

TYPED_TEST(SOT3GroupsTest, SOT3Setters)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto X = TypeParam();
    auto q = TypeParam::SO3Type::QuaternionType::UnitRandom();
    fp s = utils::random<fp>(0.1, 10);
    typename TypeParam::MatrixType Q = TypeParam::MatrixType::Identity();
    Q.template block<3, 3>(0, 0) = q.toRotationMatrix();
    Q(3, 3) = s;
    X.fromQ(Q);
    MatrixEquality(X.Q(), Q);
    MatrixEquality(X.R(), q.toRotationMatrix());
    QuaternionEquality(X.q(), q);
    ScalarEquality(X.s(), s);
  }
}

TYPED_TEST(SOT3GroupsTest, SOT3Action)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto q = TypeParam::SO3Type::QuaternionType::UnitRandom();
    fp s = utils::random<fp>(0.1, 10);
    auto X = TypeParam(q, s);
    typename TypeParam::SO3Type::VectorType x = TypeParam::SO3Type::VectorType::Random();
    MatrixEquality(X * x, s * q.toRotationMatrix() * x);
  }
}

/**
 * @brief SEn3 specific tests
 */
template <typename T>
class SEn3GroupsTest : public testing::Test
{
};
TYPED_TEST_SUITE(SEn3GroupsTest, SEn3Groups);

TYPED_TEST(SEn3GroupsTest, SEn3Constructors)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    int n = TypeParam().t().size();
    {
      auto X = TypeParam();
      MatrixEquality(X.T(), TypeParam::MatrixType::Identity());
      MatrixEquality(X.R(), TypeParam::SO3Type::MatrixType::Identity());
      QuaternionEquality(X.q(), TypeParam::SO3Type::QuaternionType::Identity());
      for (const auto& it : X.t())
      {
        MatrixEquality(it, TypeParam::SO3Type::VectorType::Zero());
      }
    }
    auto q = TypeParam::SO3Type::QuaternionType::UnitRandom();
    typename TypeParam::IsometriesType t;
    typename TypeParam::MatrixType T = TypeParam::MatrixType::Identity();
    T.template block<3, 3>(0, 0) = q.toRotationMatrix();
    for (int i = 0; i < n; ++i)
    {
      t[i] = TypeParam::SO3Type::VectorType::Random();
      T.template block<3, 1>(0, 3 + i) = t[i];
    }
    {
      auto X = TypeParam(q, t);
      MatrixEquality(X.T(), T);
      MatrixEquality(X.R(), q.toRotationMatrix());
      QuaternionEquality(X.q(), q);
      for (int i = 0; i < n; ++i)
      {
        MatrixEquality(X.t()[i], t[i]);
      }
    }
    {
      auto X = TypeParam(q.toRotationMatrix(), t);
      MatrixEquality(X.T(), T);
      MatrixEquality(X.R(), q.toRotationMatrix());
      QuaternionEquality(X.q(), q);
      for (int i = 0; i < n; ++i)
      {
        MatrixEquality(X.t()[i], t[i]);
      }
    }
    {
      auto X = TypeParam(typename TypeParam::SO3Type(q), t);
      MatrixEquality(X.T(), T);
      MatrixEquality(X.R(), q.toRotationMatrix());
      QuaternionEquality(X.q(), q);
      for (int i = 0; i < n; ++i)
      {
        MatrixEquality(X.t()[i], t[i]);
      }
    }
    {
      auto X = TypeParam(T);
      MatrixEquality(X.T(), T);
      MatrixEquality(X.R(), q.toRotationMatrix());
      QuaternionEquality(X.q(), q);
      for (int i = 0; i < n; ++i)
      {
        MatrixEquality(X.t()[i], t[i]);
      }
      auto Y = X;
      MatrixEquality(X.T(), Y.T());
      MatrixEquality(X.R(), Y.R());
      QuaternionEquality(X.q(), Y.q());
      for (int i = 0; i < n; ++i)
      {
        MatrixEquality(X.t()[i], Y.t()[i]);
      }
    }
  }
}

TYPED_TEST(SEn3GroupsTest, SEn3Setters)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    int n = TypeParam().t().size();
    auto X = TypeParam();
    auto q = TypeParam::SO3Type::QuaternionType::UnitRandom();
    typename TypeParam::IsometriesType t;
    typename TypeParam::MatrixType T = TypeParam::MatrixType::Identity();
    T.template block<3, 3>(0, 0) = q.toRotationMatrix();
    for (int i = 0; i < n; ++i)
    {
      t[i] = TypeParam::SO3Type::VectorType::Random();
      T.template block<3, 1>(0, 3 + i) = t[i];
    }
    X.fromT(T);
    MatrixEquality(X.T(), T);
    MatrixEquality(X.R(), q.toRotationMatrix());
    QuaternionEquality(X.q(), q);
    for (int i = 0; i < n; ++i)
    {
      MatrixEquality(X.t()[i], t[i]);
    }
  }
}

/**
 * @brief SE3 specific tests
 */
template <typename T>
class SE3ActionTest : public testing::Test
{
};
TYPED_TEST_SUITE(SE3ActionTest, SE3Groups);

TYPED_TEST(SE3ActionTest, SE3Action)
{
  for (int i = 0; i < N_TESTS; ++i)
  {
    auto q = TypeParam::SO3Type::QuaternionType::UnitRandom();
    typename TypeParam::SO3Type::VectorType t = TypeParam::SO3Type::VectorType::Random();
    auto X = TypeParam(q, {t});
    typename TypeParam::SO3Type::VectorType x = TypeParam::SO3Type::VectorType::Random();
    MatrixEquality(X * x, q.toRotationMatrix() * x + t);
  }
}

}  // namespace test

#endif  // TEST_GROUPS_HPP