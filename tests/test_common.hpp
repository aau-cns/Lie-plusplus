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

#ifndef TEST_COMMON_HPP
#define TEST_COMMON_HPP

#include <gtest/gtest.h>
#include <Eigen/Dense>

namespace test
{
using fp = double;
constexpr fp EPS = 1e-6;
constexpr int N_TESTS = 1000;

template <typename Derived, typename OtherDerived>
void MatrixEquality(const Eigen::MatrixBase<Derived>& A, const Eigen::MatrixBase<OtherDerived>& B, fp tol = EPS)
{
  EXPECT_TRUE(A.isApprox(B, tol) || (A - B).norm() < tol);
}

template <typename FPType>
void QuaternionEquality(const Eigen::Quaternion<FPType>& a, const Eigen::Quaternion<FPType>& b, fp tol = EPS)
{
  EXPECT_TRUE(a.coeffs().isApprox(b.coeffs(), tol) || a.coeffs().isApprox(-b.coeffs(), tol));
}

void ScalarEquality(const fp& a, const fp& b, fp tol = EPS) { EXPECT_TRUE(std::norm(a - b) < tol); }

}  // namespace test

#endif  // TEST_COMMON_HPP