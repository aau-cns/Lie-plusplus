// Copyright (C) 2023 Alessandro Fornasier.
// Control of Networked Systems, University of Klagenfurt, Austria.
// System Theory and Robotics Lab, Australian Centre for Robotic
// Vision, Australian national University, Australia.
//
// All rights reserved.
//
// This software is licensed under the terms of the Apache License, Version 2.0
// (the "License"); you may not use this file except in compliance with the
// License. You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
// License for the specific language governing permissions and limitations
// under the License.
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

constexpr fp EPS = 1e-4;
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