// Copyright (C) 2023 Alessandro Fornasier.
// Control of Networked Systems, University of Klagenfurt, Austria.
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
// You can contact the authors at <alessandro.fornasier@ieee.org>

#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <queue>
#include <random>
#include <type_traits>
#include <Eigen/Dense>

namespace utils
{
/**
 * @brief Generate random numbers
 *
 * @tparam Numeric type of random number generated
 * @tparam Generator random number generator
 * @param from Lower bound
 * @param to Upper bound
 * @return Numeric
 *
 * @note Modified from:
 * https://stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c
 */
template <typename Numeric, typename Generator = std::mt19937>
static Numeric random(Numeric from, Numeric to)
{
  thread_local static Generator gen(std::random_device{}());
  using dist_type = typename std::conditional<std::is_integral<Numeric>::value, std::uniform_int_distribution<Numeric>,
                                              std::uniform_real_distribution<Numeric>>::type;
  thread_local static dist_type dist;
  return dist(gen, typename dist_type::param_type{from, to});
}
}  // namespace utils

#endif  // TOOLS_HPP
