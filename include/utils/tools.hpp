// Copyright (C) 2023 Alessandro Fornasier.
// Control of Networked Systems, University of Klagenfurt, Austria.
//
// All rights reserved.
//
// This software is licensed under the terms of the BSD-2-Clause-License with
// no commercial use allowed, the full terms of which are made available
// in the LICENSE file. No license in patents is granted.
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
