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

#include <cstdlib>
#include <ctime>

#include "test_common.hpp"
#include "test_groups.hpp"

int main(int argc, char **argv)
{
  srand(static_cast<unsigned>(time(0)));
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
