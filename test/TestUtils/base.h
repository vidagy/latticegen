#ifndef LATTICEGEN_BASE_H
#define LATTICEGEN_BASE_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>

MATCHER_P(NearWithTolerance, tolerance, "")
{
  *result_listener << "the difference is " << std::setprecision(18) << std::get<0>(arg) - std::get<1>(arg);
  return std::fabs(std::get<0>(arg) - std::get<1>(arg)) <= tolerance;
}

#endif //LATTICEGEN_BASE_H
