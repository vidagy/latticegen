#ifndef LATTICEGEN_BASE_H
#define LATTICEGEN_BASE_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>

MATCHER_P(NearWithTolerance, tolerance, "")
{
  return std::fabs(std::get<0>(arg) - std::get<1>(arg)) <= tolerance;
}

#endif //LATTICEGEN_BASE_H
