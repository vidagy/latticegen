#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/Point2D.h>

using namespace Core;

TEST(TestCompareHelpers,equalsWithTolerance)
{
  const double x = 1.0;
  const double xpe = 1.0000000001;
  const double xme = 0.9999999999;

  EXPECT_TRUE( equalsWithTolerance(x,xme,1e-9,1e-9) );
  EXPECT_TRUE( equalsWithTolerance(x,xpe,1e-9,1e-9) );

  EXPECT_FALSE( equalsWithTolerance(x,xme,1e-9,1e-11) );
  EXPECT_FALSE( equalsWithTolerance(x,xpe,1e-9,1e-11) );
  EXPECT_FALSE( equalsWithTolerance(x,xme,1e-11,1e-9) );
  EXPECT_FALSE( equalsWithTolerance(x,xpe,1e-11,1e-9) );
  EXPECT_FALSE( equalsWithTolerance(x,xme,1e-11,1e-11) );
  EXPECT_FALSE( equalsWithTolerance(x,xpe,1e-11,1e-11) );

  const double y = 1000.0;
  const double ype = 1000.0000000001;
  const double yme =  999.9999999999;

  EXPECT_TRUE( equalsWithTolerance(y,yme,2e-10,2e-13) );
  EXPECT_TRUE( equalsWithTolerance(y,ype,2e-10,2e-13) );

  EXPECT_FALSE( equalsWithTolerance(y,yme,2e-11,2e-13) );
  EXPECT_FALSE( equalsWithTolerance(y,ype,2e-11,2e-13) );
  EXPECT_FALSE( equalsWithTolerance(y,yme,2e-10,2e-14) );
  EXPECT_FALSE( equalsWithTolerance(y,ype,2e-10,2e-14) );
  EXPECT_FALSE( equalsWithTolerance(y,yme,2e-11,2e-14) );
  EXPECT_FALSE( equalsWithTolerance(y,ype,2e-11,2e-14) );
}
