#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/ComparisonHelpers.h>

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

TEST(TestCompareHelpers,NearlyZero)
{
  const double pe =  0.0000000001;
  const double me = -0.0000000001;

  EXPECT_TRUE( nearlyZero(me,1e-9) );
  EXPECT_TRUE( nearlyZero(pe,1e-9) );

  EXPECT_FALSE( nearlyZero(me,1e-11) );
  EXPECT_FALSE( nearlyZero(pe,1e-11) );
}

TEST(TestCompareHelpers,StrictlyGreater)
{
  const double x = 1.0;
  const double pe =  0.0000000001;
  const double me = -0.0000000001;

  EXPECT_FALSE( strictlyGreater(x, x + me,1e-9) );
  EXPECT_FALSE( strictlyGreater(x, x + pe,1e-9) );

  EXPECT_TRUE ( strictlyGreater(x, x + me,1e-11) );
  EXPECT_FALSE( strictlyGreater(x, x + pe,1e-11) );
}

TEST(TestCompareHelpers,StrictlyLess)
{
  const double x = 1.0;
  const double pe =  0.0000000001;
  const double me = -0.0000000001;

  EXPECT_FALSE( strictlyLess(x, x + me,1e-9) );
  EXPECT_FALSE( strictlyLess(x, x + pe,1e-9) );

  EXPECT_FALSE( strictlyLess(x, x + me,1e-11) );
  EXPECT_TRUE ( strictlyLess(x, x + pe,1e-11) );
}

TEST(TestCompareHelpers,GreaterEqualsWithTolerance)
{
  const double x = 1.0;
  const double pe =  0.0000000001;
  const double me = -0.0000000001;

  EXPECT_TRUE( greaterEqualsWithTolerance(x, x + me,1e-9) );
  EXPECT_TRUE( greaterEqualsWithTolerance(x, x + pe,1e-9) );

  EXPECT_TRUE( greaterEqualsWithTolerance(x, x + me,1e-11) );
  EXPECT_FALSE( greaterEqualsWithTolerance(x, x + pe,1e-11) );
}

TEST(TestCompareHelpers,LessEqualsWithTolerance)
{
  const double x = 1.0;
  const double pe =  0.0000000001;
  const double me = -0.0000000001;

  EXPECT_TRUE( lessEqualsWithTolerance(x, x + me,1e-9) );
  EXPECT_TRUE( lessEqualsWithTolerance(x, x + pe,1e-9) );

  EXPECT_FALSE( lessEqualsWithTolerance(x, x + me,1e-11) );
  EXPECT_TRUE( lessEqualsWithTolerance(x, x + pe,1e-11) );
}

TEST(TestCompareHelpers,PositiveWithTolerance)
{
  const double pe =  0.0000000001;
  const double me = -0.0000000001;

  EXPECT_TRUE( positiveWithTolerance(me,1e-9) );
  EXPECT_TRUE( positiveWithTolerance(pe,1e-9) );

  EXPECT_FALSE( positiveWithTolerance(me,1e-11) );
  EXPECT_TRUE( positiveWithTolerance(pe,1e-11) );
}

TEST(TestCompareHelpers,NegativeWithTolerance)
{
  const double pe =  0.0000000001;
  const double me = -0.0000000001;

  EXPECT_TRUE( negativeWithTolerance(me,1e-9) );
  EXPECT_TRUE( negativeWithTolerance(pe,1e-9) );

  EXPECT_TRUE( negativeWithTolerance(me,1e-11) );
  EXPECT_FALSE( negativeWithTolerance(pe,1e-11) );
}

TEST(TestCompareHelpers,StrictlyPositive)
{
  const double pe =  0.0000000001;
  const double me = -0.0000000001;

  EXPECT_FALSE( strictlyPositive(me,1e-9) );
  EXPECT_FALSE( strictlyPositive(pe,1e-9) );

  EXPECT_FALSE( strictlyPositive(me,1e-11) );
  EXPECT_TRUE( strictlyPositive(pe,1e-11) );
}

TEST(TestCompareHelpers,StrictlyNegative)
{
  const double pe =  0.0000000001;
  const double me = -0.0000000001;

  EXPECT_FALSE( strictlyNegative(me,1e-9) );
  EXPECT_FALSE( strictlyNegative(pe,1e-9) );

  EXPECT_TRUE( strictlyNegative(me,1e-11) );
  EXPECT_FALSE( strictlyNegative(pe,1e-11) );
}
