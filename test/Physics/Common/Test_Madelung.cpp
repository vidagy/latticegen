#include <TestUtils/base.h>

#include <Physics/Common/MadelungConstants.h>

using namespace Physics::Common;
using namespace std::complex_literals;

TEST(TestRealMadelungConstants, Calculate)
{
  auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto madelung_constant = RealMadelungConstants(unit_cell);
  auto nothing = 0.0i;
  for (auto l = 1u; l <= 3u; ++l) {
    for (auto m = -((int) l); m <= ((int) l); ++m) {
      for (auto lp = 1u; lp <= 3u; ++lp) {
        for (auto mp = -((int) lp); mp <= ((int) lp); ++mp) {
          for (auto a = 0; a < 4; ++a) {
            for (auto b = 0; b < 4; ++b) {
              for (auto c = 0; c < 4; ++c) {
                if (!(a == 0 && b == 0 && c == 0)) {
                  auto res = madelung_constant.calculate(l, m, lp, mp, Coordinates3D(0, 0, 0), Coordinates3D(a, b, c));
                  nothing += res;
                }
              }
            }
          }
        }
      }
    }
  }
}

TEST(TestRealMadelungConstants, CalculateReduced)
{
  auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto madelung_constant = RealMadelungConstants(unit_cell);
  auto nothing = 0.0i;
  for (auto l = 1u; l <= 3u; ++l) {
    for (auto m = -((int) l); m <= ((int) l); ++m) {
      for (auto a = 0; a < 4; ++a) {
        for (auto b = 0; b < 4; ++b) {
          for (auto c = 0; c < 4; ++c) {
            if (!(a == 0 && b == 0 && c == 0)) {
              auto res = madelung_constant.calculateReduced(l, m, Coordinates3D(0, 0, 0), Coordinates3D(a, b, c));
              nothing += res;
            }
          }
        }
      }
    }
  }
}

TEST(TestReducedMadelungConstants, Calculate)
{
  auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto madelung_calculator = ReducedMadelungConstants(unit_cell, 3u);
  const auto &constants = madelung_calculator.madelung_constants;
  auto nothing = 0.0i;
  for (auto l = 0u; l <= 3u; ++l) {
    for (auto m = -((int) l); m <= ((int) l); ++m) {
      nothing += constants.at(l, m);
    }
  }
}

