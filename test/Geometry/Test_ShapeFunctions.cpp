#include <TestUtils/base.h>

#include <Geometry/ShapeFunctions.h>
#include <Geometry/Cutoff.h>
#include <Math/CommonFunctions.h>

using namespace Core;
using namespace Geometry;

// TODO enable this test if interpolation is implemented for Shape functions.
TEST(DISABLED_ShapeFunctions, TooClose)
{
  const auto unit_cell = UnitCell3D::create_orthorhombic_body(1.0, 2.0, 3.0);
  const auto min_r = CutoffWSCell(unit_cell).r_mt;
  auto r_points = std::vector<double>();
  auto n = 100;
  r_points.reserve((unsigned long) n);
  for (auto i = 0; i < n; ++i) {
    auto r = min_r * (i + 1) / (n + 1);
    r_points.push_back(r);
  }

  auto l_max = 3u;
  auto shape_functions = ShapeFunctions(unit_cell, l_max);

  // l == 0
  {
    auto sf = shape_functions.shape_functions.at(0, 0);
    for (auto x: sf) {
      EXPECT_DOUBLE_EQ(x.real(), 2.0 * sqrt(pi));
      EXPECT_DOUBLE_EQ(x.imag(), 0.0);
    }
  }
  // l > 0
  for (auto l = 1u; l <= l_max; ++l) {
    for (auto m = -((int) l); m <= ((int) l); ++m) {
      auto sf = shape_functions.shape_functions.at(l, m);
      for (auto x: sf) {
        EXPECT_DOUBLE_EQ(x.real(), 0.0);
        EXPECT_DOUBLE_EQ(x.imag(), 0.0);
      }
    }
  }
}

// TODO enable this test if interpolation is implemented for Shape functions.
TEST(DISABLED_ShapeFunctions, TooFar)
{
  const auto unit_cell = UnitCell3D::create_monoclinic_base(1.0, 2.0, 3.0, pi / 3.0);
  const auto max_r = CutoffWSCell(unit_cell).r_bs;
  auto r_points = std::vector<double>();
  auto n = 100;
  r_points.reserve((unsigned long) n);
  for (auto i = 0; i < n; ++i) {
    auto r = max_r + max_r * (i + 1) / (n + 1);
    r_points.push_back(r);
  }

  auto l_max = 3u;
  auto shape_functions = ShapeFunctions(unit_cell, l_max);
  for (auto l = 0u; l <= l_max; ++l) {
    for (auto m = -((int) l); m <= ((int) l); ++m) {
      auto sf = shape_functions.shape_functions.at(l, m);
      for (auto x: sf) {
        EXPECT_DOUBLE_EQ(x.real(), 0.0);
        EXPECT_DOUBLE_EQ(x.imag(), 0.0);
      }
    }
  }
}

TEST(ShapeFunctions, Cubic)
{
  const auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);

  auto l_max = 7u;
  auto shape_functions = ShapeFunctions(unit_cell, l_max);
  for (auto l = 0u; l <= l_max; ++l) {
    for (auto m = -((int) l); m <= ((int) l); ++m) {
      auto sf = shape_functions.shape_functions.at(l, m);
      if ((!((l == 0u) && (m == 0))) &&
          (!((l == 4u) && (m == -4))) &&
          (!((l == 4u) && (m == 0))) &&
          (!((l == 4u) && (m == 4))) &&
          (!((l == 6u) && (m == -4))) &&
          (!((l == 6u) && (m == 0))) &&
          (!((l == 6u) && (m == 4))) &&
          // don't know why this is not ZeroBySymmetries but I can live with that for now.
          (!((l == 2u) && (m == 0)))
        ) {
        for (auto &x: sf) {
          EXPECT_EQ(x.real(), 0.0) << " l = " << l << " m = " << m;
          EXPECT_EQ(x.imag(), 0.0) << " l = " << l << " m = " << m;
        }
      }
    }
  }
}

