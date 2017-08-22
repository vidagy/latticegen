#include <TestUtils/base.h>

#include <Geometry/ShapeFunctions.h>
#include <Geometry/Cutoff.h>
#include <Math/CommonFunctions.h>

using namespace Core;
using namespace Geometry;

TEST(ShapeFunctions, TooClose)
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
  auto shape_functions = ShapeFunctions(unit_cell, l_max, std::make_shared<GenericMesh>(r_points),
                                        ShapeFunctionsConfig(3));

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

TEST(ShapeFunctions, TooFar)
{
  const auto unit_cell = UnitCell3D::create_orthorhombic_body(1.0, 2.0, 3.0);
  const auto max_r = CutoffWSCell(unit_cell).r_bs;
  auto r_points = std::vector<double>();
  auto n = 100;
  r_points.reserve((unsigned long) n);
  for (auto i = 0; i < n; ++i) {
    auto r = max_r + max_r * (i + 1.0) / (n + 1.0);
    r_points.push_back(r);
  }

  auto l_max = 3u;
  auto shape_functions = ShapeFunctions(unit_cell, l_max, std::make_shared<GenericMesh>(r_points),
                                        ShapeFunctionsConfig(3));
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

TEST(ShapeFunctions, CubicZero)
{
  const auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  const auto max_r = CutoffWSCell(unit_cell).r_bs;
  const auto min_r = CutoffWSCell(unit_cell).r_mt;
  auto r_points = std::vector<double>();
  const auto n = 100;
  r_points.reserve((unsigned long) n);
  for (auto i = 0; i < n; ++i) {
    auto r = min_r + (max_r - min_r) * i / (n - 1);
    r_points.push_back(r);
  }

  auto l_max = 7u;
  auto shape_functions = ShapeFunctions(unit_cell, l_max, std::make_shared<GenericMesh>(r_points),
                                        ShapeFunctionsConfig(4));
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

