#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/BravaisLattice3D.h>
#include <Core/ComparisonHelpers.h>

#include <stdexcept>
#include <tuple>
#include <math.h>

using namespace Core::Geometry;

namespace
{
  static const double pi = 3.14159265358979323846;
}

TEST(BravaisLattice3D,GetUnitCellThrows)
{
  const Point3D p000 = Point3D(0.0,0.0,0.0);
  const Point3D p001 = Point3D(0.0,0.0,1.0);
  const Point3D p010 = Point3D(0.0,1.0,0.0);
  const Point3D p100 = Point3D(1.0,0.0,0.0);
  const Point3D p011 = Point3D(0.0,1.0,1.0);
  const Point3D p101 = Point3D(1.0,0.0,1.0);
  const Point3D p110 = Point3D(1.0,1.0,0.0);
  const Point3D p111 = Point3D(1.0,1.0,1.0);
/*
  EXPECT_THROW(BravaisLattice3D::get_unit_cell(p000, p011, p100), std::invalid_argument);
  EXPECT_THROW(BravaisLattice3D::get_unit_cell(p011, p000, p100), std::invalid_argument);
  EXPECT_THROW(BravaisLattice3D::get_unit_cell(p011, p100, p000), std::invalid_argument);

  EXPECT_THROW(BravaisLattice3D::get_unit_cell(p100, p100, p111), std::invalid_argument);
  EXPECT_THROW(BravaisLattice3D::get_unit_cell(p100, p111, p100), std::invalid_argument);
  EXPECT_THROW(BravaisLattice3D::get_unit_cell(p111, p100, p100), std::invalid_argument);
  
  EXPECT_THROW(BravaisLattice3D::get_unit_cell(p000, p000, p111), std::invalid_argument);
  EXPECT_THROW(BravaisLattice3D::get_unit_cell(p000, p111, p000), std::invalid_argument);
  EXPECT_THROW(BravaisLattice3D::get_unit_cell(p111, p000, p000), std::invalid_argument);
*/
  EXPECT_NO_THROW(BravaisLattice3D::get_unit_cell(p001, p010, p100));
  EXPECT_NO_THROW(BravaisLattice3D::get_unit_cell(p010, p100, p001));
  EXPECT_NO_THROW(BravaisLattice3D::get_unit_cell(p100, p001, p010));
}


