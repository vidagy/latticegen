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

  
}
