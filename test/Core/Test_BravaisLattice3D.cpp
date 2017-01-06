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

TEST(BravaisLattice3DUnitCell,CreateThrowsForZeroLength)
{
  const double a_ = 1.0;
  const double b_ = 2.0;
  const double c_ = 3.0;
  
  const double alpha_ = pi/2.78;
  const double beta_ = pi/3.18;
  const double gamma_ = pi/3.78;

  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(0.0,  b_,  c_,  alpha_,  beta_,  gamma_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  0.0,  c_,  alpha_,  beta_,  gamma_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  b_,  0.0,  alpha_,  beta_,  gamma_), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_primitive(0.0,  b_,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_primitive(a_,  0.0,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_primitive(a_,  b_,  0.0,  beta_), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_base(0.0,  b_,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_base(a_,  0.0,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_base(a_,  b_,  0.0,  beta_), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_primitive(0.0,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_primitive(a_,  0.0,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_primitive(a_,  b_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_base(0.0,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_base(a_,  0.0,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_base(a_,  b_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_body(0.0,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_body(a_,  0.0,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_body(a_,  b_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_face(0.0,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_face(a_,  0.0,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_face(a_,  b_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_tetragonal_primitive(0.0,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_tetragonal_primitive(a_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_tetragonal_body(0.0,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_tetragonal_body(a_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_rhombohedral(0.0,  alpha_), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_hexagonal(0.0,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_hexagonal(a_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_cubic_primitive(0.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_cubic_body(0.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_cubic_face(0.0), std::invalid_argument);
}

TEST(BravaisLattice3DUnitCell,CreateThrowsForNegativeLength)
{
  const double a_ = 1.0;
  const double b_ = 2.0;
  const double c_ = 3.0;
  
  const double alpha_ = pi/2.78;
  const double beta_ = pi/3.18;
  const double gamma_ = pi/3.78;

  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(-1.5,  b_,  c_,  alpha_,  beta_,  gamma_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  -1.5,  c_,  alpha_,  beta_,  gamma_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  b_,  -1.5,  alpha_,  beta_,  gamma_), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_primitive(-1.5,  b_,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_primitive(a_,  -1.5,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_primitive(a_,  b_,  -1.5,  beta_), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_base(-1.5,  b_,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_base(a_,  -1.5,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_base(a_,  b_,  -1.5,  beta_), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_primitive(-1.5,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_primitive(a_,  -1.5,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_primitive(a_,  b_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_base(-1.5,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_base(a_,  -1.5,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_base(a_,  b_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_body(-1.5,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_body(a_,  -1.5,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_body(a_,  b_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_face(-1.5,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_face(a_,  -1.5,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_orthorhombic_face(a_,  b_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_tetragonal_primitive(-1.5,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_tetragonal_primitive(a_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_tetragonal_body(-1.5,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_tetragonal_body(a_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_rhombohedral(-1.5,  alpha_), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_hexagonal(-1.5,  c_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_hexagonal(a_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_cubic_primitive(-1.5), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_cubic_body(-1.5), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_cubic_face(-1.5), std::invalid_argument);
}

TEST(BravaisLattice3DUnitCell,CreateThrowsForZeroAngle)
{
  const double a_ = 1.0;
  const double b_ = 2.0;
  const double c_ = 3.0;
  
  const double alpha_ = pi/2.78;
  const double beta_ = pi/3.18;
  const double gamma_ = pi/3.78;

  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  b_,  c_,     0.0,  beta_,  gamma_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  b_,  c_,  alpha_,    0.0,  gamma_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  b_,  c_,  alpha_,  beta_,     0.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_primitive(a_,  b_,  c_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_base(a_,  b_,  c_,  0.0), std::invalid_argument);

  EXPECT_THROW( BravaisLattice3D::UnitCell::create_rhombohedral(a_,  0.0), std::invalid_argument);
}

TEST(BravaisLattice3DUnitCell,CreateThrowsForNegativeAngle)
{
  const double a_ = 1.0;
  const double b_ = 2.0;
  const double c_ = 3.0;
  
  const double alpha_ = pi/2.78;
  const double beta_ = pi/3.18;
  const double gamma_ = pi/3.78;

  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  b_,  c_, -pi/4.0,  beta_,  gamma_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  b_,  c_,  alpha_,-pi/4.0,  gamma_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  b_,  c_,  alpha_,  beta_, -pi/4.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_primitive(a_,  b_,  c_, -pi/4.0), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_base(a_,  b_,  c_, -pi/4.0), std::invalid_argument);

  EXPECT_THROW( BravaisLattice3D::UnitCell::create_rhombohedral(a_, -pi/4.0), std::invalid_argument);
}

TEST(BravaisLattice3DUnitCell,CreateThrowsForObtuseAngle)
{
  const double a_ = 1.0;
  const double b_ = 2.0;
  const double c_ = 3.0;
  
  const double alpha_ = pi/2.78;
  const double beta_ = pi/3.18;
  const double gamma_ = pi/3.78;

  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  b_,  c_,  pi/1.5,  beta_,  gamma_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  b_,  c_,  alpha_, pi/1.5,  gamma_), std::invalid_argument);
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_triclinic(a_,  b_,  c_,  alpha_,  beta_,  pi/1.5), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_primitive(a_,  b_,  c_, pi/1.5), std::invalid_argument);
  
  EXPECT_THROW( BravaisLattice3D::UnitCell::create_monoclinic_base(a_,  b_,  c_, pi/1.5), std::invalid_argument);

  EXPECT_THROW( BravaisLattice3D::UnitCell::create_rhombohedral(a_, pi/1.5), std::invalid_argument);
}
