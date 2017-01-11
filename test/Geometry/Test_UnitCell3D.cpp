#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Geometry/UnitCell3D.h>
#include <Core/ComparisonHelpers.h>

#include <stdexcept>
#include <iostream>
#include <tuple>
#include <math.h>

using namespace Core;
using namespace Geometry;

namespace
{
  static const double pi = 3.14159265358979323846;
}

TEST(UnitCell3D,CreateThrowsForZeroLength)
{
  const double a_ = 1.0;
  const double b_ = 2.0;
  const double c_ = 3.0;
  
  const double alpha_ = pi/2.78;
  const double beta_ = pi/3.18;
  const double gamma_ = pi/3.78;

  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(0.0, b_, c_, alpha_, beta_, gamma_), std::invalid_argument);
  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, 0.0, c_, alpha_, beta_, gamma_), std::invalid_argument);
  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, b_, 0.0, alpha_, beta_, gamma_), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_monoclinic_primitive(0.0,  b_,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_monoclinic_primitive(a_,  0.0,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_monoclinic_primitive(a_,  b_,  0.0,  beta_), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_monoclinic_base(0.0,  b_,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_monoclinic_base(a_,  0.0,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_monoclinic_base(a_,  b_,  0.0,  beta_), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_orthorhombic_primitive(0.0,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_primitive(a_,  0.0,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_primitive(a_,  b_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_orthorhombic_base(0.0,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_base(a_,  0.0,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_base(a_,  b_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_orthorhombic_body(0.0,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_body(a_,  0.0,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_body(a_,  b_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_orthorhombic_face(0.0,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_face(a_,  0.0,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_face(a_,  b_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_tetragonal_primitive(0.0,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_tetragonal_primitive(a_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_tetragonal_body(0.0,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_tetragonal_body(a_,  0.0), std::invalid_argument);
  
  EXPECT_THROW(UnitCell3D::create_rhombohedral_centered(0.0, alpha_), std::invalid_argument);
  
  EXPECT_THROW(UnitCell3D::create_hexagonal_primitive(0.0, c_), std::invalid_argument);
  EXPECT_THROW(UnitCell3D::create_hexagonal_primitive(a_, 0.0), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_cubic_primitive(0.0), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_cubic_body(0.0), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_cubic_face(0.0), std::invalid_argument);
}

TEST(UnitCell3D,CreateThrowsForNegativeLength)
{
  const double a_ = 1.0;
  const double b_ = 2.0;
  const double c_ = 3.0;
  
  const double alpha_ = pi/2.78;
  const double beta_ = pi/3.18;
  const double gamma_ = pi/3.78;

  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(-1.5, b_, c_, alpha_, beta_, gamma_), std::invalid_argument);
  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, -1.5, c_, alpha_, beta_, gamma_), std::invalid_argument);
  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, b_, -1.5, alpha_, beta_, gamma_), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_monoclinic_primitive(-1.5,  b_,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_monoclinic_primitive(a_,  -1.5,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_monoclinic_primitive(a_,  b_,  -1.5,  beta_), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_monoclinic_base(-1.5,  b_,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_monoclinic_base(a_,  -1.5,  c_,  beta_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_monoclinic_base(a_,  b_,  -1.5,  beta_), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_orthorhombic_primitive(-1.5,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_primitive(a_,  -1.5,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_primitive(a_,  b_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_orthorhombic_base(-1.5,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_base(a_,  -1.5,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_base(a_,  b_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_orthorhombic_body(-1.5,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_body(a_,  -1.5,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_body(a_,  b_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_orthorhombic_face(-1.5,  b_,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_face(a_,  -1.5,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_orthorhombic_face(a_,  b_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_tetragonal_primitive(-1.5,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_tetragonal_primitive(a_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_tetragonal_body(-1.5,  c_), std::invalid_argument);
  EXPECT_THROW( UnitCell3D::create_tetragonal_body(a_,  -1.5), std::invalid_argument);
  
  EXPECT_THROW(UnitCell3D::create_rhombohedral_centered(-1.5, alpha_), std::invalid_argument);
  
  EXPECT_THROW(UnitCell3D::create_hexagonal_primitive(-1.5, c_), std::invalid_argument);
  EXPECT_THROW(UnitCell3D::create_hexagonal_primitive(a_, -1.5), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_cubic_primitive(-1.5), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_cubic_body(-1.5), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_cubic_face(-1.5), std::invalid_argument);
}

TEST(UnitCell3D,CreateThrowsForZeroAngle)
{
  const double a_ = 1.0;
  const double b_ = 2.0;
  const double c_ = 3.0;
  
  const double alpha_ = pi/2.78;
  const double beta_ = pi/3.18;
  const double gamma_ = pi/3.78;

  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, b_, c_, 0.0, beta_, gamma_), std::invalid_argument);
  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, b_, c_, alpha_, 0.0, gamma_), std::invalid_argument);
  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, b_, c_, alpha_, beta_, 0.0), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_monoclinic_primitive(a_,  b_,  c_,  0.0), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_monoclinic_base(a_,  b_,  c_,  0.0), std::invalid_argument);

  EXPECT_THROW(UnitCell3D::create_rhombohedral_centered(a_, 0.0), std::invalid_argument);
}

TEST(UnitCell3D,CreateThrowsForNegativeAngle)
{
  const double a_ = 1.0;
  const double b_ = 2.0;
  const double c_ = 3.0;
  
  const double alpha_ = pi/2.78;
  const double beta_ = pi/3.18;
  const double gamma_ = pi/3.78;

  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, b_, c_, -pi / 4.0, beta_, gamma_), std::invalid_argument);
  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, b_, c_, alpha_, -pi / 4.0, gamma_), std::invalid_argument);
  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, b_, c_, alpha_, beta_, -pi / 4.0), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_monoclinic_primitive(a_,  b_,  c_, -pi/4.0), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_monoclinic_base(a_,  b_,  c_, -pi/4.0), std::invalid_argument);

  EXPECT_THROW(UnitCell3D::create_rhombohedral_centered(a_, -pi / 4.0), std::invalid_argument);
}

TEST(UnitCell3D,CreateThrowsForObtuseAngle)
{
  const double a_ = 1.0;
  const double b_ = 2.0;
  const double c_ = 3.0;
  
  const double alpha_ = pi/2.78;
  const double beta_ = pi/3.18;
  const double gamma_ = pi/3.78;

  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, b_, c_, pi / 1.5, beta_, gamma_), std::invalid_argument);
  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, b_, c_, alpha_, pi / 1.5, gamma_), std::invalid_argument);
  EXPECT_THROW(UnitCell3D::create_triclinic_primitive(a_, b_, c_, alpha_, beta_, pi / 1.5), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_monoclinic_primitive(a_,  b_,  c_, pi/1.5), std::invalid_argument);
  
  EXPECT_THROW( UnitCell3D::create_monoclinic_base(a_,  b_,  c_, pi/1.5), std::invalid_argument);

  EXPECT_THROW(UnitCell3D::create_rhombohedral_centered(a_, pi / 1.5), std::invalid_argument);
}

namespace
{
  void test_get_offset(const UnitCell3D& unit_cell)
  {
    const Point3D& a = unit_cell.a;
    const Point3D& b = unit_cell.b;
    const Point3D& c = unit_cell.c;

    EXPECT_EQ( unit_cell.get_offsets( ( 0.0)*a + ( 0.0)*b + ( 0.0)*c ), std::make_tuple( 0, 0, 0) );
    
    EXPECT_EQ( unit_cell.get_offsets( ( 1.0)*a + ( 0.0)*b + ( 0.0)*c ), std::make_tuple( 1, 0, 0) );
    EXPECT_EQ( unit_cell.get_offsets( (-1.0)*a + ( 0.0)*b + ( 0.0)*c ), std::make_tuple(-1, 0, 0) );
    EXPECT_EQ( unit_cell.get_offsets( ( 0.0)*a + ( 1.0)*b + ( 0.0)*c ), std::make_tuple( 0, 1, 0) );
    EXPECT_EQ( unit_cell.get_offsets( ( 0.0)*a + (-1.0)*b + ( 0.0)*c ), std::make_tuple( 0,-1, 0) );
    EXPECT_EQ( unit_cell.get_offsets( ( 0.0)*a + ( 0.0)*b + ( 1.0)*c ), std::make_tuple( 0, 0, 1) );
    EXPECT_EQ( unit_cell.get_offsets( ( 0.0)*a + ( 0.0)*b + (-1.0)*c ), std::make_tuple( 0, 0,-1) );

    EXPECT_EQ( unit_cell.get_offsets( ( 1.0)*a + ( 1.0)*b + ( 0.0)*c ), std::make_tuple( 1, 1, 0) );
    EXPECT_EQ( unit_cell.get_offsets( (-1.0)*a + ( 1.0)*b + ( 0.0)*c ), std::make_tuple(-1, 1, 0) );
    EXPECT_EQ( unit_cell.get_offsets( ( 1.0)*a + (-1.0)*b + ( 0.0)*c ), std::make_tuple( 1,-1, 0) );
    EXPECT_EQ( unit_cell.get_offsets( (-1.0)*a + (-1.0)*b + ( 0.0)*c ), std::make_tuple(-1,-1, 0) );
    
    EXPECT_EQ( unit_cell.get_offsets( ( 1.0)*a + ( 0.0)*b + ( 1.0)*c ), std::make_tuple( 1, 0, 1) );
    EXPECT_EQ( unit_cell.get_offsets( (-1.0)*a + ( 0.0)*b + ( 1.0)*c ), std::make_tuple(-1, 0, 1) );
    EXPECT_EQ( unit_cell.get_offsets( ( 1.0)*a + ( 0.0)*b + (-1.0)*c ), std::make_tuple( 1, 0,-1) );
    EXPECT_EQ( unit_cell.get_offsets( (-1.0)*a + ( 0.0)*b + (-1.0)*c ), std::make_tuple(-1, 0,-1) );

    EXPECT_EQ( unit_cell.get_offsets( ( 0.0)*a + ( 1.0)*b + ( 1.0)*c ), std::make_tuple( 0, 1, 1) );
    EXPECT_EQ( unit_cell.get_offsets( ( 0.0)*a + (-1.0)*b + ( 1.0)*c ), std::make_tuple( 0,-1, 1) );
    EXPECT_EQ( unit_cell.get_offsets( ( 0.0)*a + ( 1.0)*b + (-1.0)*c ), std::make_tuple( 0, 1,-1) );
    EXPECT_EQ( unit_cell.get_offsets( ( 0.0)*a + (-1.0)*b + (-1.0)*c ), std::make_tuple( 0,-1,-1) );


    EXPECT_EQ( unit_cell.get_offsets( ( 1.0)*a + ( 1.0)*b + ( 1.0)*c ), std::make_tuple( 1, 1, 1) );
    EXPECT_EQ( unit_cell.get_offsets( (-1.0)*a + ( 1.0)*b + ( 1.0)*c ), std::make_tuple(-1, 1, 1) );
    EXPECT_EQ( unit_cell.get_offsets( ( 1.0)*a + (-1.0)*b + ( 1.0)*c ), std::make_tuple( 1,-1, 1) );
    EXPECT_EQ( unit_cell.get_offsets( ( 1.0)*a + ( 1.0)*b + (-1.0)*c ), std::make_tuple( 1, 1,-1) );
    EXPECT_EQ( unit_cell.get_offsets( (-1.0)*a + (-1.0)*b + ( 1.0)*c ), std::make_tuple(-1,-1, 1) );
    EXPECT_EQ( unit_cell.get_offsets( (-1.0)*a + ( 1.0)*b + (-1.0)*c ), std::make_tuple(-1, 1,-1) );
    EXPECT_EQ( unit_cell.get_offsets( ( 1.0)*a + (-1.0)*b + (-1.0)*c ), std::make_tuple( 1,-1,-1) );
    EXPECT_EQ( unit_cell.get_offsets( (-1.0)*a + (-1.0)*b + (-1.0)*c ), std::make_tuple(-1,-1,-1) );
  }
}

TEST(UnitCell3D,GetOffset)
{
  const double a_ = 1.0;
  const double b_ = 2.0;
  const double c_ = 3.0;
  
  const double alpha_ = pi/3.0;
  const double beta_ = pi/4.0;
  const double gamma_ = pi/6.0;

  test_get_offset(UnitCell3D::create_triclinic_primitive(a_, b_, c_, alpha_, beta_, gamma_));
  test_get_offset(UnitCell3D::create_monoclinic_primitive(a_,  b_,  c_,  beta_));
  test_get_offset(UnitCell3D::create_monoclinic_base(a_,  b_,  c_,  beta_));
  test_get_offset(UnitCell3D::create_orthorhombic_primitive(a_,  b_,  c_));
  test_get_offset(UnitCell3D::create_orthorhombic_base(a_,  b_,  c_));
  test_get_offset(UnitCell3D::create_orthorhombic_body(a_,  b_,  c_));
  test_get_offset(UnitCell3D::create_orthorhombic_face(a_,  b_,  c_));
  test_get_offset(UnitCell3D::create_tetragonal_primitive(a_,  c_));
  test_get_offset(UnitCell3D::create_tetragonal_body(a_,  c_));
  test_get_offset(UnitCell3D::create_rhombohedral_centered(a_, alpha_));
  test_get_offset(UnitCell3D::create_hexagonal_primitive(a_, c_));
  test_get_offset(UnitCell3D::create_cubic_primitive(a_));
  test_get_offset(UnitCell3D::create_cubic_body(a_));
  test_get_offset(UnitCell3D::create_cubic_face(a_));
}
