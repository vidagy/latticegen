#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Geometry/CrystallographicPointGroups.h>

namespace
{
  static const double pi = 3.14159265358979323846;
}

using namespace Geometry;
typedef CrystallographicPointGroup::Elements Elements;

TEST(TestCrystallographicPointGroups,TestViaInterface)
{
  const static std::unique_ptr<CrystallographicPointGroup> c1 = std::make_unique<C1>();
  EXPECT_EQ(c1->get_crystal_class(), Triclinic_Pedion);
}

namespace
{
  template<class Group, class Range = Elements>
  void testCrystallographicPointGroup(
    const CrystalClass& crystal_class,
    const CrystalSystem& crystal_system,
    const Range& generators,
    const Range& elements)
  {
    const static Group group = Group();
    EXPECT_EQ(crystal_class, group.get_crystal_class());
    EXPECT_EQ(crystal_system, get_crystal_system(group.get_crystal_class()));
    EXPECT_THAT(generators, ::testing::UnorderedElementsAreArray(group.get_generators()) );
    EXPECT_THAT(elements, ::testing::UnorderedElementsAreArray(group.get_elements()) );
  };
}

////////////////   Triclinic  ////////////////

TEST(TestCrystallographicPointGroups,C1)
{
  testCrystallographicPointGroup<C1>(
    Triclinic_Pedion,
    Triclinic,
    { CrystallographicPointGroup::Identity },
    { CrystallographicPointGroup::Identity }
  );
}

TEST(TestCrystallographicPointGroups,Ci)
{
  testCrystallographicPointGroup<Ci>(
    Triclinic_Pinacoid,
    Triclinic,
    { CrystallographicPointGroup::Inversion },
    { CrystallographicPointGroup::Inversion,
      CrystallographicPointGroup::Identity }
  );
}

////////////////  Monoclinic  ////////////////

TEST(TestCrystallographicPointGroups,C2)
{
  testCrystallographicPointGroup<C2>(
    Monoclinic_Sphenoid,
    Monoclinic,
    { CrystallographicPointGroup::Rotation2001 },
    { CrystallographicPointGroup::Rotation2001,
      CrystallographicPointGroup::Identity }
  );
}

TEST(TestCrystallographicPointGroups,Cs)
{
  testCrystallographicPointGroup<Cs>(
    Monoclinic_Dome,
    Monoclinic,
    { CrystallographicPointGroup::Reflection001 },
    { CrystallographicPointGroup::Reflection001,
      CrystallographicPointGroup::Identity }
  );
}

TEST(TestCrystallographicPointGroups,C2h)
{
  testCrystallographicPointGroup<C2h>(
    Monoclinic_Prism,
    Monoclinic,
    { CrystallographicPointGroup::Rotation2001,
      CrystallographicPointGroup::Reflection001 },
    { CrystallographicPointGroup::Identity,
      CrystallographicPointGroup::Rotation2001,
      CrystallographicPointGroup::Reflection001,
      CrystallographicPointGroup::Inversion }
  );
}

//////////////// Orthorhombic ////////////////

TEST(TestCrystallographicPointGroups,D2)
{
  testCrystallographicPointGroup<D2>(
    Orthorhombic_Didphenoid,
    Orthorhombic,
    { CrystallographicPointGroup::Rotation2100,
      CrystallographicPointGroup::Rotation2010 },
    { CrystallographicPointGroup::Identity,
      CrystallographicPointGroup::Rotation2100,
      CrystallographicPointGroup::Rotation2010,
      CrystallographicPointGroup::Rotation2001 }
  );
}

TEST(TestCrystallographicPointGroups,C2v)
{
  testCrystallographicPointGroup<C2v>(
    Orthorhombic_Pyramid,
    Orthorhombic,
    { CrystallographicPointGroup::Reflection100, CrystallographicPointGroup::Reflection010 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Reflection100,
      CrystallographicPointGroup::Reflection010, CrystallographicPointGroup::Rotation2001 }
  );
}

TEST(TestCrystallographicPointGroups,D2h)
{
  testCrystallographicPointGroup<D2h>(
    Orthorhombic_Dipyramid,
    Orthorhombic,
    { CrystallographicPointGroup::Reflection100,
      CrystallographicPointGroup::Reflection010, CrystallographicPointGroup::Reflection100 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Reflection100,
      CrystallographicPointGroup::Reflection010, CrystallographicPointGroup::Reflection100,
      CrystallographicPointGroup::Rotation2001, CrystallographicPointGroup::Rotation2010,
      CrystallographicPointGroup::Rotation2100, CrystallographicPointGroup::Inversion }
  );
}

////////////////  Tetragonal  ////////////////

TEST(TestCrystallographicPointGroups,C4)
{
  testCrystallographicPointGroup<C4>(
    Tetragonal_Pyramid,
    Tetragonal,
    { CrystallographicPointGroup::Rotationp4001 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp4001,
      CrystallographicPointGroup::Rotation2001, CrystallographicPointGroup::Rotationm4001 }
  );
}

TEST(TestCrystallographicPointGroups,S4)
{
  testCrystallographicPointGroup<S4>(
    Tetragonal_Didphenoid,
    Tetragonal,
    { CrystallographicPointGroup::ImproperRotationp4001 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::ImproperRotationp4001,
      CrystallographicPointGroup::Rotation2001, CrystallographicPointGroup::ImproperRotationm4001 }
  );
}

TEST(TestCrystallographicPointGroups,C4h)
{
  testCrystallographicPointGroup<C4h>(
    Tetragonal_Dipyramid,
    Tetragonal,
    { CrystallographicPointGroup::Rotationp4001, CrystallographicPointGroup::Reflection001 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp4001,
      CrystallographicPointGroup::Rotation2001, CrystallographicPointGroup::Rotationm4001,
      CrystallographicPointGroup::Reflection001, CrystallographicPointGroup::ImproperRotationp4001,
      CrystallographicPointGroup::Inversion, CrystallographicPointGroup::ImproperRotationm4001 }
  );
}

TEST(TestCrystallographicPointGroups,D4)
{
  testCrystallographicPointGroup<D4>(
    Tetragonal_Trapezohedron,
    Tetragonal,
    { CrystallographicPointGroup::Rotation2100, CrystallographicPointGroup::Rotation2m450 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp4001,
      CrystallographicPointGroup::Rotation2001, CrystallographicPointGroup::Rotationm4001,
      CrystallographicPointGroup::Rotation2100, CrystallographicPointGroup::Rotation2m450,
      CrystallographicPointGroup::Rotation2010, CrystallographicPointGroup::Rotation2p450 }
  );
}

TEST(TestCrystallographicPointGroups,C4v)
{
  testCrystallographicPointGroup<C4v>(
    Ditetragonal_Pyramid,
    Tetragonal,
    { CrystallographicPointGroup::Reflection100, CrystallographicPointGroup::Reflectionp450  },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp4001,
      CrystallographicPointGroup::Rotation2001, CrystallographicPointGroup::Rotationm4001,
      CrystallographicPointGroup::Reflectionp450, CrystallographicPointGroup::Reflectionp450,
      CrystallographicPointGroup::Reflection010, CrystallographicPointGroup::Reflectionm450  }
  );
}

TEST(TestCrystallographicPointGroups,D2d)
{
  testCrystallographicPointGroup<D2d>(
    Tetragonal_Scalenohedron,
    Tetragonal,
    { CrystallographicPointGroup::Rotation2100, CrystallographicPointGroup::Reflectionp450 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::ImproperRotationp4001,
      CrystallographicPointGroup::Rotation2001, CrystallographicPointGroup::ImproperRotationm4001,
      CrystallographicPointGroup::Rotation2100, CrystallographicPointGroup::Reflectionp450,
      CrystallographicPointGroup::Rotation2010, CrystallographicPointGroup::Reflectionm450}
  );
}

TEST(TestCrystallographicPointGroups,D4h)
{
  testCrystallographicPointGroup<D4h>(
    Ditetragonal_Dipyramid,
    Tetragonal,
    { CrystallographicPointGroup::Reflection100,
      CrystallographicPointGroup::Reflectionp450, CrystallographicPointGroup::Reflection001 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp4001,
      CrystallographicPointGroup::Rotation2001, CrystallographicPointGroup::Rotationm4001,
      CrystallographicPointGroup::Rotation2100, CrystallographicPointGroup::Rotation2p450,
      CrystallographicPointGroup::Rotation2010, CrystallographicPointGroup::Rotation2m450,
      CrystallographicPointGroup::Reflection100, CrystallographicPointGroup::Reflectionm450,
      CrystallographicPointGroup::Reflection010, CrystallographicPointGroup::Reflectionp450,
      CrystallographicPointGroup::ImproperRotationp4001, CrystallographicPointGroup::ImproperRotationm4001,
      CrystallographicPointGroup::Inversion, CrystallographicPointGroup::Reflection001 }
  );
}

////////////////   Trigonal   ////////////////

TEST(TestCrystallographicPointGroups,C3)
{
  testCrystallographicPointGroup<C3>(
    Trigonal_Pyramid,
    Trigonal,
    { CrystallographicPointGroup::Rotationp3001 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp3001,
      CrystallographicPointGroup::Rotationm3001 }
  );
}

TEST(TestCrystallographicPointGroups,S6)
{
  testCrystallographicPointGroup<S6>(
    Rombohedron,
    Trigonal,
    { CrystallographicPointGroup::ImproperRotationp6001 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::ImproperRotationp6001,
      CrystallographicPointGroup::ImproperRotationm6001, CrystallographicPointGroup::Rotationp3001,
      CrystallographicPointGroup::Rotationm3001, CrystallographicPointGroup::Inversion }
  );
}

TEST(TestCrystallographicPointGroups,D3)
{
  testCrystallographicPointGroup<D3>(
    Trigonal_Trapezohedron,
    Trigonal,
    { CrystallographicPointGroup::Rotation2p300, CrystallographicPointGroup::Rotation2m300  },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp3001,
      CrystallographicPointGroup::Rotationm3001, CrystallographicPointGroup::Rotation2100,
      CrystallographicPointGroup::Rotation2p600, CrystallographicPointGroup::Rotation2m600 }
  );
}

TEST(TestCrystallographicPointGroups,C3v)
{
  testCrystallographicPointGroup<C3v>(
    Ditrigonal_Pyramid,
    Trigonal,
    { CrystallographicPointGroup::Reflection100, CrystallographicPointGroup::Reflectionp600 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp3001,
      CrystallographicPointGroup::Rotationm3001, CrystallographicPointGroup::Reflection100,
      CrystallographicPointGroup::Reflectionp600, CrystallographicPointGroup::Reflectionm600 }
  );
}

TEST(TestCrystallographicPointGroups,D3d)
{
  testCrystallographicPointGroup<D3d>(
    Ditrigonal_Scalenohedron,
    Trigonal,
    { CrystallographicPointGroup::Rotation2100, CrystallographicPointGroup::Reflectionp300},
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp3001,
      CrystallographicPointGroup::Rotationm3001,CrystallographicPointGroup::Rotation2100,
      CrystallographicPointGroup::Rotation2p600, CrystallographicPointGroup::Rotation2m600,
      CrystallographicPointGroup::Reflectionp300, CrystallographicPointGroup::Reflectionm300,
      CrystallographicPointGroup::Reflection010, CrystallographicPointGroup::ImproperRotationp6001,
      CrystallographicPointGroup::ImproperRotationm6001, CrystallographicPointGroup::Inversion }
  );
}

////////////////  Hexagonal  ////////////////

TEST(TestCrystallographicPointGroups,C6)
{
  testCrystallographicPointGroup<C6>(
    Hexagonal_Pyramid,
    Hexagonal,
    { CrystallographicPointGroup::Rotationp6001 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp6001,
      CrystallographicPointGroup::Rotationm6001, CrystallographicPointGroup::Rotationp3001,
      CrystallographicPointGroup::Rotationm3001, CrystallographicPointGroup::Rotation2001 }
  );
}

TEST(TestCrystallographicPointGroups,C3h)
{
  testCrystallographicPointGroup<C3h>(
    Trigonal_Dipyramid,
    Hexagonal,
    { CrystallographicPointGroup::Rotationp3001, CrystallographicPointGroup::Reflection001 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp3001,
      CrystallographicPointGroup::Rotationm3001, CrystallographicPointGroup::Reflection001,
      CrystallographicPointGroup::ImproperRotationp3001, CrystallographicPointGroup::ImproperRotationm3001 }
  );
}

TEST(TestCrystallographicPointGroups,C6h)
{
  testCrystallographicPointGroup<C6h>(
    Hexagonal_Dipyramid,
    Hexagonal,
    { CrystallographicPointGroup::Rotationp6001, CrystallographicPointGroup::Reflection001 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp6001,
      CrystallographicPointGroup::Rotationm6001, CrystallographicPointGroup::Rotationp3001,
      CrystallographicPointGroup::Rotationm3001, CrystallographicPointGroup::Rotation2001,
      CrystallographicPointGroup::Reflection001, CrystallographicPointGroup::ImproperRotationp3001,
      CrystallographicPointGroup::ImproperRotationm3001, CrystallographicPointGroup::Inversion,
      CrystallographicPointGroup::ImproperRotationp6001, CrystallographicPointGroup::ImproperRotationm6001 }
  );
}

TEST(TestCrystallographicPointGroups,D6)
{
  testCrystallographicPointGroup<D6>(
    Hexagonal_Trapezohedron,
    Hexagonal,
    { CrystallographicPointGroup::Rotation2100, CrystallographicPointGroup::Rotation2p300 },
    { CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp3001,
      CrystallographicPointGroup::Rotationm3001, CrystallographicPointGroup::Rotationp6001,
      CrystallographicPointGroup::Rotationm6001, CrystallographicPointGroup::Rotation2001,
      CrystallographicPointGroup::Rotation2100, CrystallographicPointGroup::Rotation2p300,
      CrystallographicPointGroup::Rotation2m300, CrystallographicPointGroup::Rotation2p600,
      CrystallographicPointGroup::Rotation2m600, CrystallographicPointGroup::Rotation2010 }
  );
}

TEST(TestCrystallographicPointGroups,C6v)
{
  testCrystallographicPointGroup<C6v>(
    Dihexagonal_Pyramid,
    Hexagonal,
    { CrystallographicPointGroup::Reflection100, CrystallographicPointGroup::Reflectionp300  },
    {
      CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp6001,
      CrystallographicPointGroup::Rotationm6001, CrystallographicPointGroup::Rotationp3001,
      CrystallographicPointGroup::Rotationm3001,  CrystallographicPointGroup::Rotation2001,
      CrystallographicPointGroup::Reflection100, CrystallographicPointGroup::Reflection010,
      CrystallographicPointGroup::Reflectionp300, CrystallographicPointGroup::Reflectionm300,
      CrystallographicPointGroup::Reflectionp600, CrystallographicPointGroup::Reflectionm600 }
  );
}

TEST(TestCrystallographicPointGroups,D3h)
{
  testCrystallographicPointGroup<D3h>(
    Ditrigonal_Dipyramid,
    Hexagonal,
    { CrystallographicPointGroup::Reflection001, CrystallographicPointGroup::Reflection100,
      CrystallographicPointGroup::Reflectionp600 },
    {
      CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp3001,
      CrystallographicPointGroup::Rotationm3001, CrystallographicPointGroup::Reflection010,
      CrystallographicPointGroup::Reflectionp600, CrystallographicPointGroup::Reflectionm600,
      CrystallographicPointGroup::Reflection001, CrystallographicPointGroup::ImproperRotationp3001,
      CrystallographicPointGroup::ImproperRotationm3001, CrystallographicPointGroup::Inversion,
      CrystallographicPointGroup::ImproperRotationp6001, CrystallographicPointGroup::ImproperRotationm6001 }
  );
}

TEST(TestCrystallographicPointGroups,D6h)
{
  testCrystallographicPointGroup<D6h>(
    Dihexagonal_Dipyramid,
    Hexagonal,
    { CrystallographicPointGroup::Reflection001, CrystallographicPointGroup::Reflection100,
      CrystallographicPointGroup::Reflectionp300 },
    {
      CrystallographicPointGroup::Identity, CrystallographicPointGroup::Rotationp6001,
      CrystallographicPointGroup::Rotationm6001, CrystallographicPointGroup::Rotationp3001,
      CrystallographicPointGroup::Rotationm3001, CrystallographicPointGroup::Rotation2001,
      CrystallographicPointGroup::Reflection001, CrystallographicPointGroup::ImproperRotationp3001,
      CrystallographicPointGroup::ImproperRotationm3001, CrystallographicPointGroup::Inversion,
      CrystallographicPointGroup::ImproperRotationp6001, CrystallographicPointGroup::ImproperRotationm6001,
      CrystallographicPointGroup::Rotation2100, CrystallographicPointGroup::Rotation2p300,
      CrystallographicPointGroup::Rotation2m300, CrystallographicPointGroup::Rotation2p600,
      CrystallographicPointGroup::Rotation2m600, CrystallographicPointGroup::Rotation2010,
      CrystallographicPointGroup::Reflection100, CrystallographicPointGroup::Reflectionp600,
      CrystallographicPointGroup::Reflectionm600, CrystallographicPointGroup::Reflectionp300,
      CrystallographicPointGroup::Reflectionm300, CrystallographicPointGroup::Reflection010 }
  );
}

////////////////     Cubic    ////////////////

//TODO implement Cubic point groups
