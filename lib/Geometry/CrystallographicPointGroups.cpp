#include "CrystallographicPointGroups.h"

#include <Core/Exceptions.h>

using namespace Geometry;

CrystalSystem Geometry::get_crystal_system(CrystalClass crystal_class)
{
  switch (crystal_class) {
    case Triclinic_Pedion:
    case Triclinic_Pinacoid:
      return Triclinic;

    case Monoclinic_Sphenoid:
    case Monoclinic_Dome:
    case Monoclinic_Prism:
      return Monoclinic;

    case Orthorhombic_Disphenoid:
    case Orthorhombic_Pyramid:
    case Orthorhombic_Dipyramid:
      return Orthorhombic;

    case Tetragonal_Pyramid:
    case Tetragonal_Disphenoid:
    case Tetragonal_Dipyramid:
    case Tetragonal_Trapezohedron:
    case Ditetragonal_Pyramid:
    case Tetragonal_Scalenohedron:
    case Ditetragonal_Dipyramid:
      return Tetragonal;

    case Trigonal_Pyramid:
    case Rombohedron:
    case Trigonal_Trapezohedron:
    case Ditrigonal_Pyramid:
    case Ditrigonal_Scalenohedron:
      return Trigonal;

    case Hexagonal_Pyramid:
    case Trigonal_Dipyramid:
    case Hexagonal_Dipyramid:
    case Hexagonal_Trapezohedron:
    case Dihexagonal_Pyramid:
    case Ditrigonal_Dipyramid:
    case Dihexagonal_Dipyramid:
      return Hexagonal;

    case Tetatroid:
    case Diploid:
    case Gyroid:
    case Hexatetrahedron:
    case Hexaoctahedron:
      return Cubic;

    default:
      THROW_INVALID_ARGUMENT("Unrecognised crystal class on the input of get_crystal_system");
  }
}

std::unique_ptr<CrystallographicPointGroup> CrystallographicPointGroup::create(CrystalClass crystal_class)
{
  switch (crystal_class) {
    case Triclinic_Pedion:
      return std::make_unique<C1>();
    case Triclinic_Pinacoid:
      return std::make_unique<Ci>();
    case Monoclinic_Sphenoid:
      return std::make_unique<C2>();
    case Monoclinic_Dome:
      return std::make_unique<Cs>();
    case Monoclinic_Prism:
      return std::make_unique<C2h>();
    case Orthorhombic_Disphenoid:
      return std::make_unique<D2>();
    case Orthorhombic_Pyramid:
      return std::make_unique<C2v>();
    case Orthorhombic_Dipyramid:
      return std::make_unique<D2h>();
    case Tetragonal_Pyramid:
      return std::make_unique<C4>();
    case Tetragonal_Disphenoid:
      return std::make_unique<S4>();
    case Tetragonal_Dipyramid:
      return std::make_unique<C4h>();
    case Tetragonal_Trapezohedron:
      return std::make_unique<D4>();
    case Ditetragonal_Pyramid:
      return std::make_unique<C4v>();
    case Tetragonal_Scalenohedron:
      return std::make_unique<D2d>();
    case Ditetragonal_Dipyramid:
      return std::make_unique<D4h>();
    case Trigonal_Pyramid:
      return std::make_unique<C3>();
    case Rombohedron:
      return std::make_unique<S6>();
    case Trigonal_Trapezohedron:
      return std::make_unique<D3>();
    case Ditrigonal_Pyramid:
      return std::make_unique<C3v>();
    case Ditrigonal_Scalenohedron:
      return std::make_unique<D3d>();
    case Hexagonal_Pyramid:
      return std::make_unique<C6>();
    case Trigonal_Dipyramid:
      return std::make_unique<C3h>();
    case Hexagonal_Dipyramid:
      return std::make_unique<C6h>();
    case Hexagonal_Trapezohedron:
      return std::make_unique<D6>();
    case Dihexagonal_Pyramid:
      return std::make_unique<C6v>();
    case Ditrigonal_Dipyramid:
      return std::make_unique<D3h>();
    case Dihexagonal_Dipyramid:
      return std::make_unique<D6h>();
    case Tetatroid:
      return std::make_unique<T>();
    case Diploid:
      return std::make_unique<Th>();
    case Gyroid:
      return std::make_unique<O>();
    case Hexatetrahedron:
      return std::make_unique<Td>();
    case Hexaoctahedron:
      return std::make_unique<Oh>();

    default:
      THROW_INVALID_ARGUMENT("Unrecognised crystal class in CrystallographicPointGroup::create");
  }
}

typedef CrystallographicPointGroup::Elements Elements;

#define POPULATE_GROUP_CLASS(group_, crystal_class_) \
  template<> const CrystalClass CrystallographicPointGroupX<group_>::crystal_class = crystal_class_
#define POPULATE_GROUP_GENERATORS(group_, ...) \
  template<> const Elements CrystallographicPointGroupX<group_>::generators = Elements{__VA_ARGS__}
#define POPULATE_GROUP_ELEMENTS(group_, ...) \
  template<> const Elements CrystallographicPointGroupX<group_>::elements = Elements{__VA_ARGS__}

////////////////   Triclinic  ////////////////

POPULATE_GROUP_CLASS(C1, Triclinic_Pedion);
POPULATE_GROUP_GENERATORS(C1, Identity);
POPULATE_GROUP_ELEMENTS(C1, Identity);

POPULATE_GROUP_CLASS(Ci, Triclinic_Pinacoid);
POPULATE_GROUP_GENERATORS(Ci, Inversion);
POPULATE_GROUP_ELEMENTS(Ci, Identity, Inversion);

////////////////  Monoclinic  ////////////////

POPULATE_GROUP_CLASS(C2, Monoclinic_Sphenoid);
POPULATE_GROUP_GENERATORS(C2, Rotation2001);
POPULATE_GROUP_ELEMENTS(C2, Identity, Rotation2001);

POPULATE_GROUP_CLASS(Cs, Monoclinic_Dome);
POPULATE_GROUP_GENERATORS(Cs, Reflection001);
POPULATE_GROUP_ELEMENTS(Cs, Identity, Reflection001);

POPULATE_GROUP_CLASS(C2h, Monoclinic_Prism);
POPULATE_GROUP_GENERATORS(C2h, Rotation2001, Reflection001);
POPULATE_GROUP_ELEMENTS(C2h, Identity, Rotation2001, Reflection001, Inversion);

//////////////// Orthorhombic ////////////////

POPULATE_GROUP_CLASS(D2, Orthorhombic_Disphenoid);
POPULATE_GROUP_GENERATORS(D2, Rotation2100, Rotation2010);
POPULATE_GROUP_ELEMENTS(D2, Identity, Rotation2100, Rotation2010, Rotation2001);

POPULATE_GROUP_CLASS(C2v, Orthorhombic_Pyramid);
POPULATE_GROUP_GENERATORS(C2v, Reflection100, Reflection010);
POPULATE_GROUP_ELEMENTS(C2v, Identity, Reflection100, Reflection010, Rotation2001);

POPULATE_GROUP_CLASS(D2h, Orthorhombic_Dipyramid);
POPULATE_GROUP_GENERATORS(D2h, Reflection001, Reflection010, Reflection100);
POPULATE_GROUP_ELEMENTS(D2h,
                        Identity, Reflection001, Reflection010, Reflection100,
                        Rotation2001, Rotation2010, Rotation2100, Inversion);

////////////////  Tetragonal  ////////////////

POPULATE_GROUP_CLASS(C4, Tetragonal_Pyramid);
POPULATE_GROUP_GENERATORS(C4, Rotationp4001);
POPULATE_GROUP_ELEMENTS(C4, Identity, Rotationp4001, Rotation2001, Rotationm4001);

POPULATE_GROUP_CLASS(S4, Tetragonal_Disphenoid);
POPULATE_GROUP_GENERATORS(S4, ImproperRotationp4001);
POPULATE_GROUP_ELEMENTS(S4, Identity, ImproperRotationp4001, Rotation2001, ImproperRotationm4001);

POPULATE_GROUP_CLASS(C4h, Tetragonal_Dipyramid);
POPULATE_GROUP_GENERATORS(C4h, Rotationp4001, Reflection001);
POPULATE_GROUP_ELEMENTS(C4h,
                        Identity, Rotationp4001, Rotation2001, Rotationm4001,
                        Reflection001, ImproperRotationp4001, Inversion, ImproperRotationm4001);

POPULATE_GROUP_CLASS(D4, Tetragonal_Trapezohedron);
POPULATE_GROUP_GENERATORS(D4, Rotation2100, Rotation2p1m10);
POPULATE_GROUP_ELEMENTS(D4,
                        Identity, Rotationp4001, Rotation2001, Rotationm4001,
                        Rotation2100, Rotation2p1m10, Rotation2010, Rotation2p1p10);

POPULATE_GROUP_CLASS(C4v, Ditetragonal_Pyramid);
POPULATE_GROUP_GENERATORS(C4v, Reflection100, Reflectionp1m10);
POPULATE_GROUP_ELEMENTS(C4v,
                        Identity, Rotationp4001, Rotation2001, Rotationm4001,
                        Reflection100, Reflectionp1m10, Reflection010, Reflectionp1p10);

POPULATE_GROUP_CLASS(D2d, Tetragonal_Scalenohedron);
POPULATE_GROUP_GENERATORS(D2d, Rotation2100, Reflectionp1m10);
POPULATE_GROUP_ELEMENTS(D2d,
                        Identity, ImproperRotationp4001, Rotation2001, ImproperRotationm4001,
                        Rotation2100, Reflectionp1m10, Rotation2010, Reflectionp1p10);

POPULATE_GROUP_CLASS(D4h, Ditetragonal_Dipyramid);
POPULATE_GROUP_GENERATORS(D4h, Reflection100, Reflectionp1m10, Reflection001);
POPULATE_GROUP_ELEMENTS(D4h,
                        Identity, Rotationp4001, Rotation2001, Rotationm4001,
                        Rotation2100, Rotation2p1p10, Rotation2010, Rotation2p1m10,
                        Reflection100, Reflectionp1p10, Reflection010, Reflectionp1m10,
                        ImproperRotationp4001, ImproperRotationm4001, Inversion, Reflection001);

////////////////   Trigonal   ////////////////

POPULATE_GROUP_CLASS(C3, Trigonal_Pyramid);
POPULATE_GROUP_GENERATORS(C3, Rotationp3001);
POPULATE_GROUP_ELEMENTS(C3, Identity, Rotationp3001, Rotationm3001);

POPULATE_GROUP_CLASS(S6, Rombohedron);
POPULATE_GROUP_GENERATORS(S6, ImproperRotationp6001);
POPULATE_GROUP_ELEMENTS(S6,
                        Identity, ImproperRotationp6001, ImproperRotationm6001,
                        Rotationp3001, Rotationm3001, Inversion);

POPULATE_GROUP_CLASS(D3, Trigonal_Trapezohedron);
POPULATE_GROUP_GENERATORS(D3, Rotation2100, Rotation2p600);
POPULATE_GROUP_ELEMENTS(D3,
                        Identity, Rotationp3001, Rotationm3001,
                        Rotation2100, Rotation2p600, Rotation2m600);

POPULATE_GROUP_CLASS(C3v, Ditrigonal_Pyramid);
POPULATE_GROUP_GENERATORS(C3v, Reflection010, Reflectionm300);
POPULATE_GROUP_ELEMENTS(C3v,
                        Identity, Rotationp3001, Rotationm3001,
                        Reflection010, Reflectionm300, Reflectionp300);

POPULATE_GROUP_CLASS(D3d, Ditrigonal_Scalenohedron);
POPULATE_GROUP_GENERATORS(D3d, Rotation2100, Reflectionm600);
POPULATE_GROUP_ELEMENTS(D3d,
                        Identity, Rotationp3001, Rotationm3001,
                        Rotation2100, Rotation2p600, Rotation2m600,
                        Reflectionm600, Reflectionp600, Reflection100,
                        ImproperRotationp6001, ImproperRotationm6001, Inversion);

////////////////   Hexagonal  ////////////////

POPULATE_GROUP_CLASS(C6, Hexagonal_Pyramid);
POPULATE_GROUP_GENERATORS(C6, Rotationp6001);
POPULATE_GROUP_ELEMENTS(C6,
                        Identity, Rotationp6001, Rotationm6001,
                        Rotationp3001, Rotationm3001, Rotation2001);

POPULATE_GROUP_CLASS(C3h, Trigonal_Dipyramid);
POPULATE_GROUP_GENERATORS(C3h, ImproperRotationp3001);
POPULATE_GROUP_ELEMENTS(C3h,
                        Identity, Rotationp3001, Rotationm3001,
                        Reflection001, ImproperRotationp3001, ImproperRotationm3001);

POPULATE_GROUP_CLASS(C6h, Hexagonal_Dipyramid);
POPULATE_GROUP_GENERATORS(C6h, Rotationp6001, Reflection001);
POPULATE_GROUP_ELEMENTS(C6h,
                        Identity, Rotationp6001, Rotationm6001,
                        Rotationp3001, Rotationm3001, Rotation2001,
                        Reflection001, ImproperRotationp3001, ImproperRotationm3001,
                        Inversion, ImproperRotationp6001, ImproperRotationm6001);

POPULATE_GROUP_CLASS(D6, Hexagonal_Trapezohedron);
POPULATE_GROUP_GENERATORS(D6, Rotation2100, Rotation2p300);
POPULATE_GROUP_ELEMENTS(D6,
                        Identity, Rotationp3001, Rotationm3001,
                        Rotationp6001, Rotationm6001, Rotation2001,
                        Rotation2100, Rotation2p300, Rotation2m300,
                        Rotation2p600, Rotation2m600, Rotation2010);

POPULATE_GROUP_CLASS(C6v, Dihexagonal_Pyramid);
POPULATE_GROUP_GENERATORS(C6v, Reflection010, Reflectionm600);
POPULATE_GROUP_ELEMENTS(C6v,
                        Identity, Rotationp6001, Rotationm6001,
                        Rotationp3001, Rotationm3001, Rotation2001,
                        Reflection010, Reflection100,
                        Reflectionp300, Reflectionm300,
                        Reflectionp600, Reflectionm600);

POPULATE_GROUP_CLASS(D3h, Ditrigonal_Dipyramid);
POPULATE_GROUP_GENERATORS(D3h, Reflection001, Reflection010, Reflectionm300);
POPULATE_GROUP_ELEMENTS(D3h,
                        Identity, Rotationp3001, Rotationm3001,
                        Reflection010, Reflectionm300, Reflectionp300,
                        Reflection001, ImproperRotationp3001, ImproperRotationm3001,
                        Rotation2100, Rotation2p600, Rotation2m600);

POPULATE_GROUP_CLASS(D6h, Dihexagonal_Dipyramid);
POPULATE_GROUP_GENERATORS(D6h, Reflection001, Reflection010, Reflectionm600);
POPULATE_GROUP_ELEMENTS(D6h,
                        Identity, Rotationp6001, Rotationm6001,
                        Rotationp3001, Rotationm3001, Rotation2001,
                        Reflection001, ImproperRotationp3001, ImproperRotationm3001,
                        Inversion, ImproperRotationp6001, ImproperRotationm6001,
                        Rotation2100, Rotation2p300, Rotation2m300,
                        Rotation2p600, Rotation2m600, Rotation2010,
                        Reflection010, Reflectionp600, Reflectionm600,
                        Reflectionp300, Reflectionm300, Reflection100);

////////////////     Cubic    ////////////////

POPULATE_GROUP_CLASS(T, Tetatroid);
POPULATE_GROUP_GENERATORS(T, Rotation2001, Rotationp3p1p1p1);
POPULATE_GROUP_ELEMENTS(T,
                        Identity, Rotation2100, Rotation2010, Rotation2001,
                        Rotationp3p1p1p1, Rotationp3p1m1m1, Rotationp3m1p1m1, Rotationp3m1m1p1,
                        Rotationm3p1p1p1, Rotationm3p1m1m1, Rotationm3m1p1m1, Rotationm3m1m1p1);

POPULATE_GROUP_CLASS(Th, Diploid);
POPULATE_GROUP_GENERATORS(Th, Rotation2100, Rotationp3p1p1p1, Inversion);
POPULATE_GROUP_ELEMENTS(Th,
                        Identity, Rotation2100, Rotation2010, Rotation2001,
                        Rotationp3p1p1p1, Rotationp3p1m1m1, Rotationp3m1p1m1, Rotationp3m1m1p1,
                        Rotationm3p1p1p1, Rotationm3p1m1m1, Rotationm3m1p1m1, Rotationm3m1m1p1,
                        Inversion, Reflection001, Reflection010, Reflection100,
                        ImproperRotationp6p1p1p1, ImproperRotationp6p1m1m1, ImproperRotationp6m1p1m1,
                        ImproperRotationp6m1m1p1,
                        ImproperRotationm6p1p1p1, ImproperRotationm6p1m1m1, ImproperRotationm6m1p1m1,
                        ImproperRotationm6m1m1p1);

POPULATE_GROUP_CLASS(O, Gyroid);
POPULATE_GROUP_GENERATORS(O, Rotation2100, Rotationp3p1p1p1, Rotation2p1p10);
POPULATE_GROUP_ELEMENTS(O,
                        Identity, Rotation2100, Rotation2010, Rotation2001,
                        Rotationp3p1p1p1, Rotationp3p1m1m1, Rotationp3m1p1m1, Rotationp3m1m1p1,
                        Rotationm3p1p1p1, Rotationm3p1m1m1, Rotationm3m1p1m1, Rotationm3m1m1p1,
                        Rotationp4001, Rotationm4001, Rotationp4010, Rotationm4010, Rotationp4100, Rotationm4100,
                        Rotation2p1p10, Rotation2p1m10, Rotation2p10p1,
                        Rotation2p10m1, Rotation20p1p1, Rotation20p1m1);

POPULATE_GROUP_CLASS(Td, Hexatetrahedron);
POPULATE_GROUP_GENERATORS(Td, Rotation2100, Rotationp3p1p1p1, Reflectionp1m10);
POPULATE_GROUP_ELEMENTS(Td,
                        Identity, Rotation2100, Rotation2010, Rotation2001,
                        Rotationp3p1p1p1, Rotationp3p1m1m1, Rotationp3m1p1m1, Rotationp3m1m1p1,
                        Rotationm3p1p1p1, Rotationm3p1m1m1, Rotationm3m1p1m1, Rotationm3m1m1p1,
                        ImproperRotationp4001, ImproperRotationp4010, ImproperRotationp4100,
                        ImproperRotationm4001, ImproperRotationm4010, ImproperRotationm4100,
                        Reflectionp1p10, Reflectionp10p1, Reflection0p1p1,
                        Reflectionp1m10, Reflectionp10m1, Reflection0p1m1);

POPULATE_GROUP_CLASS(Oh, Hexaoctahedron);
POPULATE_GROUP_GENERATORS(Oh, Rotation2100, Rotationp3p1p1p1, Rotation2p1p10, Inversion);
POPULATE_GROUP_ELEMENTS(Oh,
                        Identity, Rotation2100, Rotation2010, Rotation2001,
                        Rotationp3p1p1p1, Rotationp3p1m1m1, Rotationp3m1p1m1, Rotationp3m1m1p1,
                        Rotationm3p1p1p1, Rotationm3p1m1m1, Rotationm3m1p1m1, Rotationm3m1m1p1,
                        Rotation2p1p10, Rotation2p1m10, Rotation2p10p1,
                        Rotation2p10m1, Rotation20p1p1, Rotation20p1m1,
                        Rotationp4001, Rotationp4010, Rotationp4100,
                        Rotationm4001, Rotationm4010, Rotationm4100,
                        Inversion, Reflection001, Reflection010, Reflection100,
                        ImproperRotationp6p1p1p1, ImproperRotationp6p1m1m1, ImproperRotationp6m1p1m1,
                        ImproperRotationp6m1m1p1,
                        ImproperRotationm6p1p1p1, ImproperRotationm6p1m1m1, ImproperRotationm6m1p1m1,
                        ImproperRotationm6m1m1p1,
                        Reflectionp1p10, Reflectionp10p1, Reflection0p1p1,
                        Reflectionp1m10, Reflectionp10m1, Reflection0p1m1,
                        ImproperRotationp4001, ImproperRotationp4010, ImproperRotationp4100,
                        ImproperRotationm4001, ImproperRotationm4010, ImproperRotationm4100);
