#include "CrystallographicPointGroups.h"

namespace Geometry
{
  typedef CrystallographicPointGroup::Elements Elements;

#define POPULATE_GROUP_CLASS(group_, crystal_class_) \
    template<> const CrystalClass CrystallographicPointGroupX<group_>::crystal_class = crystal_class_;
#define POPULATE_GROUP_GENERATORS(group_, ...) \
    template<> const Elements CrystallographicPointGroupX<group_>::generators = Elements{__VA_ARGS__};
#define POPULATE_GROUP_ELEMENTS(group_, ...) \
    template<> const Elements CrystallographicPointGroupX<group_>::elements = Elements{__VA_ARGS__};

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

  POPULATE_GROUP_CLASS(D2, Orthorhombic_Didphenoid);
  POPULATE_GROUP_GENERATORS(D2, Rotation2100, Rotation2010);
  POPULATE_GROUP_ELEMENTS(D2, Identity, Rotation2100, Rotation2010, Rotation2001);

  POPULATE_GROUP_CLASS(C2v, Orthorhombic_Pyramid);
  POPULATE_GROUP_GENERATORS(C2v, Reflection100, Reflection010);
  POPULATE_GROUP_ELEMENTS(C2v, Identity, Reflection100, Reflection010, Rotation2001);

  POPULATE_GROUP_CLASS(D2h, Orthorhombic_Dipyramid);
  POPULATE_GROUP_GENERATORS(D2h, Reflection100, Reflection010, Reflection100);
  POPULATE_GROUP_ELEMENTS(D2h,
                          Identity, Reflection100, Reflection010, Reflection100,
                          Rotation2001, Rotation2010, Rotation2100, Inversion);

  ////////////////  Tetragonal  ////////////////

  POPULATE_GROUP_CLASS(C4, Tetragonal_Pyramid);
  POPULATE_GROUP_GENERATORS(C4, Rotationp4001);
  POPULATE_GROUP_ELEMENTS(C4, Identity, Rotationp4001, Rotation2001, Rotationm4001);

  POPULATE_GROUP_CLASS(S4, Tetragonal_Didphenoid);
  POPULATE_GROUP_GENERATORS(S4, ImproperRotationp4001);
  POPULATE_GROUP_ELEMENTS(S4, Identity, ImproperRotationp4001, Rotation2001, ImproperRotationm4001);

  POPULATE_GROUP_CLASS(C4h, Tetragonal_Dipyramid);
  POPULATE_GROUP_GENERATORS(C4h, Rotationp4001, Reflection001);
  POPULATE_GROUP_ELEMENTS(C4h,
                          Identity, Rotationp4001, Rotation2001, Rotationm4001,
                          Reflection001, ImproperRotationp4001, Inversion, ImproperRotationm4001);

  POPULATE_GROUP_CLASS(D4, Tetragonal_Trapezohedron);
  POPULATE_GROUP_GENERATORS(D4, Rotation2100, Rotation2m450);
  POPULATE_GROUP_ELEMENTS(D4,
                          Identity, Rotationp4001, Rotation2001, Rotationm4001,
                          Rotation2100, Rotation2m450, Rotation2010, Rotation2p450);

  POPULATE_GROUP_CLASS(C4v, Ditetragonal_Pyramid);
  POPULATE_GROUP_GENERATORS(C4v, Reflection100, Reflectionp450);
  POPULATE_GROUP_ELEMENTS(C4v,
                          Identity, Rotationp4001, Rotation2001, Rotationm4001,
                          Reflectionp450, Reflectionp450, Reflection010, Reflectionm450 );

  POPULATE_GROUP_CLASS(D2d, Tetragonal_Scalenohedron);
  POPULATE_GROUP_GENERATORS(D2d, Rotation2100, Reflectionp450);
  POPULATE_GROUP_ELEMENTS(D2d,
                          Identity, ImproperRotationp4001, Rotation2001, ImproperRotationm4001,
                          Rotation2100, Reflectionp450, Rotation2010, Reflectionm450);

  POPULATE_GROUP_CLASS(D4h, Ditetragonal_Dipyramid);
  POPULATE_GROUP_GENERATORS(D4h, Reflection100, Reflectionp450, Reflection001);
  POPULATE_GROUP_ELEMENTS(D4h,
                          Identity, Rotationp4001, Rotation2001, Rotationm4001,
                          Rotation2100, Rotation2p450, Rotation2010, Rotation2m450,
                          Reflection100, Reflectionm450, Reflection010, Reflectionp450,
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
  POPULATE_GROUP_GENERATORS(D3, Rotation2p300, Rotation2m300);
  POPULATE_GROUP_ELEMENTS(D3,
                          Identity, Rotationp3001, Rotationm3001,
                          Rotation2100, Rotation2p600, Rotation2m600);

  POPULATE_GROUP_CLASS(C3v, Ditrigonal_Pyramid);
  POPULATE_GROUP_GENERATORS(C3v, Reflection100, Reflectionp600);
  POPULATE_GROUP_ELEMENTS(C3v,
                          Identity, Rotationp3001, Rotationm3001,
                          Reflection100, Reflectionp600, Reflectionm600);

  POPULATE_GROUP_CLASS(D3d, Ditrigonal_Scalenohedron);
  POPULATE_GROUP_GENERATORS(D3d, Rotation2100, Reflectionp300);
  POPULATE_GROUP_ELEMENTS(D3d,
                          Identity, Rotationp3001, Rotationm3001,
                          Rotation2100, Rotation2p600, Rotation2m600,
                          Reflectionp300, Reflectionm300, Reflection010,
                          ImproperRotationp6001, ImproperRotationm6001, Inversion);

  ////////////////   Hexagonal  ////////////////

  POPULATE_GROUP_CLASS(C6, Hexagonal_Pyramid);
  POPULATE_GROUP_GENERATORS(C6, Rotationp6001);
  POPULATE_GROUP_ELEMENTS(C6,
                          Identity, Rotationp6001, Rotationm6001,
                          Rotationp3001, Rotationm3001, Rotation2001);

  POPULATE_GROUP_CLASS(C3h, Trigonal_Dipyramid);
  POPULATE_GROUP_GENERATORS(C3h, Rotationp3001, Reflection001);
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
  POPULATE_GROUP_GENERATORS(C6v, Reflection100, Reflectionp300);
  POPULATE_GROUP_ELEMENTS(C6v,
                          Identity, Rotationp6001, Rotationm6001,
                          Rotationp3001, Rotationm3001, Rotation2001,
                          Reflection100, Reflection010,
                          Reflectionp300, Reflectionm300,
                          Reflectionp600, Reflectionm600);

  POPULATE_GROUP_CLASS(D3h, Ditrigonal_Dipyramid);
  POPULATE_GROUP_GENERATORS(D3h, Reflection001, Reflection100, Reflectionp600);
  POPULATE_GROUP_ELEMENTS(D3h,
                          Identity, Rotationp3001, Rotationm3001,
                          Reflection010, Reflectionp600, Reflectionm600,
                          Reflection001, ImproperRotationp3001, ImproperRotationm3001,
                          Inversion, ImproperRotationp6001, ImproperRotationm6001);

  POPULATE_GROUP_CLASS(D6h, Dihexagonal_Dipyramid);
  POPULATE_GROUP_GENERATORS(D6h, Reflection001, Reflection100, Reflectionp300);
  POPULATE_GROUP_ELEMENTS(D6h,
                          Identity, Rotationp6001, Rotationm6001,
                          Rotationp3001, Rotationm3001, Rotation2001,
                          Reflection001, ImproperRotationp3001, ImproperRotationm3001,
                          Inversion, ImproperRotationp6001, ImproperRotationm6001,
                          Rotation2100, Rotation2p300, Rotation2m300,
                          Rotation2p600, Rotation2m600, Rotation2010,
                          Reflection100, Reflectionp600, Reflectionm600,
                          Reflectionp300, Reflectionm300, Reflection010);

  ////////////////     Cubic    ////////////////

  //TODO implement Cubic point groups

  POPULATE_GROUP_CLASS(T, Tetatroid);
  POPULATE_GROUP_GENERATORS(T, Identity);
  POPULATE_GROUP_ELEMENTS(T, Identity);

  POPULATE_GROUP_CLASS(Td, Diploid);
  POPULATE_GROUP_GENERATORS(Td, Identity);
  POPULATE_GROUP_ELEMENTS(Td, Identity);

  POPULATE_GROUP_CLASS(Th, Gyroid);
  POPULATE_GROUP_GENERATORS(Th, Identity);
  POPULATE_GROUP_ELEMENTS(Th, Identity);

  POPULATE_GROUP_CLASS(O, Hexatetrahedron);
  POPULATE_GROUP_GENERATORS(O, Identity);
  POPULATE_GROUP_ELEMENTS(O, Identity);

  POPULATE_GROUP_CLASS(Oh, Hexaoctahedron);
  POPULATE_GROUP_GENERATORS(Oh, Identity);
  POPULATE_GROUP_ELEMENTS(Oh, Identity);

}