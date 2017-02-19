#ifndef LATTICEGEN_CRYSTALLOGRAPHICPOINTGROUPS_H
#define LATTICEGEN_CRYSTALLOGRAPHICPOINTGROUPS_H

#include <vector>
#include <memory>

namespace Geometry
{
  enum CrystalSystem
  {
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Trigonal,
    Hexagonal,
    Cubic
  };
  enum CrystalClass
  {
    Triclinic_Pedion,
    Triclinic_Pinacoid,

    Monoclinic_Sphenoid,
    Monoclinic_Dome,
    Monoclinic_Prism,

    Orthorhombic_Disphenoid,
    Orthorhombic_Pyramid,
    Orthorhombic_Dipyramid,

    Tetragonal_Pyramid,
    Tetragonal_Disphenoid,
    Tetragonal_Dipyramid,
    Tetragonal_Trapezohedron,
    Ditetragonal_Pyramid,
    Tetragonal_Scalenohedron,
    Ditetragonal_Dipyramid,

    Trigonal_Pyramid,
    Rombohedron,
    Trigonal_Trapezohedron,
    Ditrigonal_Pyramid,
    Ditrigonal_Scalenohedron,

    Hexagonal_Pyramid,
    Trigonal_Dipyramid,
    Hexagonal_Dipyramid,
    Hexagonal_Trapezohedron,
    Dihexagonal_Pyramid,
    Ditrigonal_Dipyramid,
    Dihexagonal_Dipyramid,

    Tetatroid,
    Diploid,
    Gyroid,
    Hexatetrahedron,
    Hexaoctahedron
  };

  CrystalSystem get_crystal_system(CrystalClass crystal_class);

  class CrystallographicPointGroup
  {
  public:
    enum SymmetryElement
    {
      Identity,
      Inversion,
      Rotation2100,
      Rotation2010,
      Rotation2001,
      Rotation2p300,
      Rotation2m300,
      Rotation2p600,
      Rotation2m600,
      Rotationp3001,
      Rotationm3001,
      Rotationp4001,
      Rotationm4001,
      Rotationp4010,
      Rotationm4010,
      Rotationp4100,
      Rotationm4100,
      Rotationp6001,
      Rotationm6001,
      Rotation2p1p10,
      Rotation2p1m10,
      Rotation2p10p1,
      Rotation2p10m1,
      Rotation20p1p1,
      Rotation20p1m1,
      Rotationp3p1p1p1,
      Rotationm3p1p1p1,
      Rotationp3p1m1m1,
      Rotationm3p1m1m1,
      Rotationp3m1p1m1,
      Rotationm3m1p1m1,
      Rotationp3m1m1p1,
      Rotationm3m1m1p1,
      Reflection100,
      Reflection010,
      Reflection001,
      Reflectionp1p10,
      Reflectionp1m10,
      Reflectionp10p1,
      Reflectionp10m1,
      Reflection0p1p1,
      Reflection0p1m1,
      Reflectionp300,
      Reflectionm300,
      Reflectionp600,
      Reflectionm600,
      ImproperRotationp3001,
      ImproperRotationm3001,
      ImproperRotationp4001,
      ImproperRotationm4001,
      ImproperRotationp4010,
      ImproperRotationm4010,
      ImproperRotationp4100,
      ImproperRotationm4100,
      ImproperRotationp6001,
      ImproperRotationm6001,
      ImproperRotationp6p1p1p1,
      ImproperRotationm6p1p1p1,
      ImproperRotationp6p1m1m1,
      ImproperRotationm6p1m1m1,
      ImproperRotationp6m1p1m1,
      ImproperRotationm6m1p1m1,
      ImproperRotationp6m1m1p1,
      ImproperRotationm6m1m1p1,
    };

    typedef std::vector<SymmetryElement> Elements;

    virtual CrystalClass get_crystal_class() const = 0;
    virtual Elements get_generators() const = 0;
    virtual Elements get_elements() const = 0;

    static std::unique_ptr<CrystallographicPointGroup> create(CrystalClass crystal_class);

    bool operator==(const CrystallographicPointGroup& other) const
    {
      return this->get_crystal_class() == other.get_crystal_class();
    }
  };

  template<class>
  class CrystallographicPointGroupX : public CrystallographicPointGroup
  {
  public:
    CrystalClass get_crystal_class() const override
    {
      return CrystallographicPointGroupX::crystal_class;
    }
    Elements get_generators() const override
    {
      return CrystallographicPointGroupX::generators;
    }
    Elements get_elements() const override
    {
      return CrystallographicPointGroupX::elements;
    }

  protected:
    static const CrystalClass crystal_class;
    static const Elements generators;
    static const Elements elements;
  };

  class C1 : public CrystallographicPointGroupX<C1> {};
  class Ci : public CrystallographicPointGroupX<Ci> {};

  class C2  : public CrystallographicPointGroupX<C2> {};
  class Cs  : public CrystallographicPointGroupX<Cs> {};
  class C2h : public CrystallographicPointGroupX<C2h> {};

  class D2  : public CrystallographicPointGroupX<D2> {};
  class C2v : public CrystallographicPointGroupX<C2v> {};
  class D2h : public CrystallographicPointGroupX<D2h> {};

  class C4  : public CrystallographicPointGroupX<C4> {};
  class S4  : public CrystallographicPointGroupX<S4> {};
  class C4h : public CrystallographicPointGroupX<C4h> {};
  class D4  : public CrystallographicPointGroupX<D4> {};
  class C4v : public CrystallographicPointGroupX<C4v> {};
  class D2d : public CrystallographicPointGroupX<D2d> {};
  class D4h : public CrystallographicPointGroupX<D4h> {};

  class C3  : public CrystallographicPointGroupX<C3> {};
  class S6  : public CrystallographicPointGroupX<S6> {};
  class D3  : public CrystallographicPointGroupX<D3> {};
  class C3v : public CrystallographicPointGroupX<C3v> {};
  class D3d : public CrystallographicPointGroupX<D3d> {};

  class C6 : public CrystallographicPointGroupX<C6> {};
  class C3h : public CrystallographicPointGroupX<C3h> {};
  class C6h : public CrystallographicPointGroupX<C6h> {};
  class D6 : public CrystallographicPointGroupX<D6> {};
  class C6v : public CrystallographicPointGroupX<C6v> {};
  class D3h : public CrystallographicPointGroupX<D3h> {};
  class D6h : public CrystallographicPointGroupX<D6h> {};

  class T : public CrystallographicPointGroupX<T> {};
  class Th : public CrystallographicPointGroupX<Th> {};
  class O : public CrystallographicPointGroupX<O> {};
  class Td : public CrystallographicPointGroupX<Td> {};
  class Oh : public CrystallographicPointGroupX<Oh> {};
}

#endif //LATTICEGEN_CRYSTALLOGRAPHICPOINTGROUPS_H
