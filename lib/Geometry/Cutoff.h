#ifndef LATTICEGEN_CUTOFF_H
#define LATTICEGEN_CUTOFF_H

#include <Geometry/UnitCell3D.h>

namespace Geometry
{
  ///@brief Class to define a domain in space
  class Cutoff
  {
  public:
    virtual bool is_included(const Point3D& point) const = 0;
    virtual std::tuple<size_t, size_t, size_t> steps_to_cover(const UnitCell3D& unit_cell_) const = 0;
  };

  class CutoffCube : public Cutoff
  {
  public:
    CutoffCube(const double a_);

    bool is_included(const Point3D& point) const final override;
    std::tuple<size_t, size_t, size_t> steps_to_cover(const UnitCell3D& unit_cell_) const final override;

    double a;
  };

  class CutoffSphere : public Cutoff
  {
  public:
    CutoffSphere(const double r_);

    bool is_included(const Point3D& point) const final override;
    std::tuple<size_t, size_t, size_t> steps_to_cover(const UnitCell3D& unit_cell_) const final override;

    double r;
  };

  class CutoffUnitVectors : public Cutoff
  {
  public:
    CutoffUnitVectors(
      const UnitCell3D& unit_cell_,
      const size_t a_max_, const size_t b_max_, const size_t c_max_);

    bool is_included(const Point3D& point) const final override;
    std::tuple<size_t, size_t, size_t> steps_to_cover(const UnitCell3D& unit_cell_) const final override;

    UnitCell3D unit_cell;
    size_t a_max;
    size_t b_max;
    size_t c_max;
  };

  // class CutoffWSCell{...}
}

#endif //LATTICEGEN_CUTOFF_H
