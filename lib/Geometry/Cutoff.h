#ifndef LATTICEGEN_CUTOFF_H_
#define LATTICEGEN_CUTOFF_H_

#include "Point3D.h"
#include "UnitCell3D.h"

#include <vector>
#include <tuple>
namespace Geometry
{
  class Cutoff
  {
  public:
    Cutoff(const UnitCell3D& unit_cell_)
      : unit_cell(unit_cell_)
    {}

    virtual bool is_included(const Point3D& point) const = 0;
    virtual std::tuple<size_t, size_t, size_t> steps_to_cover() const = 0;

    UnitCell3D unit_cell;
  };

  class CutoffCube : public Cutoff
  {
  public:
    CutoffCube(const UnitCell3D& unit_cell_, const double a_);

    bool is_included(const Point3D& point) const override;
    std::tuple<size_t, size_t, size_t> steps_to_cover() const override;

    double a;
  };

  class CutoffSphere : public Cutoff
  {
  public:
    CutoffSphere(const UnitCell3D& unit_cell_, const double r_);

    bool is_included(const Point3D& point) const override;
    std::tuple<size_t, size_t, size_t> steps_to_cover() const override;
  
    double r;
  };

  class CutoffUnitVectors : public Cutoff
  {
  public:
    CutoffUnitVectors(
      const UnitCell3D& unit_cell_,
      const size_t a_max_, const size_t b_max_, const size_t c_max_);

    bool is_included(const Point3D& point) const override;
    std::tuple<size_t, size_t, size_t> steps_to_cover() const override;
  
    size_t a_max;
    size_t b_max;
    size_t c_max;
  };
}

#endif // LATTICEGEN_CUTOFF_H_