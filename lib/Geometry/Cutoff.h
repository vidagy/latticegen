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

    struct StepsToCover
    {
      StepsToCover(long x_min_, long y_min_, long z_min_, long x_max_, long y_max_, long z_max_)
        : x_min(x_min_), y_min(y_min_), z_min(z_min_)
        , x_max(x_max_), y_max(y_max_), z_max(z_max_)
      {}
      long x_min;
      long y_min;
      long z_min;
      long x_max;
      long y_max;
      long z_max;
    };
    virtual StepsToCover steps_to_cover(const UnitCell3D& unit_cell_) const = 0;
  };

  class CutoffCube : public Cutoff
  {
  public:
    CutoffCube(const double a_);

    bool is_included(const Point3D& point) const final override;
    StepsToCover steps_to_cover(const UnitCell3D& unit_cell_) const final override;

    double a;
  };

  class CutoffSphere : public Cutoff
  {
  public:
    CutoffSphere(const double r_);

    bool is_included(const Point3D& point) const final override;
    StepsToCover steps_to_cover(const UnitCell3D& unit_cell_) const final override;

    double r;
  };

  class CutoffUnitVectors : public Cutoff
  {
  public:
    CutoffUnitVectors(
      const UnitCell3D& unit_cell_,
      const size_t a_max_, const size_t b_max_, const size_t c_max_, const bool positive_only_ = false);

    bool is_included(const Point3D& point) const final override;
    StepsToCover steps_to_cover(const UnitCell3D& unit_cell_) const final override;

    UnitCell3D unit_cell;
    size_t a_max;
    size_t b_max;
    size_t c_max;
    bool positive_only;
  };

  class CutoffWSCell : public Cutoff
  {
  public:
    CutoffWSCell(const UnitCell3D& unit_cell_);

    bool is_included(const Point3D& point) const final override;
    StepsToCover steps_to_cover(const UnitCell3D& unit_cell_) const final override;

  private:
    UnitCell3D unit_cell;
    std::vector<Point3D> neighbors;
  };
}

#endif //LATTICEGEN_CUTOFF_H
