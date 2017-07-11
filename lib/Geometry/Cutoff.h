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

      const long x_min;
      const long y_min;
      const long z_min;
      const long x_max;
      const long y_max;
      const long z_max;
    };

    virtual StepsToCover steps_to_cover(const Cell3D &cell_) const = 0;
  };

  class CutoffCube : public Cutoff
  {
  public:
    CutoffCube(double a_);

    bool is_included(const Point3D& point) const final override;

    StepsToCover steps_to_cover(const Cell3D &cell_) const final override;

    const double a;
  };

  class CutoffSphere : public Cutoff
  {
  public:
    CutoffSphere(double r_);

    bool is_included(const Point3D& point) const final override;

    StepsToCover steps_to_cover(const Cell3D &cell_) const final override;

    const double r;
  };

  class CutoffUnitVectors : public Cutoff
  {
  public:
    CutoffUnitVectors(
      const Cell3D &cell_,
      size_t a_max_, size_t b_max_, size_t c_max_, bool positive_only_ = false);

    bool is_included(const Point3D& point) const final override;

    StepsToCover steps_to_cover(const Cell3D &cell_) const final override;

    const Cell3D cell;
    const size_t a_max;
    const size_t b_max;
    const size_t c_max;
    const bool positive_only;
  };

  class CutoffWSCell : public Cutoff
  {
  public:
    CutoffWSCell(const Cell3D &cell_);

    bool is_included(const Point3D& point) const final override;

    StepsToCover steps_to_cover(const Cell3D &cell_) const final override;

    const Cell3D cell;
    // face_points are pairs of unit vectors and lengths: unit vector is the perpendicular direction from the plane
    // length is the distance of the plane from the origin
    const std::vector<std::pair<Point3D, double>> face_points;
    const double r_mt; /// muffin tin radius
    const double r_bs; /// radius of bounding sphere
  };
}

#endif //LATTICEGEN_CUTOFF_H
