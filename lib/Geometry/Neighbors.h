#ifndef LATTICEGEN_NEIGHBORS_H_
#define LATTICEGEN_NEIGHBORS_H_

#include "Point3D.h"

namespace Geometry
{
  class Cutoff
  {
    virtual bool is_included(const Point3D& point) const = 0;
  };

  class CutoffCube : public Cutoff
  {
  public:
    CutoffCube(const double a_);

    bool is_included(const Point3D& point) const override;
  private:
    double a;
  };

  class CutoffSphere : public Cutoff
  {
  public:
    CutoffSphere(const double r_);

    bool is_included(const Point3D& point) const override;
  private:
    double r;
  };
}

#endif // LATTICEGEN_NEIGHBORS_H_