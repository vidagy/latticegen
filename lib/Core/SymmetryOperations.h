#ifndef LATTICEGEN_SYMMETRYOPERATIONS_H_
#define LATTICEGEN_SYMMETRYOPERATIONS_H_

#include "Point3D.h"
#include <array>

namespace Core
{
  namespace Geometry
  {
    typedef Point3D Vector3D;

    class Rotation
    {
    public:
      typedef std::array< std::array<double,3>, 3> RotationMatrix;

      Rotation(const Vector3D& rotation_vector);
      Vector3D operator()(const Vector3D& vector) const;
    private: 
      RotationMatrix rotation_matrix;
    };

    inline Vector3D operator*(const Rotation& rotation, const Vector3D& vector)
    {
      return rotation(vector);
    }
  }
}
 
#endif // LATTICEGEN_SYMMETRYOPERATIONS_H_