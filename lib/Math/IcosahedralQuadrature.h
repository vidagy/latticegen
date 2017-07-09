#ifndef LATTICEGEN_ICOSAHEDRALQUADRATURE_H
#define LATTICEGEN_ICOSAHEDRALQUADRATURE_H

#include <Core/Point3D.h>
#include <vector>

namespace Math
{
  using namespace Core;

  /// @brief Generating mesh on the surface of the unit sphere by repeated subdivision and remapping of face of
  /// icosahedron
  class IcosahedralQuadrature
  {
  public:
    /// @brief n is the number of subdivisions. The total number of grid points are 20 * 2^(2*n). Grid points are
    /// located in the center of each triangle faces and weighted by the area of surrounding triangle
    static std::vector<std::pair<Point3D, double>> generate(unsigned int n);
  };
}
#endif //LATTICEGEN_ICOSAHEDRALQUADRATURE_H
