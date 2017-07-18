#ifndef LATTICEGEN_POLYHEDRALQUADRATURE_H
#define LATTICEGEN_POLYHEDRALQUADRATURE_H

#include <Core/Point3D.h>
#include <vector>

namespace Math
{
  using namespace Core;

  typedef std::vector<std::pair<Point3D, double>> Quadrature;

  /// @brief Generating mesh on the surface of the unit sphere by repeated subdivision and remapping of face of
  /// an icosahedron
  class IcosahedralQuadrature
  {
  public:
    /// @brief n is the number of subdivisions. The total number of grid points are 20 * 2^(2*n). Grid points are
    /// located in the center of each triangle faces and weighted by the area of surrounding triangle
    static Quadrature generate(unsigned int n);
  };

  /// @brief Generating mesh on the surface of the unit sphere by repeated subdivision and remapping of face of
  /// a cube
  class OctahedralQuadrature
  {
  public:
    /// @brief n is the number of subdivisions. The total number of grid points are 6 * 2^(2*n). Grid points are
    /// located in the center of each triangle faces and weighted by the area of surrounding rectangle
    static Quadrature generate(unsigned int n);
  };
}
#endif //LATTICEGEN_POLYHEDRALQUADRATURE_H
