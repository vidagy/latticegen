#include "Point3D.h"

namespace std
{
  ostream& operator<<(ostream& o, const Geometry::Point3D& p)
  {
    o << p.toString();
    return o;
  }
} 
