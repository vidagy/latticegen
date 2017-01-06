#include "Point2D.h"

namespace std
{
  ostream& operator<<(ostream& o, const Geometry::Point2D& p)
  {
    o << p.toString();
    return o;
  }
} 
