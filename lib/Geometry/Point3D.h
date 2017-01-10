#ifndef LATTICEGEN_POINT3D_H_
#define LATTICEGEN_POINT3D_H_

#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>

#include <Core/ComparisonHelpers.h>
using namespace Core;

namespace Geometry
{
  struct Point3D
  {
    Point3D(const double x_ = 0.0, const double y_ = 0.0, const double z_ = 0.0) 
      : x(x_), y(y_), z(z_)
    {}
    
    // operators
    Point3D& operator=(const Point3D& rhs)
    {	
      this->x = rhs.x;
      this->y = rhs.y;
      this->z = rhs.z;
     return *this;
    }
    Point3D& operator+=(const Point3D& rhs)
    {
      this->x += rhs.x;
      this->y += rhs.y;
      this->z += rhs.z;
      return *this;
    }
    Point3D& operator-=(const Point3D& rhs)
    {
      this->x -= rhs.x;
      this->y -= rhs.y;
      this->z -= rhs.z;
      return *this;
    }
    Point3D& operator*=(const double num)
    {
      this->x *= num;
      this->y *= num;
      this->z *= num;
      return *this;
    }
    Point3D& operator/=(const double num)
    {
      this->x /= num;
      this->y /= num;
      this->z /= num;
      return *this;
    }

    double length() const
    {
      return sqrt( x*x + y*y + z*z );
    }

	  double x;
    double y;
    double z;
  };

  typedef Point3D Vector3D;

  inline bool operator==(const Point3D& lhs, const Point3D& rhs)
  {
    return equalsWithTolerance(lhs.x, rhs.x) && 
      equalsWithTolerance(lhs.y, rhs.y) && 
      equalsWithTolerance(lhs.z, rhs.z);
  }
  inline Point3D operator+(Point3D lhs, const Point3D& rhs)
  {
    lhs += rhs;
    return lhs;
  }
  inline Point3D operator-(Point3D lhs, const Point3D& rhs)
  {
    lhs -= rhs;
    return lhs;
  }
  inline Point3D operator*(Point3D lhs, const double num)
  {
    lhs *= num;
    return lhs;
  }
  inline Point3D operator*(const double num, Point3D rhs)
  {
    rhs *= num;
    return rhs;
  }
  inline Point3D operator/(Point3D lhs, const double num)
  {
    lhs /= num;
    return lhs;
  }
  inline double operator*(Point3D lhs, const Point3D& rhs)
  {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
  }
  inline Point3D cross_product(Point3D lhs, const Point3D& rhs)
  {
    return Point3D(
        lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.z * rhs.x - lhs.x * rhs.z,
        lhs.x * rhs.y - lhs.y * rhs.x
      );
  }
  inline void swap(Point3D& lhs, Point3D& rhs)
  {
    double tmp_x = lhs.x;
    double tmp_y = lhs.y;
    double tmp_z = lhs.z;

    lhs.x = rhs.x;
    lhs.y = rhs.y;
    lhs.z = rhs.z;

    rhs.x = tmp_x;
    rhs.y = tmp_y;
    rhs.z = tmp_z;
  }
}

namespace std
{
  using namespace Geometry;
  inline std::string to_string(const Point3D& point)
  {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(std::numeric_limits< double >::max_digits10)
       << "(" << point.x << " , " << point.y << " , " << point.z << " )";
    return ss.str();
  }

  inline ostream& operator<<(ostream& o, const Geometry::Point3D& p)
  {
    o << to_string(p);
    return o;
  }
}
#endif // LATTICEGEN_POINT3D_H_
