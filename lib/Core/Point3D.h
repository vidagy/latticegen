#ifndef LATTICEGEN_POINT3D_H_
#define LATTICEGEN_POINT3D_H_

#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

#include <Core/ComparisonHelpers.h>

namespace Core
{
  struct Point3D
  {
    Point3D(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0)
      : x(x_), y(y_), z(z_)
    {}

    static Point3D create_polar(double r, double theta, double phi)
    {
      return Point3D{r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)};
    }
    
    // operators
    Point3D &operator=(const Point3D &rhs) = default;

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

    Point3D &operator*=(double num)
    {
      this->x *= num;
      this->y *= num;
      this->z *= num;
      return *this;
    }

    Point3D &operator/=(double num)
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

    double length2() const
    {
      return x * x + y * y + z * z;
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

  inline Point3D operator+(const Point3D &lhs, const Point3D &rhs)
  {
    return Point3D{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
  }

  inline Point3D operator-(const Point3D &lhs, const Point3D &rhs)
  {
    return Point3D{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
  }
  inline Point3D operator-(Point3D lhs)
  {
    lhs.x = -lhs.x;
    lhs.y = -lhs.y;
    lhs.z = -lhs.z;
    return lhs;
  }

  inline Point3D operator*(Point3D lhs, double num)
  {
    lhs *= num;
    return lhs;
  }

  inline Point3D operator*(double num, Point3D rhs)
  {
    rhs *= num;
    return rhs;
  }

  inline Point3D operator/(Point3D lhs, double num)
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

  struct Point3DComparator
  {
    explicit Point3DComparator(double abs_tol, double rel_tol)
      : x_abs_tol(abs_tol), y_abs_tol(abs_tol), z_abs_tol(abs_tol), x_rel_tol(rel_tol), y_rel_tol(rel_tol),
        z_rel_tol(rel_tol) {}

    Point3DComparator(
      double x_abs_tol_, double y_abs_tol_, double z_abs_tol_,
      double x_rel_tol_, double y_rel_tol_, double z_rel_tol_)
      : x_abs_tol(x_abs_tol_), y_abs_tol(y_abs_tol_), z_abs_tol(z_abs_tol_), x_rel_tol(x_rel_tol_),
        y_rel_tol(y_rel_tol_), z_rel_tol(z_rel_tol_) {}

    bool isEqual(const Point3D &lhs, const Point3D &rhs) const
    {
      return equalsWithTolerance(lhs.x, rhs.x, x_abs_tol, x_rel_tol) &&
             equalsWithTolerance(lhs.y, rhs.y, y_abs_tol, y_rel_tol) &&
             equalsWithTolerance(lhs.z, rhs.z, z_abs_tol, z_rel_tol);
    }

    const double x_abs_tol;
    const double y_abs_tol;
    const double z_abs_tol;
    const double x_rel_tol;
    const double y_rel_tol;
    const double z_rel_tol;
  };
}

namespace std
{
  using namespace Core;
  inline std::string to_string(const Point3D& point)
  {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(std::numeric_limits< double >::max_digits10)
       << "(" << std::setw(24) << point.x << " , " << std::setw(24) << point.y << " , " << std::setw(24) << point.z
       << " )";
    return ss.str();
  }

  inline ostream &operator<<(ostream &o, const Core::Point3D &p)
  {
    o << to_string(p);
    return o;
  }
}
#endif // LATTICEGEN_POINT3D_H_
