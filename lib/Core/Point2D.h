#ifndef LATTICEGEN_POINT2D_H_
#define LATTICEGEN_POINT2D_H_

#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>

#include <Core/ComparisonHelpers.h>

namespace Core
{
  struct Point2D
  {
    Point2D(double x_ = 0.0, double y_ = 0.0)
      : x(x_), y(y_) {}

    // operators
    Point2D &operator=(const Point2D &rhs)
    {
      this->x = rhs.x;
      this->y = rhs.y;
      return *this;
    }

    Point2D &operator+=(const Point2D &rhs)
    {
      this->x += rhs.x;
      this->y += rhs.y;
      return *this;
    }

    Point2D &operator-=(const Point2D &rhs)
    {
      this->x -= rhs.x;
      this->y -= rhs.y;
      return *this;
    }

    Point2D &operator*=(double num)
    {
      this->x *= num;
      this->y *= num;
      return *this;
    }

    double length() const
    {
      return sqrt(x * x + y * y);
    }

    double x;
    double y;
  };

  inline bool operator==(const Point2D &lhs, const Point2D &rhs)
  {
    return equalsWithTolerance(lhs.x, rhs.x) && equalsWithTolerance(lhs.y, rhs.y);
  }

  inline Point2D operator+(Point2D lhs, const Point2D &rhs)
  {
    lhs += rhs;
    return lhs;
  }

  inline Point2D operator-(Point2D lhs, const Point2D &rhs)
  {
    lhs -= rhs;
    return lhs;
  }

  inline Point2D operator*(Point2D lhs, double num)
  {
    lhs *= num;
    return lhs;
  }

  inline Point2D operator*(double num, Point2D rhs)
  {
    rhs *= num;
    return rhs;
  }

  inline double operator*(Point2D lhs, const Point2D &rhs)
  {
    return lhs.x * rhs.x + lhs.y * rhs.y;
  }

  inline bool isRectangular(const Point2D &lhs, const Point2D &rhs)
  {
    return equalsWithTolerance(lhs * rhs, 0.0);
  }

  inline void swap(Point2D &lhs, Point2D &rhs)
  {
    double tmp_x = lhs.x;
    double tmp_y = lhs.y;
    lhs.x = rhs.x;
    lhs.y = rhs.y;
    rhs.x = tmp_x;
    rhs.y = tmp_y;
  }
}

namespace std
{
  inline std::string to_string(const Core::Point2D &p)
  {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << "(" << p.x << " , "
       << p.y << " )";
    return ss.str();
  }

  inline ostream &operator<<(ostream &o, const Core::Point2D &p)
  {
    o << std::to_string(p);
    return o;
  }
}
#endif // LATTICEGEN_POINT2D_H_
