#ifndef LATTICEGEN_POINT2D_H_
#define LATTICEGEN_POINT2D_H_

#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>

#include "ComparisonHelpers.h"

namespace Core
{
  namespace Geometry
  {
    struct Point2D
    {
      Point2D(const double x_ = 0.0, const double y_ = 0.0) 
        : x(x_), y(y_)
      {}
      
      // operators
      Point2D& operator=(const Point2D& rhs)
      {	
        this->x = rhs.x;
        this->y = rhs.y;
	     return *this;
	    }
      Point2D& operator+=(const Point2D& rhs)
      {
        this->x += rhs.x;
        this->y += rhs.y;
        return *this;
      }
      Point2D& operator-=(const Point2D& rhs)
      {
        this->x -= rhs.x;
        this->y -= rhs.y;
        return *this;
      }
      Point2D& operator*=(const double num)
      {
        this->x *= num;
        this->y *= num;
        return *this;
      }

      double getLength() const
      {
        return sqrt( x*x + y*y );
      }

      std::string toString() const 
      {
        std::stringstream ss;
        ss << std::scientific << std::setprecision(std::numeric_limits< double >::max_digits10) << "(" << x << " , " << y << " )";
        return ss.str();
      }

  	  double x;
  	  double y;
    };

    inline bool operator==(const Point2D& lhs, const Point2D& rhs)
    {
      return equalsWithTolerance(lhs.x, rhs.x) && equalsWithTolerance(lhs.y, rhs.y);
    }
    inline Point2D operator+(Point2D lhs, const Point2D& rhs)
    {
      lhs += rhs;
      return lhs;
    }
    inline Point2D operator-(Point2D lhs, const Point2D& rhs)
    {
      lhs -= rhs;
      return lhs;
    }
    inline Point2D operator*(Point2D lhs, const double num)
    {
      lhs *= num;
      return lhs;
    }
    inline Point2D operator*(const double num, Point2D rhs)
    {
      rhs *= num;
      return rhs;
    }
    inline double operator*(Point2D lhs, const Point2D& rhs)
    {
      return lhs.x * rhs.x + lhs.y * rhs.y;
    }
    inline bool isRectangular(const Point2D& lhs, const Point2D& rhs)
    {
      return equalsWithTolerance(lhs * rhs, 0.0);
    }
    inline void swap(Point2D& lhs, Point2D& rhs)
    {
      double tmp_x = lhs.x;
      double tmp_y = lhs.y;
      lhs.x = rhs.x;
      lhs.y = rhs.y;
      rhs.x = tmp_x;
      rhs.y = tmp_y;
    }
  }
}

namespace std
{
  ostream& operator<<(ostream& o, const Core::Geometry::Point2D& p);
}
#endif // LATTICEGEN_POINT2D_H_
