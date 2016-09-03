#ifndef LATTICEGEN_GEOMETRY_H_
#define LATTICEGEN_GEOMETRY_H_

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
      Point2D(const double x = 0.0, const double y = 0.0) : x_(x), y_(y)
      { }
      
      // operators
      Point2D& operator=(const Point2D& rhs)
      {	
        this->x_ = rhs.x_;
        this->y_ = rhs.y_;
	     return *this;
	    }
      Point2D& operator+=(const Point2D& rhs)
      {
        this->x_ += rhs.x_;
        this->y_ += rhs.y_;
        return *this;
      }
      Point2D& operator-=(const Point2D& rhs)
      {
        this->x_ -= rhs.x_;
        this->y_ -= rhs.y_;
        return *this;
      }
      Point2D& operator*=(const double num)
      {
        this->x_ *= num;
        this->y_ *= num;
        return *this;
      }

      double getLength() const
      {
        return sqrt( x_*x_ + y_*y_ );
      }

      std::string toString() const 
      {
        std::stringstream ss;
        ss << std::scientific << std::setprecision(std::numeric_limits< double >::max_digits10) << "(" << x_ << " , " << y_ << " )";
        return ss.str();
      }

  	  double x_;
  	  double y_;
    };

    inline bool operator==(const Point2D& lhs, const Point2D& rhs)
    {
      return equalsWithTolerance(lhs.x_, rhs.x_) && equalsWithTolerance(lhs.y_, rhs.y_);
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
      return lhs.x_ * rhs.x_ + lhs.y_ * rhs.y_;
    }
    inline bool isRectangular(const Point2D& lhs, const Point2D& rhs)
    {
      return equalsWithTolerance(lhs * rhs, 0.0);
    }
    inline void swap(Point2D& lhs, Point2D& rhs)
    {
      double tmp_x = lhs.x_;
      double tmp_y = lhs.y_;
      lhs.x_ = rhs.x_;
      lhs.y_ = rhs.y_;
      rhs.x_ = tmp_x;
      rhs.y_ = tmp_y;
    }
  }
}

namespace std
{
  ostream& operator<<(ostream& o, const Core::Geometry::Point2D& p);
}
#endif // LATTICEGEN_GEOMETRY_H_
