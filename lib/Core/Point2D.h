#ifndef LATTICEGEN_GEOMETRY_H_
#define LATTICEGEN_GEOMETRY_H_

#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace Core
{
  namespace Geometry
  {
    struct Point2D
    {
  	 Point2D(const double& x, const double& y) : x_(x), y_(y)
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
      Point2D& operator*=(const double& num)
      {
        this->x_ *= num;
        this->y_ *= num;
        return *this;
      }
      Point2D& operator*=(const Point2D& rhs)
      {
        this->x_ *= rhs.x_;
        this->y_ *= rhs.y_;
        return *this;
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
    inline Point2D operator*(Point2D lhs, const double& num)
    {
      lhs *= num;
      return lhs;
    }
    inline Point2D operator*(const double& num, Point2D rhs)
    {
      rhs *= num;
      return rhs;
    }
    inline Point2D operator*(Point2D lhs, const Point2D& rhs)
    {
      lhs *= rhs;
      return lhs;
    }
  } // namespace Core
}
#endif // LATTICEGEN_GEOMETRY_H_
