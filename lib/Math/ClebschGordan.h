#ifndef LATTICEGEN_CLEBSCHGORDAN_H
#define LATTICEGEN_CLEBSCHGORDAN_H


namespace Math
{
  class ClebschGordan
  {
  public:
    /// @brief It calculates the < j1, m1 ; j2, m2 | J, M > integral as defined on wikipedia
    static double calculate(double j1, double m1, double j2, double m2, double J, double M);
  };

  class Gaunt
  {
  public:
    /// @brief This is the integral of Y{l1,m1}*Y{l2,m2}*[Y{L,M}]^* dOmega
    static double calculate(double l1, double m1, double l2, double m2, double L, double M);
  };
}


#endif //LATTICEGEN_CLEBSCHGORDAN_H
