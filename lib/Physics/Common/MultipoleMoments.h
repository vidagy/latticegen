#ifndef LATTICEGEN_MULTIPOLEMOMENTS_H
#define LATTICEGEN_MULTIPOLEMOMENTS_H

#include "ChargeDensity.h"

namespace Physics
{
  namespace Common
  {
    using namespace Geometry;

    class MultipoleMoments
    {
    public:
      ///@brief Zabloudil et al (19.21)
      MultipoleMoments(const ShapeTruncatedChargeDensity &charge_density);

      lm_vector<std::complex<double>> moments;
    };

  }
}


#endif //LATTICEGEN_MULTIPOLEMOMENTS_H
