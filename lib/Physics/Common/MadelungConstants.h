#ifndef LATTICEGEN_MADELUNGCONSTANTS_H
#define LATTICEGEN_MADELUNGCONSTANTS_H

#include <Geometry/UnitCell3D.h>
#include <complex>
#include <Geometry/Shell.h>
#include <Core/lm_vector.h>

using namespace Geometry;
using namespace Core;

namespace Physics
{
  namespace Common
  {
    /// @brief These constants are not summed up for the lattice, these are local quantities
    class RealMadelungConstants
    {
    public:
      RealMadelungConstants(const UnitCell3D &unit_cell_) : unit_cell(unit_cell_) {}

      ///@brief Zabloudil et al (19.28)
      std::complex<double> calculate(
        unsigned int l, int m, unsigned int lprime, int mprime, Coordinates3D n, Coordinates3D nprime
      ) const;

      ///@brief Zabloudil et al (19.36)
      std::complex<double> calculateReduced(unsigned int l, int m, Coordinates3D n, Coordinates3D nprime) const;

      const UnitCell3D unit_cell;
    };

    struct MadelungConstantsConfig
    {
      MadelungConstantsConfig(
        const double ewald_param_ = default_ewald_param,
        const double lattice_cutoff_scale_ = default_lattice_cutoff_scale,
        const double reciprocal_lattice_cutoff_ = default_reciprocal_lattice_cutoff
      ) : ewald_param(ewald_param_), lattice_cutoff_scale(lattice_cutoff_scale_),
          reciprocal_lattice_cutoff(reciprocal_lattice_cutoff_) {}

      const double ewald_param;
      const double lattice_cutoff_scale;
      const double reciprocal_lattice_cutoff;

      // TODO review default values here
      constexpr static double default_ewald_param = 1.0;
      constexpr static double default_lattice_cutoff_scale = 5.0;
      constexpr static double default_reciprocal_lattice_cutoff = 5.0;
    };

    /// @brief Madelung constants as defined in Zabloudil et al. (19.89)
    class ReducedMadelungConstants
    {
    public:
      ReducedMadelungConstants(
        const UnitCell3D &unit_cell_,
        const unsigned int l_max,
        const MadelungConstantsConfig &config_ = MadelungConstantsConfig()
      );

      const UnitCell3D unit_cell;
      const lm_vector<std::complex<double>> madelung_constants;
    };
  }
}
#endif //LATTICEGEN_MADELUNGCONSTANTS_H
