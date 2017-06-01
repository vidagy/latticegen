#ifndef LATTICEGEN_STRUCTURECONSTANTS_H
#define LATTICEGEN_STRUCTURECONSTANTS_H

#include <complex>
#include <Geometry/UnitCell3D.h>
#include <Geometry/Shell.h>

using namespace Geometry;

namespace Physics
{
  namespace NonRelativistic
  {
    struct StructureConstantsConfig
    {
      StructureConstantsConfig(
        const double ewald_param_ = default_ewald_param,
        const double lattice_cutoff_scale_ = default_lattice_cutoff_scale,
        const double reciprocal_lattice_cutoff_ = default_reciprocal_lattice_cutoff,
        const double integral_tolerance_ = default_integral_tolerance,
        const int max_step_count_ = default_max_step_count,
        const double steps_per_unit_ = default_steps_per_unit
      ) : ewald_param(ewald_param_), lattice_cutoff_scale(lattice_cutoff_scale_),
          reciprocal_lattice_cutoff(reciprocal_lattice_cutoff_), integral_tolerance(integral_tolerance_),
          max_step_count(max_step_count_), steps_per_unit(steps_per_unit_) {}

      const double ewald_param;
      const double lattice_cutoff_scale;
      const double reciprocal_lattice_cutoff;
      const double integral_tolerance;
      const int max_step_count;
      const double steps_per_unit;

      // TODO review default values here
      constexpr static double default_ewald_param = 1.0;
      constexpr static double default_lattice_cutoff_scale = 5.0;
      constexpr static double default_reciprocal_lattice_cutoff = 5.0;
      constexpr static double default_integral_tolerance = 1e-12;
      const static int default_max_step_count = 10000;
      constexpr static double default_steps_per_unit = 200.0;
    };

    class StructureConstants
    {
    public:
      StructureConstants(
        const UnitCell3D &unit_cell_,
        const StructureConstantsConfig config_ = StructureConstantsConfig()
      );

      std::complex<double> calculate_real_space(
        unsigned int l, unsigned int m, unsigned int lprime, unsigned int mprime,
        Coordinates3D n, Coordinates3D nprime, const std::complex<double> &z
      ) const;

      std::complex<double> calculate_reciprocal_space(
        unsigned int l, int m, unsigned int lprime, int mprime,
        const Vector3D &k, const std::complex<double> &z
      ) const;

    private:
      std::complex<double> D1(unsigned int l, int m, const Vector3D &k, const std::complex<double> &z) const;

      std::complex<double> D2(unsigned int l, int m, const Vector3D &k, const std::complex<double> &z) const;

      std::complex<double> D3(unsigned int l, int m, const Vector3D &k, const std::complex<double> &z) const;

      const UnitCell3D unit_cell;
      const std::vector<Shell> direct_shells;
      const std::vector<Shell> reciprocal_shells;
      const StructureConstantsConfig config;

      friend class TestAccessor;
    };
  }
}

#endif //LATTICEGEN_STRUCTURECONSTANTS_H
