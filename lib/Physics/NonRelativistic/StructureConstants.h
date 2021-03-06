#ifndef LATTICEGEN_STRUCTURECONSTANTS_H
#define LATTICEGEN_STRUCTURECONSTANTS_H

#include <complex>
#include <Geometry/UnitCell3D.h>
#include <Geometry/Shell.h>
#include <map>
#include <utility>

using namespace Geometry;

namespace Physics
{
  namespace NonRelativistic
  {
    class RealStructureConstants
    {
    public:
      explicit RealStructureConstants(UnitCell3D unit_cell_) : unit_cell(std::move(unit_cell_)) {}

      std::complex<double> calculate(
        unsigned int l, int m, unsigned int lprime, int mprime,
        Coordinates3D n, Coordinates3D nprime, const std::complex<double> &z
      ) const;

      const UnitCell3D unit_cell;
    };

    struct StructureConstantsConfig
    {
      explicit StructureConstantsConfig(
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

    struct IntegralCache
    {
      IntegralCache(
        const std::shared_ptr<const std::vector<Shell>> &direct_shells,
        const std::shared_ptr<const StructureConstantsConfig> &config,
        unsigned int l_max_,
        const std::complex<double> &z
      );

      std::complex<double> get(unsigned int l, unsigned int n) const;

      const unsigned int n_max;
      const unsigned int l_max;
      const std::vector<std::complex<double>> cache;
    };

    class ReciprocalStructureConstantsCalculator
    {
    public:
      std::complex<double>
      calculate(unsigned int l, int m, unsigned int lprime, int mprime, const Point3DCRef &k) const;

    private:
      ReciprocalStructureConstantsCalculator(
        std::shared_ptr<const UnitCell3D> unit_cell_,
        std::shared_ptr<const std::vector<Shell>> direct_shells_,
        std::shared_ptr<const std::vector<Shell>> reciprocal_shells_,
        std::shared_ptr<const StructureConstantsConfig> config_,
        const std::complex<double> &z_,
        IntegralCache integral_cache_
      ) : unit_cell(std::move(unit_cell_)), direct_shells(std::move(direct_shells_)), reciprocal_shells(
        std::move(reciprocal_shells_)),
          config(std::move(config_)), z(z_), integral_cache(std::move(integral_cache_)) {}

      std::complex<double> D1(unsigned int l, int m, const Point3DCRef &k) const;

      std::complex<double> D2(unsigned int l, int m, const Point3DCRef &k) const;

      std::complex<double> D3(unsigned int l, int m) const;

      const std::shared_ptr<const UnitCell3D> unit_cell;
      const std::shared_ptr<const std::vector<Shell>> direct_shells;
      const std::shared_ptr<const std::vector<Shell>> reciprocal_shells;
      const std::shared_ptr<const StructureConstantsConfig> config;
      const std::complex<double> z;
      const IntegralCache integral_cache;

      friend class ReciprocalStructureConstants;

      friend class TestAccessor;
    };

    class ReciprocalStructureConstants
    {
    public:
      explicit ReciprocalStructureConstants(
        const UnitCell3D &unit_cell_,
        const StructureConstantsConfig &config_ = StructureConstantsConfig()
      );

      ReciprocalStructureConstantsCalculator get_calculator(unsigned int l_max, const std::complex<double> &z) const;

      const std::shared_ptr<const UnitCell3D> unit_cell;
      const std::shared_ptr<const std::vector<Shell>> direct_shells;
      const std::shared_ptr<const std::vector<Shell>> reciprocal_shells;
      const std::shared_ptr<const StructureConstantsConfig> config;
    };
  }
}

#endif //LATTICEGEN_STRUCTURECONSTANTS_H
