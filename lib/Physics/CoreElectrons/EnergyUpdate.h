#ifndef LATTICEGEN_ENERGYUPDATE_H
#define LATTICEGEN_ENERGYUPDATE_H

namespace Physics
{
  namespace CoreElectrons
  {
    struct EnergyUpdate
    {
      EnergyUpdate(int required_number_of_nodes_, double energy_tolerance_)
        : required_number_of_nodes(required_number_of_nodes_), energy_tolerance(energy_tolerance_), lower(0.0),
          upper(0.0), lower_set(false), upper_set(false) {}

      double coarse(int number_of_nodes, double energy_)
      {
        double energy = energy_;
        if (number_of_nodes > required_number_of_nodes) {
          // energy is too large
          // set or update upper energy bound
          if (!upper_set) {
            upper_set = true;
            upper = energy;
          }
          if (energy < upper)
            upper = energy;
          if (lower_set)
            // if we have a lower bound then bisect
            energy = (upper + lower) / 2.0;
          else
            // else guess 10% less energy (energy is negative, less is more :) )
            energy = 1.1 * upper;
        } else {
          // energy is too small
          // set or update lower energy bound
          if (!lower_set) {
            lower_set = true;
            lower = energy;
          }
          if (energy > lower)
            lower = energy;
          if (upper_set)
            // if we have an upper bound then bisect
            energy = (upper + lower) / 2.0;
          else
            // else guess 10% more energy
            energy = 0.9 * lower;
        }
        return energy;
      }

      std::pair<double, bool> fine(double energy, double norm, double new_R, double new_dR_dr, double old_dR_dr)
      {
        double energy_diff = new_R * (new_dR_dr - old_dR_dr) / (2.0 * norm);

        if (lower_set && (energy + energy_diff < lower))
          return std::make_pair((energy + lower) / 2.0, false);
        else if (energy + energy_diff > upper)
          return std::make_pair((energy + upper) / 2.0, false);
        else if (fabs(energy_diff / energy) > energy_tolerance)
          return std::make_pair(energy + energy_diff, false);
        else
          return std::make_pair(energy, true);
      }

      const int required_number_of_nodes;
      const double energy_tolerance;

    private:
      // TODO if C++17 is available rework this to optional
      double lower;
      double upper;
      bool lower_set;
      bool upper_set;
    };
  }
}

#endif //LATTICEGEN_ENERGYUPDATE_H
