#include <boost/math/tools/roots.hpp>
#include <Math/CommonFunctions.h>
#include <Math/SphericalHarmonics.h>
#include "ShapeFunctions.h"
#include "Cutoff.h"

using namespace Geometry;
using namespace Math;
using namespace std::complex_literals;

namespace
{
  class ShapeFunctionsCalculator
  {
  public:
    ShapeFunctionsCalculator(const UnitCell3D unit_cell_, const ShapeFunctionsConfig &config_ = ShapeFunctionsConfig())
      : unit_cell(unit_cell_), cutoff(unit_cell), config(config_), mesh(generate_mesh_and_bounding_r(config, cutoff)) {}

  private:
    std::vector<std::pair<Point3D, double>>
    generate_mesh_and_bounding_r(const ShapeFunctionsConfig &config, const CutoffWSCell &cutoff) const
    {
      auto mesh = generate_mesh(config);
      auto mesh_and_bounding_r = std::vector<std::pair<Point3D, double>>();
      mesh_and_bounding_r.reserve(mesh.size());
      std::transform(mesh.begin(), mesh.end(), std::back_inserter(mesh_and_bounding_r),
                     [&cutoff, this](const Point3D &point)
                     {
                       auto bounding_r = this->get_bounding_r(point, cutoff);
                       return std::make_pair(point, bounding_r);
                     }
      );
      std::sort(mesh_and_bounding_r.begin(), mesh_and_bounding_r.end(),
                [](const std::pair<Point3D, double> &lhs, const std::pair<Point3D, double> &rhs)
                {
                  return lhs.second < rhs.second;
                }
      );
      return mesh_and_bounding_r;
    }

    std::vector<Point3D> generate_mesh(const ShapeFunctionsConfig &config) const
    {
      auto mesh = std::vector<Point3D>();
      mesh.reserve((unsigned long) (config.theta_res * config.phi_res));
      for (auto i = 0; i < config.theta_res; ++i) {
        auto theta = i * pi / (config.theta_res - 1);
        auto phi_res = 1 + std::lround(sin(theta) * config.phi_res);
        for (auto j = 0l; j < phi_res; ++j) {
          auto phi = (j + 0.5) * 2.0 * pi / phi_res;
          mesh.push_back(Point3D::create_polar(1.0, theta, phi));
        }
      }
      mesh.shrink_to_fit();
      return mesh;
    }

    double get_bounding_r(const Point3D &point, const CutoffWSCell &cutoff) const
    {
      auto objective_function = [&cutoff, &point](const double scale)
      {
        return cutoff.is_included(scale * point) ? 1.0 : -1.0;
      };
      // epsilons added/subtracted to make sure that bracketing is OK
      auto lower_guess = cutoff.r_in_for_sure - std::numeric_limits<double>::epsilon();
      auto upper_guess = cutoff.r_out_for_sure + std::numeric_limits<double>::epsilon();

      using namespace boost::math::tools;
      uintmax_t max_iter = (uintmax_t) config.bracketing_max_iter;
      auto res = bisect(objective_function, lower_guess, upper_guess, eps_tolerance<double>(12), max_iter);
      if (max_iter >= config.bracketing_max_iter)
        THROW_LOGIC_ERROR(
          "iter_num = " + std::to_string(max_iter) + " while max_iter = " + std::to_string(config.bracketing_max_iter)
        );

      return (res.first + res.second) / 2.0;
    }

  public:
    lm_vector<std::vector<std::complex<double>>>
    calculate(const unsigned int l_max, const std::vector<double> &r_points) const
    {
      auto is_sorted = std::is_sorted(r_points.begin(), r_points.end(),
                                      [](const Point3D &lhs, const Point3D &rhs)
                                      {
                                        return lhs.length() < rhs.length();
                                      }
      );
      if (!is_sorted)
        THROW_INVALID_ARGUMENT("r_points are not sorted");

      auto res = lm_vector<std::vector<std::complex<double>>>(l_max);
      for (auto l = 0u; l <= l_max; ++l) {
        for (auto m = -((int) l); m <= ((int) l); ++m) {
          auto spherical_harmonics_values = get_spherical_harmonics_values(l, m);
          res.at(l, m) = integrate(l, m, spherical_harmonics_values, r_points);
        }
      }
      return res;
    }

    const UnitCell3D unit_cell;
    const CutoffWSCell cutoff;
    const ShapeFunctionsConfig config;
    // note that mesh is ordered in increasing order of entry.second
    const std::vector<std::pair<Point3D, double>> mesh;

  private:
    std::vector<std::complex<double>> get_spherical_harmonics_values(const unsigned int l, const int m) const
    {
      auto res = std::vector<std::complex<double>>();
      res.reserve(mesh.size());
      std::transform(mesh.begin(), mesh.end(), std::back_inserter(res),
                     [l, m](const std::pair<Point3D, double> &p)
                     {
                       return std::conj(Complex::spherical_harmonic(l, m, p.first));
                     });
      return res;
    }

    std::vector<std::complex<double>> integrate(
      const unsigned int l, const int m,
      const std::vector<std::complex<double>> &spherical_harmonics_values,
      const std::vector<double> &r_points
    ) const
    {
      if (r_points.empty())
        return std::vector<std::complex<double>>();

      auto res = std::vector<std::complex<double>>(r_points.size(), 0.0);
      auto integral = 0.0i;

      auto itMesh = mesh.crbegin();
      const auto itMeshEnd = mesh.crend();
      auto itSH = spherical_harmonics_values.crbegin();
      const auto itSHEnd = spherical_harmonics_values.crend();

      auto itR = r_points.crbegin();
      const auto itREnd = r_points.crend();
      auto itRes = res.rbegin();
      const auto itResEnd = res.rend();
      for (; itR != itREnd && itRes != itResEnd; ++itR, ++itRes) {
        auto r = *itR;
        if (greaterEqualsWithTolerance(r, cutoff.r_out_for_sure)) {
          // we are too far, nothing is in the WS cell
          *itRes = 0.0;
        } else if (lessEqualsWithTolerance(r, cutoff.r_in_for_sure)) {
          // we are too close, the whole sphere is in the WS cell
          if (l == 0u && m == 0)
            *itRes = 2.0 * sqrt(pi);
          else
            *itRes = 0.0;
        } else {
          // we are in between
          while (itMesh != itMeshEnd && itSH != itSHEnd && itMesh->second >= r) {
            integral += *itSH;
            ++itMesh;
            ++itSH;
          }
          *itRes = integral / static_cast<double>(mesh.size()) * 4.0 * pi;
        }
      }
      if (itR != r_points.crend() || itRes != itResEnd)
        THROW_LOGIC_ERROR("r_points and res does not have the same length");

      return res;
    }
  };
}

ShapeFunctions::ShapeFunctions(const UnitCell3D &unit_cell_, unsigned int l_max, const std::vector<double> &r_points_,
                               const ShapeFunctionsConfig &config_)
  : unit_cell(unit_cell_), shape_functions(ShapeFunctionsCalculator(unit_cell, config_).calculate(l_max, r_points_)),
    r_points(r_points_) {}
