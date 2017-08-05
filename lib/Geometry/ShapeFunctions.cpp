#include <boost/math/tools/roots.hpp>
#include <Math/CommonFunctions.h>
#include <Math/SphericalHarmonics.h>
#include <Math/PolyhedralQuadrature.h>
#include <Math/WignerDMatrix.h>
#include "ShapeFunctions.h"
#include "Cutoff.h"
#include "SymmetryTransformationFactory.h"

using namespace Geometry;
using namespace Math;
using namespace std::complex_literals;

namespace
{
  class ShapeFunctionsCalculator
  {
  public:
    ShapeFunctionsCalculator(const UnitCell3D unit_cell_, const ShapeFunctionsConfig &config_ = ShapeFunctionsConfig())
      : cutoff(unit_cell_), config(config_), mesh(generate_mesh_and_bounding_r(config, cutoff)) {}

  private:
    struct MeshPoint
    {
      MeshPoint(const Point3D &point_, const double weight_, const double r_WS_)
        : point(point_), weight(weight_), r_WS(r_WS_) {}

      Point3D point;
      double weight;
      double r_WS;
    };

    std::vector<MeshPoint>
    generate_mesh_and_bounding_r(const ShapeFunctionsConfig &config, const CutoffWSCell &cutoff) const
    {
      auto quadrature = OctahedralQuadrature::generate(config.quadrature_order);
      auto mesh_points = std::vector<MeshPoint>();
      mesh_points.reserve(quadrature.size());
      std::transform(
        quadrature.begin(), quadrature.end(), std::back_inserter(mesh_points),
        [&cutoff, this](const std::pair<Point3D, double> &point_and_weight)
        {
          auto bounding_r = this->get_bounding_r(point_and_weight.first, cutoff);
          return MeshPoint(point_and_weight.first, point_and_weight.second, bounding_r);
        }
      );
      std::sort(
        mesh_points.begin(), mesh_points.end(),
        [](const MeshPoint &lhs, const MeshPoint &rhs)
        {
          return lhs.r_WS < rhs.r_WS;
        }
      );
      return mesh_points;
    }

    double get_bounding_r(const Point3D &point, const CutoffWSCell &cutoff) const
    {
      auto objective_function = [&cutoff, &point](const double scale)
      {
        return cutoff.is_included(scale * point) ? 1.0 : -1.0;
      };
      // epsilons added/subtracted to make sure that bracketing is OK
      auto lower_guess = cutoff.r_mt - 1e8 * std::numeric_limits<double>::epsilon();
      auto upper_guess = cutoff.r_bs + 1e8 * std::numeric_limits<double>::epsilon();

      using namespace boost::math::tools;
      uintmax_t max_iter = (uintmax_t) config.bracketing_max_iter;
      auto res = bisect(
        objective_function, lower_guess, upper_guess, eps_tolerance<double>(config.default_r_ws_bits), max_iter
      );
      if (max_iter >= config.bracketing_max_iter)
        THROW_LOGIC_ERROR(
          "iter_num = " + std::to_string(max_iter) + " while max_iter = " + std::to_string(config.bracketing_max_iter)
        );

      return (res.first + res.second) / 2.0;
    }

    // workaround because ref to boolean vector element is deleted
    enum class ZeroBySymmetries
    {
      True, False
    };

  public:
    lm_vector<std::vector<std::complex<double>>> calculate(const unsigned int l_max) const
    {
      auto zero_by_symmetries = get_zero_by_symmetries(l_max);

      auto res = lm_vector<std::vector<std::complex<double>>>(l_max);
      for (auto l = 0u; l <= l_max; ++l) {
        for (auto m = -((int) l); m <= ((int) l); ++m) {
          if (zero_by_symmetries.at(l, m) == ZeroBySymmetries::True) {
            res.at(l, m) = std::vector<std::complex<double>>(mesh.size(), 0.0);
          } else {
            auto spherical_harmonics_values = get_spherical_harmonics_values(l, m);
            res.at(l, m) = integrate(l, m, spherical_harmonics_values);
          }
        }
      }
      return res;
    }

    const CutoffWSCell cutoff;
    const ShapeFunctionsConfig config;
    /// note that mesh is ordered in increasing order of r_WS.
    const std::vector<MeshPoint> mesh;

  private:
    lm_vector<ZeroBySymmetries> get_zero_by_symmetries(const unsigned int l_max) const
    {
      // initialize result to uniform False
      auto res = lm_vector<ZeroBySymmetries>(l_max);
      for (auto l = 0u; l <= l_max; ++l) {
        for (auto m = -((int) l); m <= ((int) l); ++m) {
          res.at(l, m) = ZeroBySymmetries::False;
        }
      }

      // get symmetries
      auto symmetries = SymmetryTransformationFactory::generate(cutoff.cell);

      // if we have inversion then we can use the parity of Spherical harmonics
      if (std::find(symmetries.begin(), symmetries.end(), Inversion()) != symmetries.end()) {
        for (auto l = 1u; l <= l_max; l += 2) {
          for (auto m = -((int) l); m <= ((int) l); ++m) {
            res.at(l, m) = ZeroBySymmetries::True;
          }
        }
      }


      // get rotations
      symmetries.erase(
        std::remove_if(
          symmetries.begin(), symmetries.end(),
          [](const Transformation &t)
          {
            return t.type != Transformation::Rotation;
          }),
        symmetries.end()
      );

      // handle rotations
      for (auto l = 0u; l <= l_max; ++l) {
        for (auto &rotation: symmetries) {
          auto euler_angles = reinterpret_cast<Rotation *>(&rotation)->get_euler_angles();
          auto wignerD = WignerDMatrix::calculate(euler_angles.alpha, euler_angles.beta, euler_angles.gamma, l);
          for (auto m = 0; m <= ((int) l); ++m) {
            auto isZero = ZeroBySymmetries::False;
            // we go through multiple symmetries, only have to do the math if it is still non-zero
            if (res.at(l, m) == ZeroBySymmetries::False) {
              for (auto i = -((int) l); i <= ((int) l); ++i) {
                if (i == m) { // we should have -1.0 in the diagonal
                  if (equalsWithTolerance(wignerD(l + m, l + i), -1.0 + 0.0i))
                    isZero = ZeroBySymmetries::True;
                  else {
                    isZero = ZeroBySymmetries::False;
                    break;
                  }
                } else if (!nearlyZero(std::abs(wignerD(l + m, l + i)))) { // and zeros everywhere else
                  isZero = ZeroBySymmetries::False;
                  break;
                }
              }
            }
            // we only update if it is false
            if (res.at(l, m) != ZeroBySymmetries::True) {
              res.at(l, m) = isZero;
              // since the integrand is real => we are symmetric in m.
              res.at(l, -static_cast<int>(m)) = isZero;
            }
          }
        }
      }
      return res;
    }

    std::vector<std::complex<double>> get_spherical_harmonics_values(const unsigned int l, const int m) const
    {
      auto res = std::vector<std::complex<double>>();
      res.reserve(mesh.size());
      std::transform(
        mesh.begin(), mesh.end(), std::back_inserter(res),
        [l, m](const MeshPoint &p)
        {
          return std::conj(Complex::spherical_harmonic(l, m, p.point)) * p.weight;
        });
      return res;
    }

    std::vector<std::complex<double>> integrate(
      const unsigned int l, const int m,
      const std::vector<std::complex<double>> &spherical_harmonics_values
    ) const
    {
      auto res = std::vector<std::complex<double>>();
      res.reserve(mesh.size());
      auto integral = 0.0i;

      auto itMesh = mesh.crbegin();
      const auto itMeshEnd = mesh.crend();
      auto itSH = spherical_harmonics_values.crbegin();
      const auto itSHEnd = spherical_harmonics_values.crend();

      auto pre_r_ws = itMesh->r_WS;
      for (; itMesh != itMeshEnd && itSH != itSHEnd; ++itMesh, ++itSH) {
        if (pre_r_ws - itMesh->r_WS > 1e-8) {
          res.push_back(integral * 4.0 * pi);
          pre_r_ws = itMesh->r_WS;
        }
        integral += *itSH;
      }
      res.push_back(integral * 4.0 * pi);
      res.shrink_to_fit();
      std::reverse(res.begin(), res.end());

      return res;
    }
  };
}

// TODO we should implement interpolation from going to r_WS mesh to custom input r mesh
ShapeFunctions::ShapeFunctions(const UnitCell3D &unit_cell_, unsigned int l_max, const ShapeFunctionsConfig &config_)
  : unit_cell(unit_cell_), shape_functions(ShapeFunctionsCalculator(unit_cell, config_).calculate(l_max)) {}
