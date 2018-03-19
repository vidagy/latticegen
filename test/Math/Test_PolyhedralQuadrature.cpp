#include <TestUtils/base.h>

#include <Math/PolyhedralQuadrature.h>
#include <fstream>
#include <numeric>
#include <Math/CommonFunctions.h>
#include <Math/SphericalHarmonics.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <Math/LebedevQuadrature.h>

using namespace Math;

namespace
{
  void print_icosahedral_quadrature(const std::string &name, const Quadrature &res)
  {
    std::ofstream out_R;
    out_R.open(name + ".dat");
    for (auto r : res)
      out_R << std::setw(20) << std::setprecision(17) << std::fixed
            << r.first(0) << "\t" << r.first(1) << "\t" << r.first(2) << "\t" << r.second << "\n";
    out_R.close();
  }
}

TEST(DISABLED_IcosahedralQuadrature, Print)
{
  static const auto n = 8u;
  for (auto i = 0u; i <= n; ++i) {
    auto res = IcosahedralQuadrature::generate(i);
    print_icosahedral_quadrature("Icosahedral_" + std::to_string(i), res);
  }
}

TEST(IcosahedralQuadrature, NumberOfPoints)
{
  static const auto n = 4u;
  for (auto i = 0u; i <= n; ++i) {
    auto res = IcosahedralQuadrature::generate(i);
    EXPECT_EQ(res.size(), (unsigned int) (20 * Math::pow(4, i)));
  }
}

TEST(IcosahedralQuadrature, LengthOfPoints)
{
  static const auto n = 4u;
  for (auto i = 0u; i <= n; ++i) {
    auto res = IcosahedralQuadrature::generate(i);
    std::for_each(
      res.begin(), res.end(),
      [](const std::pair<Point3D, double> &p)
      {
        EXPECT_DOUBLE_EQ(p.first.norm(), 1.0);
      }
    );
  }
}

TEST(IcosahedralQuadrature, SumWeights)
{
  static const auto n = 4u;
  for (auto i = 0u; i <= n; ++i) {
    auto res = IcosahedralQuadrature::generate(i);
    auto sum_weights = std::accumulate(
      res.begin(), res.end(), 0.0,
      [](double sum, const std::pair<Point3D, double> &p)
      {
        return sum + p.second;
      }
    );
    EXPECT_NEAR(sum_weights, 1.0, 3e-14);
  }
}

namespace
{
  void
  integrate_spherical_harmonics(unsigned int l, int m, double tol, const std::vector<std::pair<Point3D, double>> &mesh)
  {
    using namespace std::complex_literals;
    auto integral = std::accumulate(
      mesh.begin(), mesh.end(), 0.0i,
      [l, m](std::complex<double> sum, const std::pair<Point3D, double> &p)
      {
        return sum
               + 0.5 / sqrt(pi)
                 * std::conj(Complex::spherical_harmonic(l, m, p.first))
                 * p.second;
      }
    );
    EXPECT_NEAR(std::abs(integral), 0.0, tol) << " mesh.size = " << mesh.size() << " l = " << l << " m = " << m;
  }
}

TEST(IcosahedralQuadrature, IntegrateSphericalHarmonics)
{
  static const auto n = 5u;
  for (auto i = 0u; i <= n; ++i) {
    auto mesh = IcosahedralQuadrature::generate(i);
    // up to l = 5, icosahedral mesh has the symmetry of spherical harmonics, so there is no large diff.
    // However for larger l, the errors are quite large and they decrease slowly with n
    for (auto l = 1u; l <= 5u; ++l) {
      for (auto m = -((int) l); m <= ((int) l); ++m) {
        integrate_spherical_harmonics(l, m, 2e-16, mesh);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TEST(DISABLED_OctahedralQuadrature, Print)
{
  static const auto n = 6u;
  for (auto i = 0u; i <= n; ++i) {
    auto res = OctahedralQuadrature::generate(i);
    print_icosahedral_quadrature("Octahedral_" + std::to_string(i), res);
  }
}

TEST(OctahedralQuadrature, NumberOfPoints)
{
  static const auto n = 4u;
  for (auto i = 0u; i <= n; ++i) {
    auto res = OctahedralQuadrature::generate(i);
    EXPECT_EQ(res.size(), (unsigned int) (6 * Math::pow(4, i)));
  }
}

TEST(OctahedralQuadrature, LengthOfPoints)
{
  static const auto n = 4u;
  for (auto i = 0u; i <= n; ++i) {
    auto res = OctahedralQuadrature::generate(i);
    std::for_each(
      res.begin(), res.end(),
      [](const std::pair<Point3D, double> &p)
      {
        EXPECT_DOUBLE_EQ(p.first.norm(), 1.0);
      }
    );
  }
}

TEST(OctahedralQuadrature, SumWeights)
{
  static const auto n = 4u;
  for (auto i = 0u; i <= n; ++i) {
    auto res = OctahedralQuadrature::generate(i);
    auto sum_weights = std::accumulate(
      res.begin(), res.end(), 0.0,
      [](double sum, const std::pair<Point3D, double> &p)
      {
        return sum + p.second;
      }
    );
    EXPECT_NEAR(sum_weights, 1.0, 3e-14);
  }
}

TEST(OctahedralQuadrature, IntegrateSphericalHarmonics)
{
  static const auto n = 5u;
  for (auto i = 0u; i <= n; ++i) {
    auto mesh = OctahedralQuadrature::generate(i);
    // up to l = 3, icosahedral mesh has the symmetry of spherical harmonics, so there is no large diff.
    // However for larger l, the errors are quite large and they decrease slowly with n
    for (auto l = 1u; l <= 3u; ++l) {
      for (auto m = -((int) l); m <= ((int) l); ++m) {
        integrate_spherical_harmonics(l, m, 2e-16, mesh);
      }
    }
  }
}
