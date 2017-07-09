#include <TestUtils/base.h>

#include <Math/IcosahedralQuadrature.h>
#include <fstream>
#include <Math/CommonFunctions.h>
#include <Math/SphericalHarmonics.h>

using namespace Math;

namespace
{
  void print_icosahedral_quadrature(unsigned int n)
  {
    auto res = IcosahedralQuadrature::generate(n);

    std::ofstream out_R;
    out_R.open("Icosahedral_" + std::to_string(n) + ".dat");
    for (auto r : res)
      out_R << std::setw(20) << std::setprecision(17) << std::fixed
            << r.first.x << "\t" << r.first.y << "\t" << r.first.z << "\t" << r.second << "\n";
    out_R.close();
  }
}

TEST(DISABLED_IcosahedralQuadrature, Print)
{
  static const auto n = 8u;
  for (auto i = 0u; i <= n; ++i) {
    print_icosahedral_quadrature(i);
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
        EXPECT_DOUBLE_EQ(p.first.length(), 1.0);
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
               + Complex::spherical_harmonic(0, 0, p.first)
                 * std::conj(Complex::spherical_harmonic(l, m, p.first))
                 * p.second;
      }
    );
    EXPECT_NEAR(std::abs(integral), 0.0, tol) << " l = " << l << " m = " << m;
  }
}

TEST(IcosahedralQuadrature, IntegrateSphericalHarmonics)
{
  static const auto n = 4u;
  for (auto i = 0u; i <= n; ++i) {
    auto mesh = IcosahedralQuadrature::generate(i);
    for (auto l = 1u; l <= 4u; ++l) {
      for (auto m = -((int) l); m <= ((int) l); ++m) {
        integrate_spherical_harmonics(l, m, 2e-16, mesh);
      }
    }
  }
}
