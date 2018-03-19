#include <TestUtils/base.h>

#include <Math/WignerDMatrix.h>
#include <Math/CommonFunctions.h>
#include <Math/PolyhedralQuadrature.h>
#include <Geometry/Transformations.h>
#include <Math/SphericalHarmonics.h>

using namespace Math;
using namespace Geometry;

namespace
{
  inline int c(int l, int m)
  {
    return l + m;
  }

  void check_symmetry(const Eigen::MatrixXd &d)
  {
    auto size = d.size();
    int l_int = (int) (std::lround(sqrt(size)) - 1) / 2;
    for (auto m_p = -l_int; m_p <= l_int; ++m_p) {
      for (auto m = -l_int; m <= m_p; ++m) {
        const auto c_m = c(l_int, m);
        const auto c_m_p = c(l_int, m_p);
        EXPECT_DOUBLE_EQ(d(c_m_p, c_m), sign(m - m_p) * d(c_m, c_m_p));

        const auto c_m_m = c(l_int, -m);
        const auto c_m_m_p = c(l_int, -m_p);
        EXPECT_DOUBLE_EQ(d(c_m_p, c_m), d(c_m_m, c_m_m_p));
      }
    }
  }

  const static auto tol = 1e-15;
}

TEST(WignerDMatrix, smallDSymmetry)
{
  const auto n = 10;
  for (auto i = 0; i < n; ++i) {
    const auto beta = (pi * i) / n;
    Eigen::MatrixXd d1 = WignerDMatrix::calculate_small(beta, 1);
    check_symmetry(d1);
    Eigen::MatrixXd d2 = WignerDMatrix::calculate_small(beta, 2);
    check_symmetry(d2);
    Eigen::MatrixXd d3 = WignerDMatrix::calculate_small(beta, 3);
    check_symmetry(d3);
    Eigen::MatrixXd d4 = WignerDMatrix::calculate_small(beta, 4);
    check_symmetry(d4);
  }
}

TEST(WignerDMatrix, smallD1ToReference)
{
  const auto n = 10;
  const static auto l = 1;
  for (auto i = 0; i < n; ++i) {
    const auto beta = (pi * i) / n;
    Eigen::MatrixXd d = WignerDMatrix::calculate_small(beta, l);

    EXPECT_NEAR(d(c(l, 1), c(l, 1)), (1.0 + cos(beta)) / 2.0, tol);
    EXPECT_NEAR(d(c(l, 1), c(l, 0)), -sin(beta) / sqrt(2.0), tol);
    EXPECT_NEAR(d(c(l, 1), c(l, -1)), (1.0 - cos(beta)) / 2.0, tol);
    EXPECT_NEAR(d(c(l, 0), c(l, 0)), cos(beta), tol);
  }
}

TEST(WignerDMatrix, smallD2ToReference)
{
  const auto n = 10;
  const static auto l = 2;
  for (auto i = 0; i < n; ++i) {
    const auto beta = (pi * i) / n;
    Eigen::MatrixXd d = WignerDMatrix::calculate_small(beta, l);

    EXPECT_NEAR(d(c(l, 2), c(l, 2)), Math::pow(1.0 + cos(beta), 2) / 4.0, tol);
    EXPECT_NEAR(d(c(l, 2), c(l, 1)), -sin(beta) * (1.0 + cos(beta)) / 2.0, tol);
    EXPECT_NEAR(d(c(l, 2), c(l, 0)), sqrt(3.0 / 8.0) * sin(beta) * sin(beta), tol);
    EXPECT_NEAR(d(c(l, 2), c(l, -1)), -sin(beta) * (1.0 - cos(beta)) / 2.0, tol);
    EXPECT_NEAR(d(c(l, 2), c(l, -2)), Math::pow(1.0 - cos(beta), 2) / 4.0, tol);
    EXPECT_NEAR(d(c(l, 1), c(l, 1)), (2 * cos(beta) * cos(beta) + cos(beta) - 1) / 2.0, tol);
    EXPECT_NEAR(d(c(l, 1), c(l, 0)), -sqrt(3.0 / 8.0) * sin(2 * beta), tol);
    EXPECT_NEAR(d(c(l, 1), c(l, -1)), (-2 * cos(beta) * cos(beta) + cos(beta) + 1) / 2.0, tol);
    EXPECT_NEAR(d(c(l, 0), c(l, 0)), (3 * cos(beta) * cos(beta) - 1) / 2.0, tol);
  }
}

namespace
{
  void check_spherical_harmonics_rotation(
    const Point3D &p, const Point3D &rotated_p, int l, const Eigen::MatrixXcd &d, const double tolerance = tol)
  {
    Eigen::VectorXcd Y_lm_at_rotated(2 * l + 1);
    for (auto m = -l; m <= l; ++m) {
      Y_lm_at_rotated(c(l, m)) = Complex::spherical_harmonic((unsigned) l, m, rotated_p);
    }

    Eigen::VectorXcd transformed_Y = d * Y_lm_at_rotated;

    for (auto m_p = -l; m_p <= l; ++m_p) {
      auto sh = Complex::spherical_harmonic((unsigned) l, m_p, p);
      EXPECT_NEAR(sh.real(), transformed_Y(c(l, m_p)).real(), tolerance);
      EXPECT_NEAR(sh.imag(), transformed_Y(c(l, m_p)).imag(), tolerance);
    }
  }
}

TEST(WignerDMatrix, smallDRotateSphericalHarmonics)
{
  auto mesh = IcosahedralQuadrature::generate(0);
  const auto n = 10;
  for (auto i = 1; i < n; ++i) {
    const auto beta = (pi * i) / n;
    Eigen::MatrixXcd d1 = WignerDMatrix::calculate_small(beta, 1).cast<std::complex<double>>();
    Eigen::MatrixXcd d2 = WignerDMatrix::calculate_small(beta, 2).cast<std::complex<double>>();
    Eigen::MatrixXcd d3 = WignerDMatrix::calculate_small(beta, 3).cast<std::complex<double>>();

    const auto rotation = Geometry::Rotation(Point3D(0, beta, 0));

    for (const auto &mp: mesh) {
      const auto p = mp.first;
      check_spherical_harmonics_rotation(p, rotation(p), 1, d1);
      check_spherical_harmonics_rotation(p, rotation(p), 2, d2);
      check_spherical_harmonics_rotation(p, rotation(p), 3, d3);
    }
  }
}

TEST(WignerDMatrix, fullDRotateSphericalHarmonics)
{
  auto mesh = IcosahedralQuadrature::generate(0);
  const auto n = 5;
  for (auto i = 1; i <= n; ++i) {
    for (auto j = 1; j <= n; ++j) {
      for (auto k = 1; k <= n; ++k) {
        const auto alpha = (2 * pi * i) / n;
        const auto beta = (2 * pi * j) / n;
        const auto gamma = (2 * pi * k) / n;
        Eigen::MatrixXcd d1 = WignerDMatrix::calculate(alpha, beta, gamma, 1);
        Eigen::MatrixXcd d2 = WignerDMatrix::calculate(alpha, beta, gamma, 2);
        Eigen::MatrixXcd d3 = WignerDMatrix::calculate(alpha, beta, gamma, 3);

        const auto rotation = Rotation(alpha, beta, gamma);

        for (const auto &mp: mesh) {
          const auto p = mp.first;
          const auto rotated_p = rotation(p);
          check_spherical_harmonics_rotation(p, rotated_p, 1, d1, 1e-14);
          check_spherical_harmonics_rotation(p, rotated_p, 2, d2, 1e-14);
          check_spherical_harmonics_rotation(p, rotated_p, 3, d3, 1e-14);
        }
      }
    }
  }
}
