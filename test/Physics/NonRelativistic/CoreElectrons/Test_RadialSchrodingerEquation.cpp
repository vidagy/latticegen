#include <tuple>

#include <Physics/NonRelativistic/CoreElectrons/RadialSchrodingerEquation.h>
#include <TestUtils/Utils.h>
#include <TestUtils/base.h>

using namespace Physics::NonRelativistic::CoreElectrons;

namespace
{
  struct CoulombReferenceSolutions
  {
    CoulombReferenceSolutions(const std::shared_ptr<const ExponentialMesh> &mesh, double Q_)
      : z(std::vector<double>(mesh->points.size(), Q_), mesh), Q(Q_)
    {
      reference_R10 = this->generate([this](double r) { return R10(r); });
      reference_dR_dr_10 = this->generate([this](double r) { return dR_dr_10(r); });
      reference_R20 = this->generate([this](double r) { return R20(r); });
      reference_dR_dr_20 = this->generate([this](double r) { return dR_dr_20(r); });
      reference_R21 = this->generate([this](double r) { return R21(r); });
      reference_dR_dr_21 = this->generate([this](double r) { return dR_dr_21(r); });
      reference_R30 = this->generate([this](double r) { return R30(r); });
      reference_dR_dr_30 = this->generate([this](double r) { return dR_dr_30(r); });
      reference_R31 = this->generate([this](double r) { return R31(r); });
      reference_dR_dr_31 = this->generate([this](double r) { return dR_dr_31(r); });
      reference_R32 = this->generate([this](double r) { return R32(r); });
      reference_dR_dr_32 = this->generate([this](double r) { return dR_dr_32(r); });
    }

    double energy(int n) const { return -Q * Q / 2.0 / n / n; }

    double R10(double r) const { return 2.0 * pow(Q, 3.0 / 2.0) * r * exp(-Q * r); }

    double dR_dr_10(double r) const { return 2.0 * pow(Q, 3.0 / 2.0) * (exp(-Q * r) - Q * r * exp(-Q * r)); }

    double R20(double r) const
    {
      return 1.0 / sqrt(2.0) * pow(Q, 3.0 / 2.0) * r * exp(-Q * r / 2.0) * (1.0 - Q * r / 2.0);
    }

    double dR_dr_20(double r) const
    {
      return 1.0 / sqrt(2.0) * pow(Q, 3.0 / 2.0) *
        (
          -Q / 2.0 * exp(-Q * r / 2.0) * (r - Q * r * r / 2.0)
          + exp(-Q * r / 2.0) * (1.0 - Q * r)
        );
    }

    double R21(double r) const { return 1.0 / (2.0 * sqrt(6.0)) * pow(Q, 5.0 / 2.0) * r * r * exp(-Q * r / 2.0); }

    double dR_dr_21(double r) const
    {
      return 1.0 / (2.0 * sqrt(6.0)) * pow(Q, 5.0 / 2.0) *
        (2.0 * r * exp(-Q * r / 2.0) - Q / 2.0 * r * r * exp(-Q * r / 2.0));
    }

    double R30(double r) const
    {
      return 2.0 / 3.0 / sqrt(3.0) * pow(Q, 3.0 / 2.0) * r * exp(-Q * r / 3.0) *
        (1.0 - Q * r * 2.0 / 3.0 + 2.0 / 27.0 * Q * Q * r * r);
    }

    double dR_dr_30(double r) const
    {
      return 2.0 / 3.0 / sqrt(3.0) * pow(Q, 3.0 / 2.0) *
        (
          -Q / 3.0 * exp(-Q * r / 3.0) *
          (r - Q * r * r * 2.0 / 3.0 + 2.0 / 27.0 * Q * Q * r * r * r)
          + exp(-Q * r / 3.0) *
            (1.0 - Q * r * 4.0 / 3.0 + 2.0 / 9.0 * Q * Q * r * r)
        );
    }

    double R31(double r) const
    {
      return 8.0 / 27.0 / sqrt(6.0) * pow(Q, 5.0 / 2.0) * r * r * exp(-Q * r / 3.0) * (1.0 - Q * r / 6.0);
    }

    double dR_dr_31(double r) const
    {
      return 8.0 / 27.0 / sqrt(6.0) * pow(Q, 5.0 / 2.0) *
        (
          -Q / 3.0 * exp(-Q * r / 3.0) * (r * r - Q * r * r * r / 6.0)
          + exp(-Q * r / 3.0) * (2.0 * r - Q * r * r / 2.0)
        );
    }

    double R32(double r) const { return 4.0 / 81.0 / sqrt(30.0) * pow(Q, 7.0 / 2.0) * r * r * r * exp(-Q * r / 3.0); }

    double dR_dr_32(double r) const
    {
      return 4.0 / 81.0 / sqrt(30.0) * pow(Q, 7.0 / 2.0) *
        (3.0 * r * r * exp(-Q * r / 3.0) - Q / 3.0 * r * r * r * exp(-Q * r / 3.0));
    }

    std::vector<double> generate(std::function<double(double)> f) const
    {
      std::vector<double> res;
      res.reserve(z.r->points.size());
      for (auto r: z.r->points) {
        res.push_back(f(r));
      }
      return res;
    }

    std::vector<double> get_R(int n, int l) const
    {
      switch (n) {
        case 1:
          switch (l) {
            case 0:
              return reference_R10;
            default:
              throw std::invalid_argument("not implemented");
          }
        case 2:
          switch (l) {
            case 0:
              return reference_R20;
            case 1:
              return reference_R21;
            default:
              throw std::invalid_argument("not implemented");
          }
        case 3:
          switch (l) {
            case 0:
              return reference_R30;
            case 1:
              return reference_R31;
            case 2:
              return reference_R32;
            default:
              throw std::invalid_argument("not implemented");
          }
        default:
          throw std::invalid_argument("not implemented");
      }
    }

    std::vector<double> get_dR_dr(int n, int l) const
    {
      switch (n) {
        case 1:
          switch (l) {
            case 0:
              return reference_dR_dr_10;
            default:
              throw std::invalid_argument("not implemented");
          }
        case 2:
          switch (l) {
            case 0:
              return reference_dR_dr_20;
            case 1:
              return reference_dR_dr_21;
            default:
              throw std::invalid_argument("not implemented");
          }
        case 3:
          switch (l) {
            case 0:
              return reference_dR_dr_30;
            case 1:
              return reference_dR_dr_31;
            case 2:
              return reference_dR_dr_32;
            default:
              throw std::invalid_argument("not implemented");
          }
        default:
          throw std::invalid_argument("not implemented");
      }
    }

    const EffectiveCharge z;
    const double Q;
    std::vector<double> reference_R10;
    std::vector<double> reference_dR_dr_10;
    std::vector<double> reference_R20;
    std::vector<double> reference_dR_dr_20;
    std::vector<double> reference_R21;
    std::vector<double> reference_dR_dr_21;
    std::vector<double> reference_R30;
    std::vector<double> reference_dR_dr_30;
    std::vector<double> reference_R31;
    std::vector<double> reference_dR_dr_31;
    std::vector<double> reference_R32;
    std::vector<double> reference_dR_dr_32;
  };
}

namespace Physics
{
  namespace NonRelativistic
  {
    namespace CoreElectrons
    {
      class TestAccessor
      {
      public:
        static void adams_moulton_method(
          const AdamsIntegrator &adamsIntegrator, std::vector<double> &R, std::vector<double> &dR_dr, int from, int to
        )
        {
          return adamsIntegrator.adams_moulton_method(R, dR_dr, from, to);
        }

        static std::vector<double> get_adams_parameters(int quadrature)
        {
          return AdamsIntegrator::get_adams_parameters(quadrature);
        }
      };
    }
  }
}

TEST(TestAdamsIntegrator, AdamsMoultonMethod)
{
  auto mesh = std::make_shared<const ExponentialMesh>(0.01, 50, 100, 1.0);
  auto reference = CoulombReferenceSolutions(mesh, 1.0);

  auto R = reference.reference_R10;
  auto dR_dR = reference.reference_dR_dr_10;
  auto z = std::vector<double>(mesh->points.size(), 1.0);

  auto ai = AdamsIntegrator(mesh, z, reference.energy(1), 0, 97, 77, AdamsIntegratorConfig());

  auto from = 20;
  auto to = 70;
  Physics::NonRelativistic::CoreElectrons::TestAccessor::adams_moulton_method(ai, R, dR_dR, from, to);

  for (auto i = 0; i < from; ++i) {
    EXPECT_EQ(R[i], reference.reference_R10[i]);
    EXPECT_EQ(dR_dR[i], reference.reference_dR_dr_10[i]);
  }
  for (auto i = from; i <= to; ++i) {
    EXPECT_NEAR(R[i], reference.reference_R10[i], 1.6e-8);
    EXPECT_NEAR(dR_dR[i], reference.reference_dR_dr_10[i], 1e-8);
  }
  for (auto i = to + 1u; i < R.size(); ++i) {
    EXPECT_EQ(R[i], reference.reference_R10[i]);
    EXPECT_EQ(dR_dR[i], reference.reference_dR_dr_10[i]);
  }
}

namespace
{
  static const std::vector<double> adams_params_1{0.5, 0.5};
  static const std::vector<double> adams_params_2{-1.0 / 12.0, 8.0 / 12.0, 5.0 / 12.0};
  static const std::vector<double> adams_params_3{1.0 / 24.0, -5.0 / 24.0, 19.0 / 24.0, 9.0 / 24.0};
  static const std::vector<double> adams_params_4{-19.0 / 720.0, 106.0 / 720.0, -264.0 / 720.0, 646.0 / 720.0,
                                                  251.0 / 720.0};
  static const std::vector<double> adams_params_5{27.0 / 1440.0, -173.0 / 1440.0, 482.0 / 1440.0, -798.0 / 1440.0,
                                                  1427.0 / 1440.0, 475.0 / 1440.0};
  static const std::vector<double> adams_params_6{-863.0 / 60480.0, 6312 / 60480.0, -20211.0 / 60480.0,
                                                  37504.0 / 60480.0, -46461.0 / 60480.0, 65112.0 / 60480.0,
                                                  19087.0 / 60480.0};
  static const std::vector<double> adams_params_7{1375.0 / 120960.0, -11351.0 / 120960.0, 41499.0 / 120960.0,
                                                  -88547.0 / 120960.0,
                                                  123133.0 / 120960.0, -121797.0 / 120960.0, 139849.0 / 120960.0,
                                                  36799.0 / 120960.0};
  static const std::vector<double> adams_params_8{-33953.0 / 3628800.0, 312874.0 / 3628800.0, -1291214.0 / 3628800.0,
                                                  3146338.0 / 3628800.0, -5033120.0 / 3628800.0, 5595358.0 / 3628800.0,
                                                  -4604594.0 / 3628800.0, 4467094.0 / 3628800.0, 1070017.0 / 3628800.0};
}

TEST(TestAdamsIntegrator, AdamsParameters)
{
  auto eps = std::numeric_limits<double>::epsilon();
  EXPECT_THAT(TestAccessor::get_adams_parameters(1), ::testing::Pointwise(NearWithTolerance(eps), adams_params_1));
  EXPECT_THAT(TestAccessor::get_adams_parameters(2), ::testing::Pointwise(NearWithTolerance(eps), adams_params_2));
  EXPECT_THAT(TestAccessor::get_adams_parameters(3), ::testing::Pointwise(NearWithTolerance(eps), adams_params_3));
  EXPECT_THAT(TestAccessor::get_adams_parameters(4), ::testing::Pointwise(NearWithTolerance(eps), adams_params_4));
  EXPECT_THAT(TestAccessor::get_adams_parameters(5), ::testing::Pointwise(NearWithTolerance(eps), adams_params_5));
  EXPECT_THAT(TestAccessor::get_adams_parameters(6), ::testing::Pointwise(NearWithTolerance(eps), adams_params_6));
  EXPECT_THAT(TestAccessor::get_adams_parameters(7), ::testing::Pointwise(NearWithTolerance(eps), adams_params_7));
  EXPECT_THAT(TestAccessor::get_adams_parameters(8),
              ::testing::Pointwise(NearWithTolerance(2.0 * eps), adams_params_8));
}

namespace
{
  void compare_to_reference(
    double Z, const std::shared_ptr<const ExponentialMesh> &mesh,
    int n, int l, double tol = std::numeric_limits<double>::epsilon(), bool print_logs = false
  )
  {
    auto reference = CoulombReferenceSolutions(mesh, Z);
    auto sch = RadialSchrodingerEquation(EffectiveCharge(std::vector<double>(mesh->points.size(), Z), mesh));
    auto solution = sch.solve(n, l, reference.energy(n));

    if (print_logs) {
      Utils::log(solution.R,
                 "ScaledSolution_R" + std::to_string(n) + std::to_string(l) + "_" + std::to_string(mesh->scale));
      Utils::log(solution.dR_dr,
                 "ScaledSolution_dR_dr" + std::to_string(n) + std::to_string(l) + "_" + std::to_string(mesh->scale));
    }

    auto reference_R = reference.get_R(n, l);
    auto reference_dR_dr = reference.get_dR_dr(n, l);

    if (print_logs) {
      Utils::log(reference_R, "ScaledReferenceSolution_reference_R" + std::to_string(n) + std::to_string(l) + "_" +
                              std::to_string(mesh->scale));
      Utils::log(reference_dR_dr,
                 "ScaledReferenceSolution_reference_dR_dr" + std::to_string(n) + std::to_string(l) + "_" +
                 std::to_string(mesh->scale));
    }

    EXPECT_LE(solution.number_of_iteration, 2) << "for n = " << n << " l = " << l;
    EXPECT_NEAR(solution.E, reference.energy(n), tol) << "for n = " << n << " l = " << l;

    EXPECT_THAT(solution.R, ::testing::Pointwise(NearWithTolerance(tol), reference_R))
            << "for n = " << n << " l = " << l;
    EXPECT_THAT(solution.dR_dr, ::testing::Pointwise(NearWithTolerance(tol), reference_dR_dr))
            << "for n = " << n << " l = " << l;

    if (print_logs) {
      double max_R = 0.0;
      double max_dR = 0.0;
      for (auto i = 0u; i < reference_R.size(); ++i) {
        if (fabs(solution.R[i] - reference_R[i]) > max_R)
          max_R = fabs(solution.R[i] - reference_R[i]);
        if (fabs(solution.R[i] - reference_R[i]) > max_dR)
          max_dR = fabs(solution.dR_dr[i] - reference_dR_dr[i]);
      }
      std::cout << "for n = " << n << " l = " << l << " scale " << mesh->scale << " max_R = " << std::setprecision(16)
                << std::scientific
                << std::setw(21) << max_R << " msx_dr " << max_dR << std::endl;
    }
  }
}

TEST(TestRadialSchrodingerEquation, CompareToCoulombReference)
{
  auto mesh = std::make_shared<const ExponentialMesh>(0.00001, 200, 200, 1.0);
  double Z = 1.0;
  compare_to_reference(Z, mesh, 1, 0, 3e-8);
  compare_to_reference(Z, mesh, 2, 0, 2e-7);
  compare_to_reference(Z, mesh, 2, 1, 1e-7);
  compare_to_reference(Z, mesh, 3, 0, 5e-7);
  compare_to_reference(Z, mesh, 3, 1, 5e-7);
  compare_to_reference(Z, mesh, 3, 2, 5e-7);
}

TEST(TestRadialSchrodingerEquation, CompareToCoulombReferenceScaledMesh)
{
  auto mesh = std::make_shared<const ExponentialMesh>(0.00001, 200, 200, 0.7);
  double Z = 1.0;
  Utils::log(mesh->points, "NotScaledSolution_r");
  compare_to_reference(Z, mesh, 1, 0, 9e-10);
  compare_to_reference(Z, mesh, 2, 0, 7e-9);
  compare_to_reference(Z, mesh, 2, 1, 4e-9);
  compare_to_reference(Z, mesh, 3, 0, 4e-8);
  compare_to_reference(Z, mesh, 3, 1, 3e-8);
  compare_to_reference(Z, mesh, 3, 2, 9e-9);
}

TEST(DISABLED_TestRadialSchrodingerEquation, GenerateInputForParametrizedExponentialMesh)
{
  for (auto scale = 0.1; scale < 1.05; scale += 0.1) {
    auto mesh = std::make_shared<const ExponentialMesh>(0.00001, 200, 200, scale);
    Utils::log(mesh->points, "ScaledSolution_r_" + std::to_string(scale));
    double Z = 1.0;
    compare_to_reference(Z, mesh, 1, 0, 3e-1, true);
    compare_to_reference(Z, mesh, 2, 0, 2e-1, true);
    compare_to_reference(Z, mesh, 2, 1, 1e-1, true);
    compare_to_reference(Z, mesh, 3, 0, 5e-1, true);
    compare_to_reference(Z, mesh, 3, 1, 5e-1, true);
    compare_to_reference(Z, mesh, 3, 2, 5e-1, true);
  }
}
