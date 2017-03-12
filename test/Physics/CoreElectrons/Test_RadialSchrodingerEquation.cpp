#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Physics/CoreElectrons/RadialSchrodingerEquation.h>
#include <Physics/CoreElectrons/AdamsIntegrator.h>
#include "Utils.h"

using namespace Physics::CoreElectrons;

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
      reference_R21 = this->generate([this](double r) { return R21(r); });
      reference_R30 = this->generate([this](double r) { return R30(r); });
      reference_R31 = this->generate([this](double r) { return R31(r); });
      reference_R32 = this->generate([this](double r) { return R32(r); });
    }

    double energy(int n) const { return -Q * Q / 2.0 / n / n; }

    double R10(double r) const { return 2.0 * pow(Q, 3.0 / 2.0) * r * exp(-Q * r); }

    double dR_dr_10(double r) const { return 2.0 * pow(Q, 3.0 / 2.0) * (exp(-Q * r) - Q * r * exp(-Q * r)); }

    double R20(double r) const
    {
      return 1.0 / sqrt(2.0) * pow(Q, 3.0 / 2.0) * r * exp(-Q * r / 2.0) * (1.0 - Q * r / 2.0);
    }

    double R21(double r) const { return 1.0 / (2.0 * sqrt(6.0)) * pow(Q, 5.0 / 2.0) * r * r * exp(-Q * r / 2.0); }

    double R30(double r) const
    {
      return 2.0 / 3.0 / sqrt(3.0) * pow(Q, 3.0 / 2.0) * r * exp(-Q * r / 3.0) *
             (1.0 - Q * r * 2.0 / 3.0 + 2.0 / 27.0 * Q * Q * r * r);
    }

    double R31(double r) const
    {
      return 8.0 / 27.0 / sqrt(6.0) * pow(Q, 5.0 / 2.0) * r * r * exp(-Q * r / 3.0) * (1.0 - Q * r / 6.0);
    }

    double R32(double r) const { return 4.0 / 81.0 / sqrt(30.0) * pow(Q, 7.0 / 2.0) * r * r * r * exp(-Q * r / 3.0); }

    std::vector<double> generate(std::function<double(double)> f) const
    {
      std::vector<double> res;
      res.reserve(z.r->points.size());
      for (auto r: z.r->points) {
        res.push_back(f(r));
      }
      return res;
    }

    const EffectiveCharge z;
    const double Q;
    std::vector<double> reference_R10;
    std::vector<double> reference_dR_dr_10;
    std::vector<double> reference_R20;
    std::vector<double> reference_R21;
    std::vector<double> reference_R30;
    std::vector<double> reference_R31;
    std::vector<double> reference_R32;
  };
}

TEST(DISABLED_TestRadialSchrodingerEquation, ReferenceSolution)
{
  auto mesh = std::make_shared<const ExponentialMesh>(0.00001, 50, 100);
  auto reference = CoulombReferenceSolutions(mesh, 1.0);

  Utils::log(mesh->points, "ReferenceSolution_r");
  Utils::log(reference.reference_R10, "ReferenceSolution_reference_R10");
  Utils::log(reference.reference_dR_dr_10, "ReferenceSolution_reference_dR_dr_10");
  Utils::log(reference.reference_R20, "ReferenceSolution_reference_R20");
  Utils::log(reference.reference_R21, "ReferenceSolution_reference_R21");
  Utils::log(reference.reference_R30, "ReferenceSolution_reference_R30");
  Utils::log(reference.reference_R31, "ReferenceSolution_reference_R31");
  Utils::log(reference.reference_R32, "ReferenceSolution_reference_R32");
}

namespace Physics
{
  namespace CoreElectrons
  {
    class TestAccessor
    {
    public:
      static void adams_moulton_method(
        const AdamsIntegrator &adamsIntegrator, std::vector<double> &R, std::vector<double> &dR_dr, unsigned long from,
        unsigned long to)
      {
        return adamsIntegrator.adams_moulton_method(R, dR_dr, from, to);
      }
    };
  }
}

TEST(TestAdamsIntegrator, AdamsMoultonMethod)
{
  auto mesh = std::make_shared<const ExponentialMesh>(0.01, 50, 100);
  auto reference = CoulombReferenceSolutions(mesh, 1.0);

  auto R = reference.reference_R10;
  auto dR_dR = reference.reference_dR_dr_10;
  auto z = std::vector<double>(mesh->points.size(), 1.0);

  auto ai = AdamsIntegrator(mesh, z, reference.energy(1), 0, 97, 77, 8, 0);

  auto from = 20;
  auto to = 70;
  Physics::CoreElectrons::TestAccessor::adams_moulton_method(ai, R, dR_dR, from, to);

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

TEST(DISABLED_TestRadialSchrodingerEquation, ctor)
{
  auto mesh = std::make_shared<const ExponentialMesh>(0.00001, 50, 100);
  auto Z = 1.0;
  auto n = 1;
  auto energy = -Z * Z / 2.0 / n / n;
  auto l = 0;

  auto sch = RadialSchrodingerEquation(EffectiveCharge(std::vector<double>(mesh->points.size(), Z), mesh));
  auto res = sch.solve(n, l, energy, 8);
}


