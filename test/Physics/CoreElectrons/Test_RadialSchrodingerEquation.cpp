#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Physics/CoreElectrons/RadialSchrodingerEquation.h>

using namespace Physics::CoreElectrons;

namespace
{
  std::vector<double> V(const double Z, const std::shared_ptr<const ExponentialMesh> &mesh)
  {
    std::vector<double> res;
    res.reserve(mesh->points.size());
    for (auto x: mesh->points) {
      res.push_back(-Z / x);
    }
    return res;
  }
}

TEST(TestRadialSchrodingerEquation, ctor)
{
  auto mesh = std::make_shared<const ExponentialMesh>(0.00001, 50, 500);
  auto potential = V(1.0, mesh);

  auto sch = RadialSchrodingerEquation(EffectiveCharge(potential, mesh));
}
