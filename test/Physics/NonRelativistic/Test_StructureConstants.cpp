#include <tuple>

#include <Physics/NonRelativistic/StructureConstants.h>
#include <TestUtils/base.h>
#include <Geometry/IrreducibleWedge.h>

using namespace Physics::NonRelativistic;
using namespace std::complex_literals;

namespace
{
}

TEST(TestStructureConstants, Real)
{
  auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto structure_constant = StructureConstants(unit_cell);
  auto nothing = 0i;
  for (auto rez = 0.01; rez < 1.0; rez += 0.1) {
    for (auto n = 0; n < 4; ++n) {
      for (auto m = 0; m < 4; ++m) {
        for (auto l = 0; l < 4; ++l) {
          if (!(m == 0 && n == 0 && l == 0)) {
            auto res = structure_constant.calculate_real_space(
              0, 0, 0, 0, Coordinates3D(0, 0, 0), Coordinates3D(n, m, l), rez + 0.1i
            );
            // std::cout << unit_cell.at(Coordinates3D(0, 0, 0) - Coordinates3D(n, m, l)).length()
            //           << "  " << res << std::endl;
            nothing += res;
          }
        }
      }
    }
  }
}


namespace Physics
{
  namespace NonRelativistic
  {
    class TestAccessor
    {
    public:
      static std::complex<double> D1(
        const StructureConstants &structure_constants,
        unsigned int l, int m, const Vector3D &k, const std::complex<double> &z
      )
      {
        return structure_constants.D1(l, m, k, z);
      }

      static std::complex<double> D2(
        const StructureConstants &structure_constants,
        unsigned int l, int m, const Vector3D &k, const std::complex<double> &z
      )
      {
        return structure_constants.D2(l, m, k, z);
      }

      static std::complex<double> D3(
        const StructureConstants &structure_constants,
        unsigned int l, int m, const Vector3D &k, const std::complex<double> &z
      )
      {
        return structure_constants.D3(l, m, k, z);
      }
    };
  }
}

TEST(TestStructureConstants, Reciprocal)
{
  auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto reciprocal_unit_cell = ReciprocalUnitCell3D(unit_cell);
  auto irreduc_wedge = IrreducibleWedge::get_irreducible_wedge(reciprocal_unit_cell, 5);

  // auto k = irreduc_wedge[irreduc_wedge.size()/2];

  // auto res = structure_constant.calculate_reciprocal_space(0, 0, 0, 0, k, 0.1 + 0.1i);

  for (auto k: irreduc_wedge) {
    auto params = StructureConstantsConfig(
      10.0,
      3.0
    );
    auto structure_constant = StructureConstants(unit_cell, params);
    auto d1 = TestAccessor::D1(structure_constant, 0, 0, k, 0.1 + 0.1i);
    auto d2 = TestAccessor::D2(structure_constant, 0, 0, k, 0.1 + 0.1i);
    auto d3 = TestAccessor::D3(structure_constant, 0, 0, k, 0.1 + 0.1i);
    auto res = structure_constant.calculate_reciprocal_space(0, 0, 0, 0, k, 0.1 + 0.1i);
    std::cout
      //<< "scale " << std::setw(10) << scale
      << " k = " << std::scientific << std::setprecision(14) << std::to_string(k)
      << " d1 = " << d1
      << " d2 = " << d2
      << " d3 = " << d3
      << " res = " << res
      << std::endl;
  }
}
