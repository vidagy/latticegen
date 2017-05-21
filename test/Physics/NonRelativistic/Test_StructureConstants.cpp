#include <tuple>

#include <Physics/NonRelativistic/StructureConstants.h>
#include <TestUtils/base.h>

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

TEST(TestStructureConstants, Reciprocal)
{
  auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto structure_constant = StructureConstants(unit_cell);
  auto nothing = 0i;
  auto k = Vector3D(0.144, 0.42, 0.137);
  for (auto rez = 0.01; rez < 0.1; rez += 0.1) {
//    for (auto n = 0; n < 4; ++n) {
//      for (auto m = 0; m < 4; ++m) {uuu
//        for (auto l = 0; l < 4; ++l) {
//          if (!(m == 0 && n == 0 && l == 0)) {
    auto res = structure_constant.calculate_reciprocal_space(
      0, 0, 0, 0, k, rez + 0.1i
    );
    std::cout << k.length()
              << "  " << res << std::endl;
    nothing += res;
//          }
//        }
//      }
//    }
  }
}
