#include <tuple>

#include <Physics/NonRelativistic/StructureConstants.h>
#include <TestUtils/base.h>
#include <Geometry/IrreducibleWedge.h>
#include <fstream>
#include <TestUtils/Utils.h>

using namespace Physics::NonRelativistic;
using namespace std::complex_literals;

namespace
{
}

TEST(TestStructureConstants, Real)
{
  auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto structure_constant = RealStructureConstants(unit_cell);
  auto nothing = 0i;
  for (auto rez = 0.01; rez < 1.0; rez += 0.1) {
    for (auto n = 0; n < 4; ++n) {
      for (auto m = 0; m < 4; ++m) {
        for (auto l = 0; l < 4; ++l) {
          if (!(m == 0 && n == 0 && l == 0)) {
            auto res = structure_constant.calculate(
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
        const ReciprocalStructureConstantsCalculator &structure_constants, unsigned int l, int m, const Vector3D &k
      )
      {
        return structure_constants.D1(l, m, k);
      }

      static std::complex<double> D2(
        const ReciprocalStructureConstantsCalculator &structure_constants, unsigned int l, int m, const Vector3D &k
      )
      {
        return structure_constants.D2(l, m, k);
      }

      static std::complex<double> D3(
        const ReciprocalStructureConstantsCalculator &structure_constants, unsigned int l, int m
      )
      {
        return structure_constants.D3(l, m);
      }
    };
  }
}

TEST(TestStructureConstants, Reciprocal)
{
  auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto reciprocal_unit_cell = ReciprocalUnitCell3D(unit_cell);
  auto irreduc_wedge = IrreducibleWedge::get_irreducible_wedge(reciprocal_unit_cell, 5, true);
  auto nothing = 0i;

  for (auto k: irreduc_wedge) {
    auto params = StructureConstantsConfig(
      10.0,
      3.0
    );
    auto structure_constant = ReciprocalStructureConstants(unit_cell, params);
    auto calculator = structure_constant.get_calculator(0, 0.1 + 0.1i);
    auto res = calculator.calculate(0, 0, 0, 0, k);
    nothing += res;

//std::cout
//<< "scale " << std::setw(10) << scale
//      << " k = " << std::scientific << std::setprecision(14) << std::to_string(k)
//      << " d1 = " << d1
//      << " d2 = " << d2
//      << " d3 = " << d3
//      << " res = " << res
//      << std::endl;
  }
}

namespace
{
  void log(const std::vector<std::pair<Point3D, std::complex<double>>> &R, const std::string &filename)
  {
    std::ofstream out_R;
    out_R.open(filename + ".dat");
    for (auto p : R)
      out_R << std::setfill(' ') << std::setw(20) << std::setprecision(17) << std::fixed
            << p.first.x << '\t' << p.first.y << '\t' << p.first.z << '\t'
            << p.second.real() << '\t' << p.second.imag() << "\n";
    out_R.close();
  }
}

TEST(DISABLED_TestStructureConstants, GenerateSlice)
{
  auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto reciprocal_unit_cell = ReciprocalUnitCell3D(UnitCell3D::create_cubic_primitive(4.0));
  std::cout << "start to generate irreduc wedge" << std::endl;
  // TODO cannot use replicate that simply because spherical harmonics have their own transformations
  //auto irreduc_wedge = IrreducibleWedge::get_irreducible_wedge(reciprocal_unit_cell, 10, true);
  auto irreduc_wedge = CubicMesh(reciprocal_unit_cell.v1.length() / 10).generate(CutoffWSCell(reciprocal_unit_cell),
                                                                                 true);
  Utils::log(irreduc_wedge, "irreduc_wedge");
  std::cout << "start to generate irreduc wedge... done" << std::endl;

  auto params = StructureConstantsConfig(10.0, 3.0);
  std::cout << "start to generate ReciprocalStructureConstants" << std::endl;
  auto structure_constant = ReciprocalStructureConstants(unit_cell, params);
  std::cout << "start to generate ReciprocalStructureConstants... done" << std::endl;

  std::cout << "start to generate calculator" << std::endl;
  auto l_max = 2u;
  auto calculator = structure_constant.get_calculator(l_max, 0.1 + 0.1i);
  std::cout << "start to generate calculator... done" << std::endl;

  auto size = irreduc_wedge.size();

  for (auto l = 0u; l <= l_max; ++l) {
    for (int m = -l; m <= static_cast<int>(l); ++m) {
      std::cout << "start to calculate struct consts l = " << l << " m = " << m << std::endl;
      auto i = 0;
      auto irreduc_res = std::vector<std::pair<Point3D, std::complex<double>>>();
      irreduc_res.reserve(irreduc_wedge.size() / 24);
      for (auto k: irreduc_wedge) {
        if ((i++ % 1000) == 0) {
          std::cout << std::fixed << std::setw(6) << i << " / " << size << std::endl;
        }
        auto res = calculator.calculate(l, m, 0, 0, k);
        irreduc_res.push_back(std::make_pair(k, res));
      }
      std::cout << "start to calculate struct consts... done" << std::endl;
      auto transformations =
        SymmetryTransformationFactory::get(
          CrystallographicPointGroup::create(unit_cell.get_point_group())->get_elements());

      std::cout << "start to replicate" << std::endl;
      //auto replicated_res = IrreducibleWedge::replicate(irreduc_res, transformations);
      std::cout << "start to replicate... done" << std::endl;
      auto filename = "rec_struct_z_l" + std::to_string(l) + "_m" + std::to_string(m);
      std::cout << "start to write to file " << filename << std::endl;
      log(irreduc_res, filename);
      std::cout << "start to write to file " << filename << "... done" << std::endl;
    }
  }
}
