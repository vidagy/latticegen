#include <TestUtils/base.h>

#include <Math/LebedevQuadrature.h>
#include <fstream>

using namespace Math;

namespace
{
  const static std::vector<LebedevQuadrature::Order> all_quadratures = {
    LebedevQuadrature::Order::LD0006, LebedevQuadrature::Order::LD0014, LebedevQuadrature::Order::LD0026,
    LebedevQuadrature::Order::LD0038, LebedevQuadrature::Order::LD0050, LebedevQuadrature::Order::LD0074,
    LebedevQuadrature::Order::LD0086, LebedevQuadrature::Order::LD0110, LebedevQuadrature::Order::LD0146,
    LebedevQuadrature::Order::LD0170, LebedevQuadrature::Order::LD0194, LebedevQuadrature::Order::LD0230,
    LebedevQuadrature::Order::LD0266, LebedevQuadrature::Order::LD0302, LebedevQuadrature::Order::LD0350,
    LebedevQuadrature::Order::LD0434, LebedevQuadrature::Order::LD0590, LebedevQuadrature::Order::LD0770,
    LebedevQuadrature::Order::LD0974, LebedevQuadrature::Order::LD1202, LebedevQuadrature::Order::LD1454,
    LebedevQuadrature::Order::LD1730, LebedevQuadrature::Order::LD2030, LebedevQuadrature::Order::LD2354,
    LebedevQuadrature::Order::LD2702, LebedevQuadrature::Order::LD3074, LebedevQuadrature::Order::LD3470,
    LebedevQuadrature::Order::LD3890, LebedevQuadrature::Order::LD4334, LebedevQuadrature::Order::LD4802,
    LebedevQuadrature::Order::LD5294, LebedevQuadrature::Order::LD5810
  };

  void print_lebedev_quadrature(LebedevQuadrature::Order order)
  {
    auto res = LebedevQuadrature::generate(order);

    std::ofstream out_R;
    out_R.open("Lebedev_" + std::to_string(LebedevQuadrature::order_to_uint(order)) + ".dat");
    for (auto r : res)
      out_R << std::setw(20) << std::setprecision(17) << std::fixed
            << r.first.x << "\t" << r.first.y << "\t" << r.first.z << "\t" << r.second << "\n";
    out_R.close();
  }
}

TEST(DISABLED_LebedevQuadrature, Print)
{
  std::for_each(all_quadratures.begin(), all_quadratures.end(), print_lebedev_quadrature);
}

TEST(LebedevQuadrature, NumberOfPoints)
{
  std::for_each(
    all_quadratures.begin(), all_quadratures.end(),
    [](LebedevQuadrature::Order order)
    {
      EXPECT_EQ(LebedevQuadrature::generate(order).size(), LebedevQuadrature::order_to_uint(order));
    }
  );
}

TEST(LebedevQuadrature, SumWeights)
{
  std::for_each(
    all_quadratures.begin(), all_quadratures.end(),
    [](LebedevQuadrature::Order order)
    {
      auto res = LebedevQuadrature::generate(order);
      auto sum_weights = std::accumulate(
        res.begin(), res.end(), 0.0,
        [](double sum, const std::pair<Point3D, double> &p)
        {
          return sum + p.second;
        }
      );
      EXPECT_NEAR(sum_weights, 1.0, 3e-14);
    }
  );
}