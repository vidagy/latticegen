
#include <Core/Point3D.h>
#include <vector>
#include "PolyhedralQuadrature.h"

using namespace Core;
namespace Math
{
  class LebedevQuadrature
  {
  public:

    enum class Order
    {
      LD0006 = 0, LD0014, LD0026, LD0038, LD0050, LD0074, LD0086, LD0110, LD0146, LD0170, LD0194, LD0230, LD0266,
      LD0302, LD0350, LD0434, LD0590, LD0770, LD0974, LD1202, LD1454, LD1730, LD2030, LD2354, LD2702, LD3074, LD3470,
      LD3890, LD4334, LD4802, LD5294, LD5810
    };

    static Quadrature generate(Order order);

    static unsigned int order_to_uint(Order order);
  };
}

