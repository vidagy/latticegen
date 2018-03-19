#include <TestUtils/base.h>

#include <Geometry/SymmetryTransformationFactory.h>

using namespace Geometry;

TEST(TestSymmetryTransformationFactory,Identity)
{
  EXPECT_EQ(
    Identity(),
    SymmetryTransformationFactory::get(CrystallographicPointGroup::Identity)
  );
}

TEST(TestSymmetryTransformationFactory,Ci)
{
  const static SymmetryTransformationFactory::Transformations reference{Identity(), Inversion()};
  EXPECT_THAT(
    reference,
    ::testing::UnorderedElementsAreArray(SymmetryTransformationFactory::get(Ci().get_elements())) );
}

namespace
{
  template<class Group>
  void testCrystallographicPointGroup()
  {
    const SymmetryTransformationFactory::Transformations elements =
      SymmetryTransformationFactory::get(Group().get_elements());
    const CrystallographicPointGroup::Elements generators = Group().get_generators();
    const SymmetryTransformationFactory::Transformations generated_elements =
      SymmetryTransformationFactory::generate(generators);

    EXPECT_EQ(elements.size(), generated_elements.size());

    std::for_each(
      generated_elements.begin(),
      generated_elements.end(),
      [&elements](const Transformation &generated_element)
      {
        auto it = std::find_if(
          elements.begin(),
          elements.end(),
          [&generated_element](const Transformation &element)
          {
            return element.transformation_matrix.isApprox(generated_element.transformation_matrix);
          });

        if (it == elements.end())
          std::cout << "generated element is NOT in elements " << std::endl << generated_element << std::endl;

        EXPECT_NE(it, elements.end());
      }
    );
  }
}

TEST(TestSymmetryTransformationFactory,GeneratorsC1)
{
  testCrystallographicPointGroup<C1>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsCi)
{
  testCrystallographicPointGroup<Ci>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsC2)
{
  testCrystallographicPointGroup<C2>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsCs)
{
  testCrystallographicPointGroup<Cs>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsC2h)
{
  testCrystallographicPointGroup<C2h>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsD2)
{
  testCrystallographicPointGroup<D2>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsC2v)
{
  testCrystallographicPointGroup<C2v>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsD2h)
{
  testCrystallographicPointGroup<D2h>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsC4)
{
  testCrystallographicPointGroup<C4>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsS4)
{
  testCrystallographicPointGroup<S4>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsC4h)
{
  testCrystallographicPointGroup<C4h>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsD4)
{
  testCrystallographicPointGroup<D4>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsC4v)
{
  testCrystallographicPointGroup<C4v>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsD2d)
{
  testCrystallographicPointGroup<D2d>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsD4h)
{
  testCrystallographicPointGroup<D4h>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsC3)
{
  testCrystallographicPointGroup<C3>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsS6)
{
  testCrystallographicPointGroup<S6>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsD3)
{
  testCrystallographicPointGroup<D3>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsC3v)
{
  testCrystallographicPointGroup<C3v>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsD3d)
{
  testCrystallographicPointGroup<D3d>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsC6)
{
  testCrystallographicPointGroup<C6>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsC3h)
{
  testCrystallographicPointGroup<C3h>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsC6h)
{
  testCrystallographicPointGroup<C6h>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsD6)
{
  testCrystallographicPointGroup<D6>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsC6v)
{
  testCrystallographicPointGroup<C6v>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsD3h)
{
  testCrystallographicPointGroup<D3h>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsD6h)
{
  testCrystallographicPointGroup<D6h>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsT)
{
  testCrystallographicPointGroup<T>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsTh)
{
  testCrystallographicPointGroup<Th>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsO)
{
  testCrystallographicPointGroup<O>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsTd)
{
  testCrystallographicPointGroup<Td>();
}

TEST(TestSymmetryTransformationFactory,GeneratorsOh)
{
  testCrystallographicPointGroup<Oh>();
}
