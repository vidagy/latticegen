#include <TestUtils/base.h>

#include <Core/StructuralCompare.h>

using namespace Core;

TEST(TestStructuralCompare, BaseTypes)
{
  // int
  EXPECT_EQ(structural_compare(2, 3), Cmp_Result::Cmp_Less);
  EXPECT_EQ(structural_compare(2, 2), Cmp_Result::Cmp_Eq);
  EXPECT_EQ(structural_compare(4, 3), Cmp_Result::Cmp_Greater);

  // long
  EXPECT_EQ(structural_compare(2l, 3l), Cmp_Result::Cmp_Less);
  EXPECT_EQ(structural_compare(2l, 2l), Cmp_Result::Cmp_Eq);
  EXPECT_EQ(structural_compare(4l, 3l), Cmp_Result::Cmp_Greater);

  // double
  EXPECT_EQ(structural_compare(2.0, 3.0), Cmp_Result::Cmp_Less);
  EXPECT_EQ(structural_compare(2.0, 2.0), Cmp_Result::Cmp_Eq);
  EXPECT_EQ(structural_compare(4.0, 3.0), Cmp_Result::Cmp_Greater);
}

namespace
{
  struct A : public StructuralComparable<A>
  {
    explicit A(int a) : a(a) {};

    static constexpr auto get_members()
    {
      return std::make_tuple(&A::a);
    }

    int a;
  };
}

TEST(TestStructuralCompare, SimpleStruct)
{
  auto a = A(1);
  auto b = A(2);
  EXPECT_EQ(structural_compare(a, b), Cmp_Result::Cmp_Less);
  EXPECT_EQ(structural_compare(b, a), Cmp_Result::Cmp_Greater);
  EXPECT_EQ(structural_compare(a, a), Cmp_Result::Cmp_Eq);
  EXPECT_EQ(structural_compare(b, b), Cmp_Result::Cmp_Eq);
}

namespace
{
  struct B : public StructuralComparable<B>
  {
    B(double b, A a) : b(b), a(a) {};

    static constexpr auto get_members()
    {
      return std::make_tuple(&B::b, &B::a);
    }

    double b;
    A a;
  };

  class C : public StructuralComparable<C>
  {
  public:
    C(bool c, B b) : c(c), b(b) {};


    static constexpr auto get_members()
    {
      return std::make_tuple(&C::c, &C::b);
    }

  private:
    bool c;
    B b;
  };

}

TEST(TestStructuralCompare, CompositeStruct)
{
  auto a1 = A(1);
  auto b1 = B(2.0, a1);
  auto c1 = C(true, b1);

  auto a2 = A(3);
  auto b2 = B(2.0, a2);
  auto c2 = C(true, b2);

  EXPECT_EQ(structural_compare(c1, c2), Cmp_Result::Cmp_Less);
  EXPECT_EQ(structural_compare(c2, c2), Cmp_Result::Cmp_Eq);
  EXPECT_EQ(structural_compare(c1, c1), Cmp_Result::Cmp_Eq);
}

TEST(TestStructuralCompare, CompositeInSharedPtr)
{
  auto a1 = A(1);
  auto b1 = B(2.0, a1);
  auto c1p = std::make_shared<C>(true, b1);

  auto a2 = A(3);
  auto b2 = B(2.0, a2);
  auto c2p = std::make_shared<C>(true, b2);

  EXPECT_EQ(structural_compare(c1p, c2p), Cmp_Result::Cmp_Less);
  EXPECT_EQ(structural_compare(c2p, c2p), Cmp_Result::Cmp_Eq);
  EXPECT_EQ(structural_compare(c1p, c1p), Cmp_Result::Cmp_Eq);
}

TEST(TestStructuralCompare, CompositeInUniquePtr)
{
  auto a1 = A(1);
  auto b1 = B(2.0, a1);
  auto c1p = std::make_unique<C>(true, b1);

  auto a2 = A(3);
  auto b2 = B(2.0, a2);
  auto c2p = std::make_unique<C>(true, b2);

  EXPECT_EQ(structural_compare(c1p, c2p), Cmp_Result::Cmp_Less);
  EXPECT_EQ(structural_compare(c2p, c2p), Cmp_Result::Cmp_Eq);
  EXPECT_EQ(structural_compare(c1p, c1p), Cmp_Result::Cmp_Eq);
}

TEST(TestStructuralCompare, Vector)
{
  auto a1 = std::vector<int>{3, 1, 2};
  auto a2 = std::vector<int>{3, 1, 3};
  EXPECT_EQ(structural_compare(a1, a2), Cmp_Result::Cmp_Less);

  auto b1 = std::vector<int>{3, 1, 2};
  auto b2 = std::vector<int>{3, 1};
  EXPECT_EQ(structural_compare(b1, b2), Cmp_Result::Cmp_Greater);
  EXPECT_EQ(structural_compare(b2, b1), Cmp_Result::Cmp_Less);

  auto c1 = std::vector<A>{A(3), A(1), A(2)};
  auto c2 = std::vector<A>{A(3), A(1), A(3)};
  EXPECT_EQ(structural_compare(c1, c2), Cmp_Result::Cmp_Less);

  EXPECT_EQ(structural_compare(c1, c1), Cmp_Result::Cmp_Eq);
}


TEST(TestStructuralCompare, Pair)
{
  auto a1 = A(1);
  auto b1 = B(2.0, a1);

  auto a2 = A(3);
  auto b2 = B(2.0, a2);

  auto p1 = std::make_pair(a1, b1);
  auto p2 = std::make_pair(a1, b2);
  auto p3 = std::make_pair(a2, b1);
  EXPECT_EQ(structural_compare(p1, p2), Cmp_Result::Cmp_Less);
  EXPECT_EQ(structural_compare(p3, p2), Cmp_Result::Cmp_Greater);
}

TEST(TestStructuralCompare, Map)
{
  auto a1 = A(1);
  auto a2 = A(2);
  auto a3 = A(3);
  auto b1 = B(3.0, a1);
  auto b2 = B(2.0, a2);
  auto b3 = B(1.0, a3);
  auto map1 = std::map<int, B>{{1, b1},
                               {2, b2},
                               {3, b3}};
  auto map2 = std::map<int, B>{{2, b1},
                               {1, b2},
                               {3, b3}};
  EXPECT_EQ(structural_compare(map1, map2), Cmp_Result::Cmp_Greater);
}