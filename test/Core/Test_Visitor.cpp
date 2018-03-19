#include <TestUtils/base.h>

#include <Core/Visitor.h>

using namespace Core;

namespace
{
  struct A
  {
    enum class Type
    {
      B, C
    };
    Type type;
  protected:
    explicit A(Type type) : type(type) {}

    virtual ~A() {}
  };

  struct B : A
  {
    B() : A(Type::B) {};
  };

  struct C : A
  {
    C() : A(Type::C) {};
  };
}

namespace Core
{
  template<typename V>
  struct VisitorAdapter<V, A>
  {
    static typename V::result_type apply_visitor(const V &v, const A &a)
    {
      switch (a.type) {
        case A::Type::B:
          return Core::apply_visitor(v, static_cast<const B &>(a));
        case A::Type::C:
          return Core::apply_visitor(v, static_cast<const C &>(a));
      }
      throw std::logic_error("Cannot get here, switch case above is exhaustive.");
    }
  };
}

namespace
{
  struct GetName
  {
    typedef std::string result_type;

    LATTICEGEN_MUTE_BEGIN
    LATTICEGEN_MUTE_UNUSED_VAR
    result_type operator()(const B &b) const
    {
      LATTICEGEN_MUTE_END
      return "B";
    }

    LATTICEGEN_MUTE_BEGIN
    LATTICEGEN_MUTE_UNUSED_VAR
    result_type operator()(const C &c) const
    {
      LATTICEGEN_MUTE_END
      return "C";
    }
  };
}

TEST(Visitor, EnumBased)
{
  auto b = B();
  auto c = C();

  EXPECT_EQ(Core::apply_visitor(GetName(), b), "B");
  EXPECT_EQ(Core::apply_visitor(GetName(), c), "C");

  A &ba = b;
  A &ca = c;

  EXPECT_EQ(Core::apply_visitor(GetName(), ba), "B");
  EXPECT_EQ(Core::apply_visitor(GetName(), ca), "C");

  auto bs = std::make_shared<B>();
  auto cs = std::make_shared<C>();

  EXPECT_EQ(Core::apply_visitor(GetName(), bs), "B");
  EXPECT_EQ(Core::apply_visitor(GetName(), cs), "C");
}

namespace
{
  struct D
  {
    virtual void foo() = 0; // so that D cannot be instantiated
  };

  struct E : D
  {
    E() = default;

    E(const E &) = delete;

    E(E &&) = delete;

    E &operator=(const E) = delete;

    E &operator=(E &&) = delete;

    virtual void foo() override {}
  };

  struct F : D
  {
    F() = default;

    F(const F &) = delete;

    F(F &&) = delete;

    F &operator=(const F)  = delete;

    F &operator=(F &&) = delete;

    virtual void foo() override {}
  };
}

namespace Core
{
  template<typename V>
  struct VisitorAdapter<V, D>
  {
    static typename V::result_type apply_visitor(const V &v, const D &d)
    {
      if (dynamic_cast<const E *>(&d))
        return Core::apply_visitor(v, static_cast<const E &>(d));
      if (dynamic_cast<const F *>(&d))
        return Core::apply_visitor(v, static_cast<const F &>(d));

      throw std::logic_error("cannot check for exhaustiveness");
    }
  };
}

namespace
{
  struct GetName2
  {
    typedef std::string result_type;

    LATTICEGEN_MUTE_BEGIN
    LATTICEGEN_MUTE_UNUSED_VAR
    result_type operator()(const E &e) const
    {
      LATTICEGEN_MUTE_END
      return "E";
    }

    LATTICEGEN_MUTE_BEGIN
    LATTICEGEN_MUTE_UNUSED_VAR
    result_type operator()(const F &f) const
    {
      LATTICEGEN_MUTE_END
      return "F";
    }
  };
}

TEST(Visitor, DynamicCastBased)
{
  E e{};
  F f{};

  EXPECT_EQ(Core::apply_visitor(GetName2(), e), "E");
  EXPECT_EQ(Core::apply_visitor(GetName2(), f), "F");

  D &ed = e;
  D &fd = f;

  EXPECT_EQ(Core::apply_visitor(GetName2(), ed), "E");
  EXPECT_EQ(Core::apply_visitor(GetName2(), fd), "F");

  auto es = std::make_shared<E>();
  auto fs = std::make_shared<F>();

  EXPECT_EQ(Core::apply_visitor(GetName2(), es), "E");
  EXPECT_EQ(Core::apply_visitor(GetName2(), fs), "F");
}

namespace
{
  struct GetTwoNames
  {
    typedef std::string result_type;

    result_type operator()(const E &, const E &) const
    {
      return "EE";
    }

    result_type operator()(const E &, const F &) const
    {
      return "EF";
    }

    result_type operator()(const F &, const E &) const
    {
      return "FE";
    }

    result_type operator()(const F &, const F &) const
    {
      return "FF";
    }
  };
}

TEST(Visitor, Multiargument)
{
  E e{};
  F f{};

  EXPECT_EQ(Core::apply_visitor(GetTwoNames(), e, f), "EF");
  EXPECT_EQ(Core::apply_visitor(GetTwoNames(), e, e), "EE");
  EXPECT_EQ(Core::apply_visitor(GetTwoNames(), f, e), "FE");
  EXPECT_EQ(Core::apply_visitor(GetTwoNames(), f, f), "FF");

  D &ed = e;
  D &fd = f;

  EXPECT_EQ(Core::apply_visitor(GetTwoNames(), ed, ed), "EE");
  EXPECT_EQ(Core::apply_visitor(GetTwoNames(), ed, fd), "EF");
  EXPECT_EQ(Core::apply_visitor(GetTwoNames(), fd, ed), "FE");
  EXPECT_EQ(Core::apply_visitor(GetTwoNames(), fd, fd), "FF");

  auto es = std::make_shared<E>();
  auto fs = std::make_shared<F>();

  EXPECT_EQ(Core::apply_visitor(GetTwoNames(), es, es), "EE");
  EXPECT_EQ(Core::apply_visitor(GetTwoNames(), es, fs), "EF");
  EXPECT_EQ(Core::apply_visitor(GetTwoNames(), fs, es), "FE");
  EXPECT_EQ(Core::apply_visitor(GetTwoNames(), fs, fs), "FF");
}