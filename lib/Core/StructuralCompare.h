#ifndef LATTICEGEN_STRUCTURALCOMPARE_H_
#define LATTICEGEN_STRUCTURALCOMPARE_H_

#include <type_traits>
#include <memory>
#include <set>
#include <vector>
#include <list>
#include <map>
#include <tuple>

namespace Core
{
  enum class Cmp_Result
  {
    Cmp_Less, Cmp_Eq, Cmp_Greater
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // free standing structural_compare functions

  namespace detail
  {
    // helper functions to dispatch structural compare for classes that implement it and to intrinsic types that don't

    // specialization for structs/classes implementing "Cmp_Result structural_compare(const T& other)"
    template<typename F, typename T>
    inline typename std::enable_if<std::is_same<F, std::true_type>::value, Cmp_Result>::type
    structural_compare(const T &lhs, const T &rhs)
    {
      // memory-wise equivalency
      if (&lhs == &rhs)
        return Cmp_Result::Cmp_Eq;

      return lhs.structural_compare(rhs);
    }

    // specializations for intrinsic types
    template<typename F, typename T>
    inline typename std::enable_if<std::is_same<F, std::false_type>::value, Cmp_Result>::type
    structural_compare(const T &lhs, const T &rhs)
    {
      if (lhs < rhs)
        return Cmp_Result::Cmp_Less;
      else if (lhs > rhs)
        return Cmp_Result::Cmp_Greater;
      else
        return Cmp_Result::Cmp_Eq;
    }

    // predicate on having structural_compare implemented by T
    template<typename T>
    struct has_structural_compare
    {
    private:
      template<typename TT>
      static constexpr auto check(int) ->
      typename std::is_same<decltype(std::declval<TT>().structural_compare(std::declval<TT>())), Cmp_Result>::type;

      template<typename>
      static constexpr auto check(long) -> std::false_type;

    public:
      typedef decltype(check<T>(0)) type;
      static constexpr bool value = type::value;
    };
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // primary entry point for structural_scompare
  template<typename T>
  inline Cmp_Result structural_compare(const T &lhs, const T &rhs)
  {
    return detail::structural_compare<typename detail::has_structural_compare<T>::type>(lhs, rhs);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // specialization for managed pointers
  template<typename T>
  inline Cmp_Result structural_compare(const std::shared_ptr<T> &lhs, const std::shared_ptr<T> &rhs)
  {
    return lhs->structural_compare(*rhs);
  }

  template<typename T>
  inline Cmp_Result structural_compare(const std::unique_ptr<T> &lhs, const std::unique_ptr<T> &rhs)
  {
    return lhs->structural_compare(*rhs);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // specialization for containers
  template<typename T1, typename T2>
  inline Cmp_Result structural_compare(const std::pair<T1, T2> &lhs, const std::pair<T1, T2> &rhs)
  {
    auto res = structural_compare(lhs.first, rhs.first);
    if (res != Cmp_Result::Cmp_Eq)
      return res;
    else
      return structural_compare(lhs.second, rhs.second);
  }

  namespace detail
  {
    template<typename ItT>
    inline Cmp_Result structural_compare_range(ItT lhs_it, const ItT &lhs_end, ItT rhs_it, const ItT &rhs_end)
    {
      for (; lhs_it != lhs_end && rhs_it != rhs_end; ++lhs_it, ++rhs_it) {
        auto res = Core::structural_compare(*lhs_it, *rhs_it);
        if (res != Cmp_Result::Cmp_Eq)
          return res;
      }
      if (lhs_it != lhs_end)
        return Cmp_Result::Cmp_Greater;

      if (rhs_it != rhs_end)
        return Cmp_Result::Cmp_Less;

      return Cmp_Result::Cmp_Eq;
    }
  }

  template<typename T>
  inline Cmp_Result structural_compare(const std::set<T> &lhs, const std::set<T> &rhs)
  {
    return detail::structural_compare_range(lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend());
  }

  template<typename T>
  inline Cmp_Result structural_compare(const std::vector<T> &lhs, const std::vector<T> &rhs)
  {
    return detail::structural_compare_range(lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend());
  }

  template<typename T>
  inline Cmp_Result structural_compare(const std::list<T> &lhs, const std::list<T> &rhs)
  {
    return detail::structural_compare_range(lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend());
  }

  template<typename K, typename V>
  inline Cmp_Result structural_compare(const std::map<K, V> &lhs, const std::map<K, V> &rhs)
  {
    return detail::structural_compare_range(lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend());
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // StructuralComparable interface
  template<typename T>
  struct StructuralComparable
  {

    inline Cmp_Result structural_compare(const T &other) const
    {
      auto members = T::get_members();
      return structural_compare_impl(other, members);
    }

    inline bool operator==(const T &other) const
    {
      return Cmp_Result::Cmp_Eq == structural_compare(other);
    }

    inline bool operator!=(const T &other) const
    {
      return Cmp_Result::Cmp_Eq != structural_compare(other);
    }

    inline bool operator<(const T &other) const
    {
      return Cmp_Result::Cmp_Less == structural_compare(other);
    }

    inline bool operator<=(const T &other) const
    {
      return Cmp_Result::Cmp_Greater != structural_compare(other);
    }

    inline bool operator>(const T &other) const
    {
      return Cmp_Result::Cmp_Greater == structural_compare(other);
    }

    inline bool operator>=(const T &other) const
    {
      return Cmp_Result::Cmp_Less != structural_compare(other);
    }

  protected:
    // helper function #1 to unwrap the member tuple
    template<typename... MV, typename Indices = std::make_index_sequence<std::tuple_size<std::tuple<MV...>>::value>>
    inline Cmp_Result structural_compare_impl(const T &right, const std::tuple<MV...> &tuple) const
    {
      return structural_compare_impl(right, tuple, Indices{});
    }

    // helper function #2 to unwrap the member tuple
    template<typename... MV, size_t... Is>
    inline Cmp_Result structural_compare_impl(
      const T &right, const std::tuple<MV...> &tuple, std::index_sequence<Is...>
    ) const
    {
      return structural_compare_members(right, std::get<Is>(tuple)...);
    }

    // actual implementation of comparison
    template<typename... MVs>
    inline Cmp_Result structural_compare_members(MVs...) const
    {
      return Cmp_Result::Cmp_Eq;
    }

    template<typename R, typename... MVs>
    inline Cmp_Result structural_compare_members(const T &right, R T::*mv, MVs... rest) const
    {
      R lhs = static_cast<const T *>(this)->*mv;
      R rhs = right.*mv;

      Cmp_Result res = Core::structural_compare(lhs, rhs);

      if (res != Cmp_Result::Cmp_Eq)
        return res;
      else {
        return structural_compare_members(right, rest...);
      }
    }
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Structural comparator operators

  template<typename T>
  struct StructuralLessOp
  {
    bool operator()(const T &lhs, const T &rhs) const
    {
      return structural_compare(lhs, rhs) == Cmp_Result::Cmp_Less;
    }
  };

  template<typename T>
  struct StructuralEqOp
  {
    bool operator()(const T &lhs, const T &rhs) const
    {
      return structural_compare(lhs, rhs) == Cmp_Result::Cmp_Eq;
    }
  };

  template<typename T> using StructuralSet = std::set<T, StructuralLessOp<T>>;
  template<typename K, typename V> using StructuralMap = std::map<K, V, StructuralLessOp<K>>;
}

#endif
