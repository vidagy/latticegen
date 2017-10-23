#ifndef LATTICEGEN_VISITOR_H
#define LATTICEGEN_VISITOR_H

#include <memory>
#include <utility>

namespace Core
{
  template<typename V, typename T, typename... Ts>
  struct VisitorAdapter;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// @brief free standing apply_visitor -- primary entry point for visitors
  ///
  /// This is a peculiar implementation of the Visitor pattern. We decoupled the application of the visitor from the
  /// process of dispatching its arguments to specific implementations.
  ///
  /// Core::apply_visitor(v, t) can be used for applying a visitor on a single visitable. The implementor of the class
  /// hierarchy T should specialize Core::VisitorAdapter<V, T>, so that its operator(const T&) can recognize the
  /// concrete T types and can call V.operator() with the specific types, T1, T2, etc.
  ///
  /// You can call apply_visitor with multiple arguments, Core::apply_visitor(v, t1, t2). Here, we have a default
  /// implementation for the corresponding multi-argument VisitorAdapter. It will dispatch each visitable one by one
  /// using the one-argument VisitorAdapters. So if the user of this pattern implemented VisitorAdapter once for a
  /// single argument, visitors with arbitrary number of arguments will also be dispatched out of the box.
  template<typename V, typename... T>
  typename V::result_type apply_visitor(const V &v, const T &... t)
  {
    return Core::VisitorAdapter<V, T...>::apply_visitor(v, t...);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // detail for dispatching multiple visitables
  namespace detail
  {
    struct EndType
    {
    };

    template<typename V, typename... Ts>
    struct VisitorAdapterVisitor
    {
      typedef typename V::result_type result_type;

      VisitorAdapterVisitor(const V &v, const Ts &... arguments)
        : v(v), arguments(std::make_tuple(std::cref(arguments)...)) {}

    private:
      template<typename VV, typename TT, typename T_ignored, typename... TTs, size_t... Is>
      static VisitorAdapterVisitor<VV, TTs..., TT> create_visitor_without_first(
        const VV &v,
        const std::tuple<std::reference_wrapper<const T_ignored>, std::reference_wrapper<const TTs>...> &arguments,
        const TT &t,
        std::index_sequence<Is...>
      )
      {
        return VisitorAdapterVisitor<V, TTs..., TT>(v, (std::get<Is + 1>(arguments).get())..., t);
      };

      template<typename VV, typename TT, typename T_ignored, typename... TTs,
        typename Indices = std::make_index_sequence<std::tuple_size<std::tuple<TTs...>>::value>>
      static VisitorAdapterVisitor<VV, TTs..., TT> create_visitor_without_first(
        const VV &v,
        const std::tuple<std::reference_wrapper<const T_ignored>, std::reference_wrapper<const TTs>...> &arguments,
        const TT &t
      )
      {
        return VisitorAdapterVisitor<V, Ts...>::create_visitor_without_first(v, arguments, t, Indices{});
      };

    public:
      template<typename T>
      result_type operator()(const T &t) const
      {
        auto visitor = VisitorAdapterVisitor<V, Ts...>::create_visitor_without_first(v, arguments, t);
        return Core::apply_visitor(visitor, std::get<0>(arguments).get());
      }

    private:
      const V &v;
      std::tuple<std::reference_wrapper<const Ts>...> arguments;
    };

    template<typename V, typename... TDs>
    struct VisitorAdapterVisitor<V, EndType, TDs...>
    {
      typedef typename V::result_type result_type;

      VisitorAdapterVisitor(const V &v, EndType, const TDs &... arguments)
        : v(v), arguments(std::make_tuple(std::cref(arguments)...)) {}

      // all arguments dispatched to the concrete types
      template<typename T, typename Indices = std::make_index_sequence<std::tuple_size<std::tuple<TDs...>>::value>>
      result_type operator()(const T &t) const
      {
        return (*this)(t, Indices{});
      }

    private:
      template<typename T, size_t... Is>
      result_type operator()(const T &t, std::index_sequence<Is...>) const
      {
        return v((std::get<Is>(arguments).get())..., t);
      }

      const V &v;
      std::tuple<std::reference_wrapper<const TDs>...> arguments;
    };
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // generic implementation for dispatching multiple visitable
  template<typename V, typename T, typename... Ts>
  struct VisitorAdapter
  {
    static typename V::result_type apply_visitor(const V &v, const T &t, const Ts &... ts)
    {
      // dispatching the first argument to the
      return Core::apply_visitor(
        Core::detail::VisitorAdapterVisitor<V, Ts..., Core::detail::EndType>(v, ts..., Core::detail::EndType{}), t);
    }
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // generic implementation for dispatching a single visitable, this should be specialized for each visitable hierarchy.
  template<typename V, typename T>
  struct VisitorAdapter<V, T>
  {
    static typename V::result_type apply_visitor(const V &v, const T &t)
    {
      return v(t);
    }
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // specialization for managed_ptrs
  template<typename V, typename T>
  struct VisitorAdapter<V, std::unique_ptr<T>>
  {
    static typename V::result_type apply_visitor(const V &v, const std::unique_ptr<T> &t)
    {
      assert(t);
      return Core::apply_visitor(v, *t);
    }
  };

  template<typename V, typename T>
  struct VisitorAdapter<V, std::shared_ptr<T>>
  {
    static typename V::result_type apply_visitor(const V &v, const std::shared_ptr<T> &t)
    {
      assert(t);
      return Core::apply_visitor(v, *t);
    }
  };
}

#endif //LATTICEGEN_VISITOR_H
