#pragma once

#include <concepts>
#include <ranges>
#include <tuple>
#include <type_traits>

namespace jr_numeric::concepts {

namespace implementation {

template <typename T>
concept ClassType = std::is_class_v<T>;

template <typename T>
concept HasCallableOperator = requires { &T::operator(); };

template <typename Type>
struct FunctionArgs {
  using args = void;
};

template <typename Ret, typename... Ts>
struct FunctionArgs<Ret(Ts...)> {
  using args = std::tuple<Ts...>;
};

template <typename Ret, typename... Ts>
struct FunctionArgs<Ret (*)(Ts...)> {
  using args = std::tuple<Ts...>;
};

template <ClassType ClassType, typename Ret, typename... Ts>
struct FunctionArgs<Ret (ClassType::*)(Ts...)> {
  using args = std::tuple<Ts...>;
};

template <ClassType ClassType, typename Ret, typename... Ts>
struct FunctionArgs<auto(ClassType::*)(Ts...) const->Ret> {
  using args = std::tuple<Ts...>;
};

template <typename Functor>
  requires HasCallableOperator<Functor>
struct FunctionArgs<Functor> {
  using args = typename FunctionArgs<decltype(&Functor::operator())>::args;
};

template <typename Function>
using FunctionArgsT = typename FunctionArgs<std::remove_cvref_t<Function>>::args;

template <typename T>
static constexpr auto kFunctionArgsCount = std::tuple_size_v<FunctionArgsT<std::remove_cvref_t<T>>>;

template <typename T>
static constexpr auto kArraySizeFromCallable =
    std::tuple_size_v<std::remove_cvref_t<std::tuple_element_t<0, typename FunctionArgs<T>::args>>>;

template <std::size_t N, typename T>
using NthFunctionParam = std::remove_cvref_t<std::tuple_element_t<N, typename FunctionArgs<T>::args>>;

template <typename T>
constexpr bool kIsTuple = false;

template <typename... Types>
constexpr bool kIsTuple<std::tuple<Types...>> = true;

}  // namespace implementation

template <typename T>
concept Tuple = implementation::kIsTuple<std::remove_cvref_t<T>>;

template <typename T>
concept FloatingPoint = std::floating_point<std::remove_cvref_t<T>>;

namespace implementation {

template <typename T>
constexpr bool kFloatingPointParamPack = false;

template <FloatingPoint... Types>
constexpr bool kFloatingPointParamPack<std::tuple<Types...>> = true;
}  // namespace implementation

template <typename T>
concept AllFloatingPointTup = Tuple<T> && implementation::kFloatingPointParamPack<std::remove_cvref_t<T>>;

template <typename... Args>
concept AllFloatingPoint = (FloatingPoint<Args> && ...);

template <typename Function, std::size_t N>
concept ScalarField =
    AllFloatingPointTup<implementation::FunctionArgsT<Function>> && implementation::kFunctionArgsCount<Function> ==
N;

template <typename T>
concept R1RealFunction = concepts::ScalarField<T, 1>;

template <typename T>
concept Integral = requires(T integral) {
                     { integral.low_ } -> FloatingPoint;
                     { integral.high_ } -> FloatingPoint;
                     { integral.function_ } -> R1RealFunction;
                   };

template <typename Iter, typename Type = typename Iter::value_type>
concept ReadOnlyRange = std::ranges::input_range<Iter>;

}  // namespace jr_numeric::concepts
