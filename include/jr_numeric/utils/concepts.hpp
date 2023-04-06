#pragma once

#include <concepts>
#include <tuple>
#include <type_traits>

namespace concepts {

namespace implementation {

template <typename T>
concept ClassType = std::is_class_v<T>;

template <typename T>
concept HasCallableOperator = requires { &T::operator(); };

template <typename Type>
struct FunctionArgs : std::false_type {};

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

template <typename Type>
using FunctionArgsT = typename FunctionArgs<Type>::args;

template <typename T>
static constexpr auto kFunctionArgsCount = std::tuple_size_v<typename FunctionArgs<T>::args>;

template <typename T>
static constexpr auto kArraySizeFromCallable =
    std::tuple_size_v<std::remove_cvref_t<std::tuple_element_t<0, typename FunctionArgs<T>::args>>>;

template <typename T>
using FirstFunctionParam = std::remove_cvref_t<std::tuple_element_t<0, typename FunctionArgs<T>::args>>;

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

}  // namespace concepts
