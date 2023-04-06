#pragma once

#include <concepts>
#include <type_traits>

#include "jr_numeric/utils/concepts.hpp"

namespace meta {

template <
    concepts::R1RealFunction FunctionType,
    concepts::FloatingPoint Param = concepts::implementation::NthFunctionParam<0, FunctionType>>
using RealFunctionResult = std::invoke_result_t<FunctionType, Param>;

template <
    concepts::Integral IntegralType,
    concepts::FloatingPoint Param = concepts::implementation::NthFunctionParam<0, decltype(IntegralType::function_)>>
using IntegralFunctionResult = RealFunctionResult<decltype(IntegralType::function_), Param>;

}  // namespace meta
