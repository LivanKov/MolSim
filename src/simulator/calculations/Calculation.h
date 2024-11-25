#include <concepts>
#include <type_traits>

#pragma once

struct AbstractPolicy {
  virtual ~AbstractPolicy() = 0;
};

template <typename T>
concept CalculationPolicy = std::is_base_of_v<AbstractPolicy, T>;

template <CalculationPolicy Policy> struct Calculation : Policy {};