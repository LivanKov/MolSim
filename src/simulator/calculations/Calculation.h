#include <concepts>
#include <type_traits>

#pragma once

/**
 * @struct AbstractPolicy
 * @brief Abstract struct, enforces the policy pattern for various calculation
 * types.
 */

enum OPTIONS { NONE, LINKED_CELLS_BOUNDARY };

struct AbstractPolicy {
  virtual ~AbstractPolicy() = 0;
};

template <typename T>
concept CalculationPolicy = std::is_base_of_v<AbstractPolicy, T>;

/**
 * @struct Calculation
 * @brief Dummy struct that enforces the policy pattern for various calculation
 * types.
 */
template <CalculationPolicy Policy> struct Calculation : Policy {};