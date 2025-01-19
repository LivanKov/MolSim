#include <concepts>
#include <type_traits>

#pragma once

/**
 * @struct AbstractPolicy
 * @brief Abstract struct, enforces the policy pattern for various calculation
 * types.
 */

enum OPTIONS { DIRECT_SUM, LINKED_CELLS };

struct AbstractPolicy {
  virtual ~AbstractPolicy() = 0;
};

/**
 * @struct Calculation
 * @brief Dummy struct that enforces the policy pattern for various calculation
 * types.
 */
template <typename Policy> struct Calculation : Policy {};