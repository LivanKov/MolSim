#include <concepts>
#include <type_traits>

#pragma once

/**
 * @enum OPTIONS
 * @brief Enum class for calculation options.
 */
enum OPTIONS { DIRECT_SUM, LINKED_CELLS };

/**
 * @struct AbstractPolicy
 * @brief Abstract struct, enforces the policy pattern for various calculation
 * types.
 */

struct AbstractPolicy {
  /**
   * @brief Virtual destructor. Abstract class.
   */
  virtual ~AbstractPolicy() = 0;
};

/**
 * @struct Calculation
 * @brief Dummy struct that enforces the policy pattern for various calculation
 * types.
 */
template <typename Policy> struct Calculation : Policy {};