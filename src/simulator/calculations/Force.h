#include "../particle/container/DirectSumContainer.h"
#include "../particle/container/LinkedCellContainer.h"
#include "Calculation.h"

#pragma once

/**
 * @enum ForceType
 * @brief Enum class for the force calculation type.
 */
enum ForceType { LENNARD_JONES, GRAVITATIONAL };

enum OMPSTRATEGY {FORK_JOIN, TASKING};

inline std::string to_string(OMPSTRATEGY strategy) {
  switch (strategy) {
  case OMPSTRATEGY::FORK_JOIN:
    return "FORK_JOIN";
  case OMPSTRATEGY::TASKING:
    return "TASKING";
  default:
    return "UNKNOWN";
  }
}

/***
 * @struct Force
 * @brief Struct, that provides functions for force calculation.
 **/
struct Force : AbstractPolicy {
  static void run(LinkedCellContainer &particles, ForceType type,
                  OPTIONS OPTION);

private:
  /***
   * @brief Lennard-Jones force calculation.
   * @param particles DirectSumContainer reference.
   **/
  static void lennard_jones(LinkedCellContainer &particles, OPTIONS OPTION);
  /**
   * @brief Gravitational force calculation.
   * @param particles DirectSumContainer reference.
   */
  static void gravitational(LinkedCellContainer &particles, OPTIONS OPTION);
};