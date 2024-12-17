#include "../particle/container/LinkedCellContainer.h"
#include "Calculation.h"

#pragma once
/**
 * @struct Position
 * @brief Struct, that provides functions for position calculation.
 */
struct BoundaryConditions : AbstractPolicy {
  static void run(LinkedCellContainer &particles);
};