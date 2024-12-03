#include "../particle/LinkedCellContainer.h"
#include "Calculation.h"

/**
 * @struct Position
 * @brief Struct, that provides functions for position calculation.
 */
struct Position : AbstractPolicy {
  static void run(LinkedCellContainer &particles, double time_delta);
};