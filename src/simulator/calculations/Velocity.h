#include "../particle/container/LinkedCellContainer.h"
#include "Calculation.h"

/**
 * @struct Position
 * @brief Struct, that provides functions for position calculation.
 */
struct Velocity : AbstractPolicy {

  /**
   * @brief Executes the velocity calculation.
   * @param particles LinkedCellContainer reference.
   * @param time_delta Time step.
   */
  static void run(LinkedCellContainer &particles, double time_delta);
};