#include "../particle/container/LinkedCellContainer.h"
#include "Calculation.h"

/**
 * @struct Position
 * @brief Struct, that provides functions for position calculation.
 */
struct Position : AbstractPolicy {

  /**
   * @brief Executes the position calculation.
   * @param particles LinkedCellContainer reference.
   * @param time_delta Time step.
   * @param option Calculation option (e.g., DIRECT_SUM or LINKED_CELL).
   */
  static void run(LinkedCellContainer &particles, double time_delta,
                  OPTIONS option);
};