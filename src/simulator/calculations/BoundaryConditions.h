#include "../particle/container/LinkedCellContainer.h"
#include "Calculation.h"

#pragma once
/**
 * @struct Position
 * @brief Struct, that provides functions for position calculation.
 */
struct BoundaryConditions : AbstractPolicy {
  static void run(LinkedCellContainer &particles);

private:
  static void handle_reflect_conditions(int particle_id, int cell_index,
                                        LinkedCellContainer &particles);
  static void handle_periodic_conditions(LinkedCellContainer &particles);
  static void handle_outflow_conditions(int particle_id, int cell_index,
                                        LinkedCellContainer &particles);
};