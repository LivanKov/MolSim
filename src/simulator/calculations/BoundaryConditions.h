#include "../particle/container/LinkedCellContainer.h"
#include "Calculation.h"

#pragma once
/**
 * @struct Position
 * @brief Struct, that provides functions for position calculation.
 */
struct BoundaryConditions : AbstractPolicy {
  /**
   * @brief Runs the boundary conditions for particles.
   * @param particles DirectSumContainer reference.
   */
  static void run(LinkedCellContainer &particles);

private:
  /**
   * @brief Handles the reflecting boundary conditions for particles.
   * @param particle_id Particle ID.
   * @param cell_index Cell index.
   * @param particles LinkedCellContainer reference.
   */
  static void handle_reflect_conditions(int particle_id, int cell_index,
                                        LinkedCellContainer &particles);
   /**
   * @brief Handles the periodic boundary conditions for particles.
   * @param particle_id Particle ID.
   * @param cell_index Cell index.
   * @param particles LinkedCellContainer reference.
   */
  static void handle_periodic_conditions(int particle_id, int cell_index,
                                         LinkedCellContainer &particles);

  /**
   * @brief Handles the outflow boundary conditions for particles.
   * @param particle_id Particle ID.
   * @param cell_index Cell index.
   * @param particles LinkedCellContainer reference.
   */
  static void handle_outflow_conditions(int particle_id, int cell_index,
                                        LinkedCellContainer &particles);
};