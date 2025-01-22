#include "../particle/container/DirectSumContainer.h"
#include "../particle/container/LinkedCellContainer.h"
#include "Calculation.h"

#pragma once

/**
 * @enum ForceType
 * @brief Enum class for the force calculation type.
 */
enum ForceType { LENNARD_JONES, GRAVITATIONAL };

/**
 * @enum OMPSTRATEGY
 * @brief Enum class for OpenMP parallelization strategies.
 */
enum OMPSTRATEGY { FORK_JOIN, TASKING };

/**
 * @brief Converts OMPSTRATEGY enum to string.
 * @param strategy The OMPSTRATEGY value.
 * @return String representation of the strategy.
 */
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
  /**
   * @brief Executes the selected force calculation.
   * @param particles LinkedCellContainer reference.
   * @param type Type of force (e.g., LENNARD_JONES or GRAVITATIONAL).
   * @param OPTION Calculation option (e.g., DIRECT_SUM or LINKED_CELL).
   */
  static void run(LinkedCellContainer &particles, ForceType type,
                  OPTIONS OPTION);

private:
  /***
   * @brief Lennard-Jones force calculation.
   * @param particles LinkedCellContainer reference.
   * @param OPTION Calculation option (e.g., DIRECT_SUM or LINKED_CELL).
   **/
  static void lennard_jones(LinkedCellContainer &particles, OPTIONS OPTION);
  /**
   * @brief Gravitational force calculation.
   * @param particles LinkedCellContainer reference.
   * @param OPTION Calculation option (e.g., DIRECT_SUM or LINKED_CELL).
   */
  static void gravitational(LinkedCellContainer &particles, OPTIONS OPTION);

  // ----------------- Helper Functions -----------------------------------

  static void compute_direct_sum(LinkedCellContainer &particles);
  static void compute_serial(LinkedCellContainer &particles);
  static void compute_parallel_fork_join(LinkedCellContainer &particles);
  static void compute_parallel_tasking(LinkedCellContainer &particles);
  static void compute_ghost_cell_forces(LinkedCellContainer &particles);
  static std::array<double, 3>
  compute_lj_force(Particle *p1, Particle *p2, const std::array<double, 3> &r12,
                   double distance);
  static void apply_gravity(LinkedCellContainer &particles);
};