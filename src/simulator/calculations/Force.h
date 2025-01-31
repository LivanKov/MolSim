#include "../particle/container/LinkedCellContainer.h"
#include "../particle/container/ParticleContainer.h"
#include "Calculation.h"

#pragma once
#define TWO_SQRT 1.41421356237
#define MAGIC_NUMBER 1.12246204831

/**
 * @enum ForceType
 * @brief Enum class for the force calculation type.
 */
enum ForceType { LENNARD_JONES, GRAVITATIONAL, MEMBRANE };

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

  // ----------------- Helper Functions -----------------------------------

  /**
   * @brief Computes the direct sum force between all particles in the
   * container.
   * @param particles LinkedCellContainer reference.
   */
  static void compute_direct_sum(LinkedCellContainer &particles);

  /**
   * @brief Computes the force between all particles in the container using
   * serial execution.
   * @param particles LinkedCellContainer reference.
   */
  static void compute_serial(LinkedCellContainer &particles);

  /**
   * @brief Computes the force between all particles in the container using
   * parallel execution. Utilizes the fork-join policy.
   * @param particles LinkedCellContainer reference.
   */
  static void compute_parallel_fork_join(LinkedCellContainer &particles);

  /**
   * @brief Computes the force between all particles in the container using
   * parallel execution using tasking. Utilizes the tasking policy.
   * @param particles LinkedCellContainer reference.
   */
  static void compute_parallel_tasking(LinkedCellContainer &particles);

  /**
   * @brief Computes forces affected by the ghost cells.
   * @param particles LinkedCellContainer reference.
   */
  static void compute_ghost_cell_forces(LinkedCellContainer &particles);

  /**
   * @brief Computes the Lennard-Jones force between two particles.
   * @param p1 First particle.
   * @param p2 Second particle.
   * @param r12 Vector representing the distance between the particles.
   * @param distance Distance between the particles.
   * @return The Lennard-Jones force between the particles.
   */
  static std::array<double, 3>
  compute_lj_force(Particle *p1, Particle *p2, const std::array<double, 3> &r12,
                   double distance);

  /**
   * @brief Applies gravity to all particles in the container.
   * @param particles LinkedCellContainer reference.
   */
  static void apply_gravity(LinkedCellContainer &particles);

  /**
   * @brief Computes the gravitational force between all particles in the
   * container.
   * @param particles LinkedCellContainer reference.
   * @param OPTION Gravitational force option.
   */
  static void gravitational(LinkedCellContainer &particles, OPTIONS OPTION);

  /**
   * @brief Computes the force within a connected membrane of particles in the
   * container.
   * @param particles LinkedCellContainer reference.
   */
  static void membrane(LinkedCellContainer &particles);
};