#include "../particle/ParticleContainer.h"
#include "../particle/LinkedCellContainer.h"
#include "Calculation.h"

#define EPSILON 5.0
#define SIGMA 1.0

/**
 * @enum ForceType
 * @brief Enum class for the force calculation type.
 */
enum ForceType { LENNARD_JONES, VERLET };

/***
 * @struct Force
 * @brief Struct, that provides functions for force calculation.
 **/
struct Force : AbstractPolicy {
  static void run(LinkedCellContainer &particles, ForceType type, OPTIONS OPTION);

private:
  /***
   * @brief Lennard-Jones force calculation.
   * @param particles ParticleContainer reference.
   **/
  static void lennard_jones(LinkedCellContainer &particles, OPTIONS OPTION);
  /**
   * @brief Verlet force calculation.
   * @param particles ParticleContainer reference.
   */
  static void verlet(LinkedCellContainer &particles, OPTIONS OPTION);
};