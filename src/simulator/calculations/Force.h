#include "../particle/container/LinkedCellContainer.h"
#include "../particle/container/ParticleContainer.h"
#include "Calculation.h"

#define TWO_SQRT 1.41421356237
#define MAGIC_NUMBER 1.12246204831 

/**
 * @enum ForceType
 * @brief Enum class for the force calculation type.
 */
enum ForceType { LENNARD_JONES, GRAVITATIONAL, MEMBRANE };

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
   * @param particles ParticleContainer reference.
   **/
  static void lennard_jones(LinkedCellContainer &particles, OPTIONS OPTION);
  /**
   * @brief Gravitational force calculation.
   * @param particles ParticleContainer reference.
   */
  static void gravitational(LinkedCellContainer &particles, OPTIONS OPTION);

  static void membrane(LinkedCellContainer &particles);
};