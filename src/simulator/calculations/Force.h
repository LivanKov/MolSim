#include "../particle/ParticleContainer.h"
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
  static void run(ParticleContainer &particles, ForceType type);

private:
  /*** 
   * @brief Lennard-Jones force calculation.
   * @param particles ParticleContainer reference.
  **/
  static void lennard_jones(ParticleContainer &particles);
  /**
   * @brief Verlet force calculation.
   * @param particles ParticleContainer reference.
  */
  static void verlet(ParticleContainer &particles);
};