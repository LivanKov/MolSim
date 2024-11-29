#include "../particle/ParticleContainer.h"
#include "Calculation.h"

/**
 * @struct Position
 * @brief Struct, that provides functions for position calculation.
 */
struct Position : AbstractPolicy {
  static void run(ParticleContainer &particles, double time_delta);
};