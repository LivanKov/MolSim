#include "../particle/ParticleContainer.h"
#include "Calculation.h"

struct Position : AbstractPolicy {
  static void run(ParticleContainer &particles, double time_delta);
};