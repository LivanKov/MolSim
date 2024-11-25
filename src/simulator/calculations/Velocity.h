#include "../particle/ParticleContainer.h"
#include "Calculation.h"

struct Velocity : AbstractPolicy {
  static void run(ParticleContainer &particles, double time_delta);
};