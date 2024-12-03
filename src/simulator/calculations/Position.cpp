#include "Position.h"
#include "../particle/ParticleContainer.h"
#include "utils/ArrayUtils.h"
#include <cmath>

void Position::run(ParticleContainer &particles, double time_delta) {
  for (auto &p : particles)
    p.updateX(p.getX() + time_delta * p.getV() +
              pow(time_delta, 2) * p.getF() / (2 * p.getM()));
}