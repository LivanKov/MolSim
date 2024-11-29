#include "Position.h"
#include "../particle/ParticleContainer.h"
#include "utils/ArrayUtils.h"
#include <cmath>

void Position::run(ParticleContainer &particles, double time_delta) {
  for (auto &p : particles) {
    auto x = p.getX();
    auto v = p.getV();
    auto f = p.getF();
    double m = p.getM();

    // Velocity-Störmer-Verlet formula (8)
    x = x + time_delta * v + pow(time_delta, 2) * f / (2 * m);

    p.updateX(x);
  }
}