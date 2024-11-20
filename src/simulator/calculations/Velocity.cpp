#include "Velocity.h"
#include "utils/ArrayUtils.h"
#include <cmath>
#include "../particle/ParticleContainer.h"   

void Velocity::run(ParticleContainer& particles, double time_delta) {
  for (auto &p : particles) {
    auto v = p.getV();
    auto old_f = p.getOldF();
    auto new_f = p.getF();
    double m = p.getM();

    // Velocity-Störmer-Verlet formula (9)
    v = v + time_delta * (old_f + new_f) / (2 * m);

    p.updateV(v);
  }
}