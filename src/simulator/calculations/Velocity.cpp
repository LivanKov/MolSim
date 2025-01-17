#include "Velocity.h"
#include "io/input/cli/SimParams.h"
#include "utils/ArrayUtils.h"
#include <cmath>

void Velocity::run(LinkedCellContainer &particles, double time_delta) {
  for (auto &p : particles.particles) {
    auto v = p.getV() + time_delta * (p.getOldF() + p.getF()) / (2 * p.getM());
    p.updateV(v);
  }
}