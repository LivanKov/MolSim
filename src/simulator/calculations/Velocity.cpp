#include "Velocity.h"
#include "utils/ArrayUtils.h"
#include <cmath>

void Velocity::run(LinkedCellContainer &particles, double time_delta) {
  for (auto &p : particles)
    p.updateV(p.getV() +
              time_delta * (p.getOldF() + p.getF()) / (2 * p.getM()));
}