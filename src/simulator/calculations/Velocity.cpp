#include "Velocity.h"
#include "io/input/cli/SimParams.h"
#include "utils/ArrayUtils.h"
#include <cmath>

void Velocity::run(LinkedCellContainer &particles, double time_delta) {
  for (auto &p : particles.particles) {
    if (SimParams::enable_v_threshold) {
      auto v =
          p.getV() + time_delta * (p.getOldF() + p.getF()) / (2 * p.getM());
      p.updateV(v[0] > SimParams::v_threshold ? SimParams::v_threshold : v[0],
                v[1] > SimParams::v_threshold ? SimParams::v_threshold : v[1],
                v[2] > SimParams::v_threshold ? SimParams::v_threshold : v[2]);
    } else {
      p.updateV(p.getV() +
                time_delta * (p.getOldF() + p.getF()) / (2 * p.getM()));
    }
  }
}