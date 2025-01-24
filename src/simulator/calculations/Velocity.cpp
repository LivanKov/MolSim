#include "Velocity.h"
#include "io/input/cli/SimParams.h"
#include "utils/ArrayUtils.h"
#include <cmath>

void Velocity::run(LinkedCellContainer &particles, double time_delta) {
  for (auto &p : particles.particles) {
    if(!p.fixed) {
    auto v = p.getV() + time_delta * (p.getOldF() + p.getF()) / (2 * p.getM());

    if (SimParams::enable_v_threshold) {
      // Clamp velocity based on absolute threshold (positive and negative)
      for (int i = 0; i < 3; ++i) {
        if (v[i] > SimParams::v_threshold) {
          v[i] = SimParams::v_threshold; // Positive threshold
        } else if (v[i] < -SimParams::v_threshold) {
          v[i] = -SimParams::v_threshold; // Negative threshold
        }
      }
    }

    p.updateV(v);
  }
}
}