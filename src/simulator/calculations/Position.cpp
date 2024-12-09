#include "Position.h"
#include "utils/ArrayUtils.h"
#include <cmath>

void Position::run(LinkedCellContainer &particles, double time_delta) {
  for (auto &p : particles.particles){
    auto old_x = p.getX();
    p.updateX(p.getX() + time_delta * p.getV() +
              pow(time_delta, 2) * p.getF() / (2 * p.getM()));
    particles.update_particle_location(std::make_shared<Particle>(p), old_x);
  }
}