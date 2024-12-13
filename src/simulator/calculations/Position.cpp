#include "Position.h"
#include "utils/ArrayUtils.h"
#include <cmath>

void Position::run(LinkedCellContainer &particles, double time_delta,
                   OPTIONS option) {
  for (auto &p : particles.particles) {
    if (!p.left_domain) {
      auto old_x = p.getX();
      p.updateX(p.getX() + time_delta * p.getV() +
                pow(time_delta, 2) * p.getF() / (2 * p.getM()));
      if (option == OPTIONS::LINKED_CELLS) {
        particles.update_particle_location(p.getType(), old_x);
      }
    }
  }
}