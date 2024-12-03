#include "Force.h"
#include "../particle/ParticleContainer.h"
#include "utils/ArrayUtils.h"

void Force::run(LinkedCellContainer &particles, ForceType type,
                OPTIONS OPTION) {
  switch (type) {
  case LENNARD_JONES:
    lennard_jones(particles, OPTION);
    break;
  case VERLET:
    verlet(particles, OPTION);
    break;
  }
}

void Force::lennard_jones(LinkedCellContainer &particles, OPTIONS OPTION) {
  if (OPTION == OPTIONS::NONE) {

    for (auto &p : particles) {
      p.updateOldF(p.getF());
      p.updateF(0, 0, 0);
    }

    for (auto it = particles.pair_begin(); it != particles.pair_end(); ++it) {
      auto r12 = it->second->getX() - it->first->getX();
      double distance = ArrayUtils::L2Norm(r12);

      double totalForce;
      double term = SIGMA / distance;
      double term6 = pow(term, 6);
      double term12 = pow(term, 12);
      totalForce = 24 * EPSILON * (term6 - 2 * term12) / distance;

      auto force = (totalForce / distance) * r12;

      it->first->updateF(it->first->getF() + force);
      it->second->updateF(it->second->getF() - force);
    }
  } else {
    for (auto &p : particles) {
      p.updateOldF(p.getF());
      p.updateF(0, 0, 0);
    }

    for (auto &particle : particles) {
      for (auto neighbour : particles.get_neighbours(particle)) {
        auto r12 = neighbour->getX() - particle.getX();
        double distance = ArrayUtils::L2Norm(r12);

        double totalForce;
        double term = SIGMA / distance;
        double term6 = pow(term, 6);
        double term12 = pow(term, 12);
        totalForce = 24 * EPSILON * (term6 - 2 * term12) / distance;

        auto force = (totalForce / distance) * r12;

        particle.updateF(particle.getF() + force);
        neighbour->updateF(neighbour->getF() - force);
      }
    }
  }
}

void Force::verlet(LinkedCellContainer &particles, OPTIONS OPTION) {
  // store the current force as the old force and reset current to 0

  if (OPTION == OPTIONS::NONE) {
    for (auto &p : particles) {
      p.updateOldF(p.getF());
      p.updateF(0, 0, 0);
    }

    // Iterate each pair
    for (auto it = particles.pair_begin(); it != particles.pair_end(); ++it) {
      ParticlePair &pair = *it;
      Particle &p1 = *(pair.first);
      Particle &p2 = *(pair.second);
      auto r12 = p2.getX() - p1.getX();
      // distance ||x_i - x_j ||
      double distance = ArrayUtils::L2Norm(r12);

      double totalForce;
      totalForce = p1.getM() * p2.getM() / pow(distance, 2);
      auto force = (totalForce / distance) * r12;
      p1.updateF(p1.getF() + force);
      // Newton's third law
      p2.updateF(p2.getF() - force);
    }
  } else {
    for (auto &p : particles) {
      p.updateOldF(p.getF());
      p.updateF(0, 0, 0);
    }

    for (auto &particle : particles) {
      for (auto neighbour : particles.get_neighbours(particle)) {
        auto r12 = neighbour->getX() - particle.getX();
        double distance = ArrayUtils::L2Norm(r12);

        double totalForce;
        totalForce = particle.getM() * neighbour->getM() / pow(distance, 2);
        auto force = (totalForce / distance) * r12;
        particle.updateF(particle.getF() + force);
        neighbour->updateF(neighbour->getF() - force);
      }
    }
  }
}