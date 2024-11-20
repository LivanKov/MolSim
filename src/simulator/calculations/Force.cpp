#include "Force.h"
#include "../particle/ParticleContainer.h"
#include "utils/ArrayUtils.h"


void Force::run_lennard_jones(ParticleContainer &particles){
    // store the current force as the old force and reset current to 0
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

    // this needs to be removed
    if (distance > 1e-5) {
      double totalForce;
        double term = SIGMA / distance;
        double term6 = pow(term, 6);
        double term12 = pow(term, 12);
        totalForce = 24 * EPSILON * (term6 - 2 * term12) / distance;
      
      auto force = (totalForce / distance) * r12;

      p1.updateF(p1.getF() + force);
      // Newton's third law
      p2.updateF(p2.getF() - force);
    }
  }
}

void Force::run_verlet(ParticleContainer &particles) {
  // store the current force as the old force and reset current to 0
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

    // this needs to be removed
    if (distance > 1e-5) {
      double totalForce;
      totalForce = p1.getM() * p2.getM() / pow(distance, 2);
      auto force = (totalForce / distance) * r12;
      p1.updateF(p1.getF() + force);
      // Newton's third law
      p2.updateF(p2.getF() - force);
    }
  }
}