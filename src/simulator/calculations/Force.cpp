#include "Force.h"
#include "../particle/container/DirectSumContainer.h"
#include "io/input/cli/SimParams.h"
#include "utils/ArrayUtils.h"
#include <iostream>
#include <omp.h>

void Force::run(LinkedCellContainer &particles, ForceType type,
                OPTIONS OPTION) {
  switch (type) {
  case LENNARD_JONES:
    lennard_jones(particles, OPTION);
    break;
  case GRAVITATIONAL:
    gravitational(particles, OPTION);
    break;
  }
}

void Force::lennard_jones(LinkedCellContainer &particles, OPTIONS OPTION) {

  for (auto &p : particles.particles) {
    p.updateOldF(p.getF());
    p.updateF(0, 0, 0);
  }

  if (OPTION == OPTIONS::DIRECT_SUM) {
    for (auto it = particles.particles.pair_begin();
         it != particles.particles.pair_end(); ++it) {

      auto r12 = it->second->getX() - it->first->getX();
      double distance = ArrayUtils::L2Norm(r12);

      if (distance > 1e-5) {
        double sigma = (it->first->getSigma() + it->second->getSigma()) / 2;
        double term = sigma / distance;
        double term6 = pow(term, 6);
        double term12 = pow(term, 12);
        double epsilon =
            sqrt(it->first->getEpsilon() * it->second->getEpsilon());

        double totalForce = 24 * epsilon * (term6 - 2 * term12) / distance;
        auto force = (totalForce / distance) * r12;
        it->first->updateF(it->first->getF() + force);
        it->second->updateF(it->second->getF() - force);
      }
    }
  } else {

    if (!SimParams::enable_omp) {
      for (size_t i = 0; i < particles.particles.size(); ++i) {
        auto &particle = particles.particles[i];
        for (auto &neighbour : particles.get_neighbours(particle.getType())) {
          if (*neighbour != particle) {
            auto r12 = neighbour->getX() - particle.getX();
            double distance = ArrayUtils::L2Norm(r12);

            if (distance > 1e-5) {
              double sigma = (particle.getSigma() + neighbour->getSigma()) / 2;
              double term = sigma / distance;
              double term6 = pow(term, 6);
              double term12 = pow(term, 12);
              double epsilon =
                  sqrt(particle.getEpsilon() * neighbour->getEpsilon());

              double totalForce =
                  24 * epsilon * (term6 - 2 * term12) / distance;
              auto force = (totalForce / distance) * r12;
              particle.updateF(particle.getF() + force);
              neighbour->updateF(neighbour->getF() - force);
            }
          }
        }
      }
    } else {
      if (SimParams::ompstrategy == OMPSTRATEGY::FORK_JOIN) {
        std::vector<std::array<double, 3>> thread_local_forces(
            particles.particles.size(), {0.0, 0.0, 0.0});

#pragma omp parallel default(none) shared(particles, thread_local_forces)
        {
#pragma omp for
          for (size_t i = 0; i < particles.particles.size(); ++i) {
            auto &particle = particles.particles[i];
            for (auto &neighbour :
                 particles.get_neighbours(particle.getType())) {
              if (*neighbour != particle) {
                auto r12 = neighbour->getX() - particle.getX();
                double distance = ArrayUtils::L2Norm(r12);

                if (distance > 1e-5) {
                  double sigma =
                      (particle.getSigma() + neighbour->getSigma()) / 2;
                  double term = sigma / distance;
                  double term6 = pow(term, 6);
                  double term12 = pow(term, 12);
                  double epsilon =
                      sqrt(particle.getEpsilon() * neighbour->getEpsilon());

                  double totalForce =
                      24 * epsilon * (term6 - 2 * term12) / distance;
                  auto force = (totalForce / distance) * r12;

                  // Accumulate forces in thread-local storage
                  thread_local_forces[i] = thread_local_forces[i] + force;
                  thread_local_forces[neighbour->getType()] =
                      thread_local_forces[neighbour->getType()] - force;
                }
              }
            }
          }
// Combine thread-local forces into shared forces
#pragma omp for
          for (size_t i = 0; i < particles.particles.size(); ++i) {
            particles.particles[i].updateF(particles.particles[i].getF() +
                                           thread_local_forces[i]);
          }
        }
      } else {

        std::vector<std::array<double, 3>> local_forces(
            particles.particles.size(), {0.0, 0.0, 0.0});

#pragma omp parallel
#pragma omp single
        {
          for (size_t i = 0; i < particles.particles.size(); ++i) {
            auto &particle = particles.particles[i];

#pragma omp task firstprivate(i) shared(particles, local_forces)
            {
              for (auto &neighbour :
                   particles.get_neighbours(particle.getType())) {
                if (*neighbour != particle) {
                  auto r12 = neighbour->getX() - particle.getX();
                  double distance = ArrayUtils::L2Norm(r12);

                  if (distance > 1e-5) {
                    double sigma =
                        (particle.getSigma() + neighbour->getSigma()) / 2;
                    double term = sigma / distance;
                    double term6 = pow(term, 6);
                    double term12 = pow(term, 12);
                    double epsilon =
                        sqrt(particle.getEpsilon() * neighbour->getEpsilon());

                    double totalForce =
                        24 * epsilon * (term6 - 2 * term12) / distance;
                    auto force = (totalForce / distance) * r12;

                    // Accumulate forces locally
                    local_forces[i] = local_forces[i] + force;
                    local_forces[neighbour->getType()] =
                        local_forces[neighbour->getType()] - force;
                  }
                }
              }
            }
          }
#pragma omp taskwait
        }

// Update particles with accumulated forces
#pragma omp parallel for
        for (size_t i = 0; i < particles.particles.size(); ++i) {
          particles.particles[i].updateF(particles.particles[i].getF() +
                                         local_forces[i]);
        }
      }
    }

    for (size_t i = 0; i < particles.particles.size(); ++i) {
      auto &particle = particles.particles[i];
      std::array<double, 3> local_force = {0.0, 0.0, 0.0};

      for (auto &neighbour :
           particles.get_additional_neighbour_indices(particle.getType())) {
        if (*(neighbour.ptr) != particle) {
          auto r12 = neighbour.position - particle.getX();
          double distance = ArrayUtils::L2Norm(r12);

          if (distance > 1e-5) {
            double sigma = (particle.getSigma() + neighbour.sigma) / 2;
            double term = sigma / distance;
            double term6 = pow(term, 6);
            double term12 = pow(term, 12);
            double epsilon = sqrt(particle.getEpsilon() * neighbour.epsilon);

            double totalForce = 24 * epsilon * (term6 - 2 * term12) / distance;
            auto force = (totalForce / distance) * r12;
            local_force = local_force + force;

            neighbour.ptr->updateF(neighbour.ptr->getF() - force);
          }
        }
      }
      particle.updateF(particle.getF() + local_force);
    }
  }

  if (SimParams::enable_gravity) {
    for (auto &particle : particles.particles) {
      double gravitational_force_y = particle.getM() * SimParams::gravity;
      particle.updateF(particle.getF()[0],
                       particle.getF()[1] + gravitational_force_y,
                       particle.getF()[2]);
    }
  }
}

void Force::gravitational(LinkedCellContainer &particles, OPTIONS OPTION) {
  // store the current force as the old force and reset current to 0

  if (OPTION == OPTIONS::DIRECT_SUM) {
    for (auto &p : particles.particles) {
      p.updateOldF(p.getF());
      p.updateF(0, 0, 0);
    }

    // Iterate each pair
    for (auto it = particles.particles.pair_begin();
         it != particles.particles.pair_end(); ++it) {
      ParticlePair &pair = *it;
      Particle &p1 = *(pair.first);
      Particle &p2 = *(pair.second);
      auto r12 = p2.getX() - p1.getX();
      // distance ||x_i - x_j ||
      double distance = ArrayUtils::L2Norm(r12);

      if (distance > 1e-5) {
        double totalForce;
        totalForce = p1.getM() * p2.getM() / pow(distance, 2);
        auto force = (totalForce / distance) * r12;
        p1.updateF(p1.getF() + force);
        // Newton's third law
        p2.updateF(p2.getF() - force);
      }
    }
  } else {
    for (auto &p : particles.particles) {
      p.updateOldF(p.getF());
      p.updateF(0, 0, 0);
    }
    for (auto &particle : particles.particles) {
      for (auto neighbour : particles.get_neighbours(particle.getType())) {
        auto r12 = neighbour->getX() - particle.getX();
        double distance = ArrayUtils::L2Norm(r12);
        if (distance > 1e-5) {
          double totalForce;
          totalForce = particle.getM() * neighbour->getM() / pow(distance, 2);
          auto force = (totalForce / distance) * r12;
          particle.updateF(particle.getF() + force);
          neighbour->updateF(neighbour->getF() - force);
        }
      }
    }
  }
}