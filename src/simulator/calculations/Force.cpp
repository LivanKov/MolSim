#include "Force.h"
#include "../particle/container/ParticleContainer.h"
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
  case MEMBRANE:
    membrane(particles);
    break;
  }
}

void Force::lennard_jones(LinkedCellContainer &particles, OPTIONS OPTION) {
  for (auto &particle : particles.particles) {
    particle.updateOldF(particle.getF());
    particle.updateF(0.0, 0.0, 0.0);
  }

  if (OPTION == OPTIONS::DIRECT_SUM) {
    compute_direct_sum(particles);
  } else {
    if (!SimParams::enable_omp) {
      compute_serial(particles);
    } else {
      if (SimParams::ompstrategy == OMPSTRATEGY::FORK_JOIN) {
        compute_parallel_fork_join(particles);
      } else {
        compute_parallel_tasking(particles);
      }
    }
    compute_ghost_cell_forces(particles);

    if (SimParams::enable_gravity) {
      apply_gravity(particles);
    }
  }
}

void Force::gravitational(LinkedCellContainer &particles, OPTIONS OPTION) {
  for (auto &particle : particles.particles) {
    particle.updateOldF(particle.getF());
    particle.updateF(0.0, 0.0, 0.0);
  }

  if (OPTION == OPTIONS::DIRECT_SUM) {
    compute_direct_sum(particles);
  } else {
    compute_serial(particles);
  }
}

// ----------------- Helper Functions -----------------------------------

void Force::compute_direct_sum(LinkedCellContainer &particles) {
  for (auto it = particles.particles.pair_begin();
       it != particles.particles.pair_end(); ++it) {
    auto r12 = it->second->getX() - it->first->getX();
    double distance = ArrayUtils::L2Norm(r12);

    if (distance > 1e-5) {
      auto force =
          compute_lj_force(it->first.get(), it->second.get(), r12, distance);
      it->first->updateF(it->first->getF() + force);
      it->second->updateF(it->second->getF() - force);
    }
  }
}

void Force::compute_serial(LinkedCellContainer &particles) {
  for (auto &particle : particles.particles) {
    for (auto &neighbour : particles.get_neighbours(particle.getId())) {
      if (*neighbour != particle) {
        auto r12 = neighbour->getX() - particle.getX();
        double distance = ArrayUtils::L2Norm(r12);

        if (distance > 1e-5) {
          auto force =
              compute_lj_force(&particle, neighbour.get(), r12, distance);
          particle.updateF(particle.getF() + force);
          if (!neighbour->is_fixed()) {
            neighbour->updateF(neighbour->getF() - force);
          }
        }
      }
    }
  }
}

void Force::compute_parallel_fork_join(LinkedCellContainer &particles) {

  std::vector<std::vector<std::array<double, 3>>> thread_local_forces(
      omp_get_max_threads());
  for (auto &forces : thread_local_forces) {
    forces.resize(particles.particles.size(), {0.0, 0.0, 0.0});
  }

#pragma omp parallel for
  for (size_t i = 0; i < particles.particles.size(); ++i) {
    int thread_id = omp_get_thread_num();
    auto &particle = particles.particles[i];
    for (auto &neighbour : particles.get_neighbours(particle.getId())) {
      if (*neighbour != particle) {
        auto r12 = neighbour->getX() - particle.getX();
        double distance = ArrayUtils::L2Norm(r12);

        if (distance > 1e-5) {
          auto force =
              compute_lj_force(&particle, neighbour.get(), r12, distance);

          thread_local_forces[thread_id][particle.getId()] =
              thread_local_forces[thread_id][particle.getId()] + force;
          thread_local_forces[thread_id][neighbour->getId()] =
              thread_local_forces[thread_id][neighbour->getId()] - force;
        }
      }
    }
  }

  for (size_t i = 0; i < particles.particles.size(); ++i) {
    auto &particle = particles.particles[i];
    for (int thread_id = 0; thread_id < omp_get_max_threads(); ++thread_id) {
      particle.updateF(particle.getF() + thread_local_forces[thread_id][i]);
    }
  }
}

void Force::compute_parallel_tasking(LinkedCellContainer &particles) {
  std::vector<std::array<double, 3>> local_forces(particles.particles.size(),
                                                  {0.0, 0.0, 0.0});

#pragma omp parallel
#pragma omp single
  {
    for (size_t i = 0; i < particles.particles.size(); ++i) {
      auto &particle = particles.particles[i];

#pragma omp task firstprivate(i) shared(particles, local_forces)
      {
        for (auto &neighbour : particles.get_neighbours(particle.getId())) {
          if (*neighbour != particle) {
            auto r12 = neighbour->getX() - particle.getX();
            double distance = ArrayUtils::L2Norm(r12);

            if (distance > 1e-5) {
              auto force =
                  compute_lj_force(&particle, neighbour.get(), r12, distance);
              local_forces[i] = local_forces[i] + force;
              if (!neighbour->is_fixed()) {
                local_forces[neighbour->getId()] =
                    local_forces[neighbour->getId()] - force;
              }
            }
          }
        }
      }
#pragma omp taskwait
    }

#pragma omp parallel for
    for (size_t i = 0; i < particles.particles.size(); ++i) {
      particles.particles[i].updateF(particles.particles[i].getF() +
                                     local_forces[i]);
    }
  }
}

void Force::compute_ghost_cell_forces(LinkedCellContainer &particles) {
  for (auto &particle : particles.particles) {
    for (auto &neighbour :
         particles.get_periodic_neighbours(particle.getId())) {
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
          particle.updateF(particle.getF() + force);
          if (!neighbour.ptr->is_fixed()) {
            neighbour.ptr->updateF(neighbour.ptr->getF() - force);
          }
        }
      }
    }
  }
}

std::array<double, 3> Force::compute_lj_force(Particle *p1, Particle *p2,
                                              const std::array<double, 3> &r12,
                                              double distance) {
  double sigma = (p1->getSigma() + p2->getSigma()) / 2;
  double epsilon = sqrt(p1->getEpsilon() * p2->getEpsilon());
  double term = sigma / distance;
  double term6 = pow(term, 6);
  double term12 = pow(term, 12);

  double totalForce = 24 * epsilon * (term6 - 2 * term12) / distance;
  return (totalForce / distance) * r12;
}

void Force::apply_gravity(LinkedCellContainer &particles) {
  for (auto &particle : particles.particles) {
    double gravitational_force_y = particle.getM() * SimParams::gravity;
    particle.updateF(particle.getF()[0],
                     particle.getF()[1] + gravitational_force_y,
                     particle.getF()[2]);
  }
}

void Force::membrane(LinkedCellContainer &particles) {
  for (auto &p : particles) {
    p.updateOldF(p.getF());
    p.updateF(0, 0, 0);
  }

  for (auto &p : particles) {
    for (auto neighbour : p.membrane_neighbours) {
      auto r12 = p.getX() - neighbour->getX();
      auto r21 = neighbour->getX() - p.getX();
      double distance = ArrayUtils::L2Norm(r12);
      auto totalForce = SimParams::membrane_stiffness *
                        (distance - SimParams::membrane_bond_length) *
                        (r21 / distance);
      p.updateF(p.getF() + totalForce);
      neighbour->updateF(neighbour->getF() - totalForce);
      double _sigma = (p.getSigma() + neighbour->getSigma()) / 2;
      double min_distance = MAGIC_NUMBER * _sigma;

      // put this in a separate function later
      if (distance < min_distance) {
        double term = _sigma / distance;
        double term6 = pow(term, 6);
        double term12 = pow(term, 12);
        double _epsilon;
        if (p.getEpsilon() == neighbour->getEpsilon()) {
          _epsilon = p.getEpsilon();
        } else {
          _epsilon = sqrt(p.getEpsilon() * neighbour->getEpsilon());
        }
        double lj_force = 24 * _epsilon * (term6 - 2 * term12) / distance;
        auto force = (lj_force / distance) * r12;
        p.updateF(p.getF() + force);
        neighbour->updateF(neighbour->getF() - force);
      }
    }

    for (auto neighbour : p.diagonal_membrane_neighbours) {
      auto r12 = p.getX() - neighbour->getX();
      auto r21 = neighbour->getX() - p.getX();
      double distance = ArrayUtils::L2Norm(r12);
      auto totalForce =
          SimParams::membrane_stiffness *
          (distance - SimParams::membrane_bond_length * TWO_SQRT) *
          (r21 / distance);
      p.updateF(p.getF() + totalForce);
      neighbour->updateF(neighbour->getF() - totalForce);
      double _sigma = p.getSigma() + neighbour->getSigma();
      double min_distance = MAGIC_NUMBER * _sigma;

      // put this in a separate function later
      if (distance < min_distance) {
        double term = _sigma / distance;
        double term6 = pow(term, 6);
        double term12 = pow(term, 12);
        double _epsilon;
        if (p.getEpsilon() == neighbour->getEpsilon()) {
          _epsilon = p.getEpsilon();
        } else {
          _epsilon = sqrt(p.getEpsilon() * neighbour->getEpsilon());
        }
        double lj_force = 24 * _epsilon * (term6 - 2 * term12) / distance;
        auto force = (lj_force / distance) * r12;
        p.updateF(p.getF() + force);
        neighbour->updateF(neighbour->getF() - force);
      }
    }

    if (SimParams::enable_z_gravity) {
      double gravitational_force_z = p.getM() * SimParams::z_gravity;
      p.updateF(p.getF()[0], p.getF()[1], p.getF()[2] + gravitational_force_z);
    }

    if (p.isApplyFZup() && SimParams::apply_fzup) {
      p.updateF(p.getF()[0], p.getF()[1],
                p.getF()[2] + SimParams::additional_force_zup);
    }
  }
}
