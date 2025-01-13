/*
 * Created by sebastianpse on 11/9/24.
 */

#include "ParticleGenerator.h"
#include "../Simulation.h"
#include "container/LinkedCellContainer.h"
#include "io/input/cli/SimParams.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include <cmath>
#include <random>

#include "utils/logger/Logger.h"

ParticleGenerator::ParticleGenerator() = default;

void ParticleGenerator::insertCuboid(
    const std::array<double, 3> &lowerLeftFrontCorner,
    const std::array<size_t, 3> &dimensions, double h, double mass,
    const std::array<double, 3> &initialVelocity,
    LinkedCellContainer &particle_container, double epsilon, double sigma,
    bool is_membrane,
    std::vector<std::array<double, 3>> additional_force_coordinates) {
  for (size_t i = 0; i < dimensions[2]; ++i) {
    for (size_t j = 0; j < dimensions[1]; ++j) {
      for (size_t k = 0; k < dimensions[0]; ++k) {
        std::array position = {lowerLeftFrontCorner[0] + k * h,
                               lowerLeftFrontCorner[1] + j * h,
                               lowerLeftFrontCorner[2] + i * h};

        std::array<double, 3> velocity = initialVelocity;

        Particle particle(position, velocity, mass,
                          particle_container.particle_id, epsilon, sigma);

        if(std::find(additional_force_coordinates.begin(), additional_force_coordinates.end(), position) != additional_force_coordinates.end())
          SimParams::additional_force_particle_ids.insert(particle_container.particle_id);
        particle_container.particle_id++;
        Logger::getInstance().trace("New Particle generated");
        particle_container.insert(particle, true);
        Logger::getInstance().trace("New Particle inserted into container");
      }
    }
  }
  if (!SimParams::fixed_Domain) {
    particle_container.readjust();
  }
  if (is_membrane) {
    for (size_t i = 0; i < particle_container.size(); i++)
      generate_membrane(i, particle_container, dimensions);
  }
  Logger::getInstance().info("New cuboid generated");
}

void ParticleGenerator::insertDisc(const std::array<double, 3> &center,
                                   const std::array<double, 3> &initialVelocity,
                                   size_t radius, double h, double mass,
                                   LinkedCellContainer &particles,
                                   double epsilon, double sigma,
                                   bool is_membrane) {

  // start on the point with the leftest x point then go the rightest x point
  for (double x = center[0] - radius * h; x <= center[0] + radius * h; x += h) {
    // same for y with lowest and highest point
    for (double y = center[1] - radius * h; y <= center[1] + radius * h;
         y += h) {

      double d_x = x - center[0];
      double d_y = y - center[1];
      // checks if particle is in the permitted radius
      if (d_x * d_x + d_y * d_y <= radius * radius * h * h) {

        // set z-coordinate of new particle to z-coordinate of center particle
        // because we create disc along the XY-plane
        double z = center[2];

        std::array position = {x, y, z};
        std::array<double, 3> velocity = initialVelocity;

        Particle particle(position, velocity, mass, particles.particle_id,
                          epsilon, sigma);
        particles.particle_id++;
        Logger::getInstance().trace("New Particle generated");
        particles.insert(particle, true);
        Logger::getInstance().trace("New Particle inserted into container");
      }
    }
  }
  if (!SimParams::fixed_Domain) {
    particles.readjust();
  }
  if (is_membrane) {
    Logger::getInstance().info("Membrane generated");
  }
  Logger::getInstance().info("New disk generated");
}

void ParticleGenerator::insertSingleMolecule(
    const std::array<double, 3> &position,
    const std::array<double, 3> &velocity, double mass,
    LinkedCellContainer &particles) {
  Particle particle(position, velocity, mass, particles.particle_id);
  particles.particle_id++;
  Logger::getInstance().trace("New Particle generated");
  particles.insert(particle, true);
  Logger::getInstance().trace("New Particle inserted into container");
  if (!SimParams::fixed_Domain) {
    particles.readjust();
  }
  Logger::getInstance().info("New single molecule generated");
}

void ParticleGenerator::generate_membrane(
    size_t particle_index, LinkedCellContainer &particle_container,
    std::array<size_t, 3> dimensions) {
  size_t i = particle_index % dimensions[0];
  size_t j = particle_index / dimensions[0];
  size_t k = particle_index / (dimensions[0] * dimensions[1]);

  auto particle = particle_container.particles[particle_index];
  int particle_id = particle.getId();

  for (int di = -1; di <= 1; ++di) {
    for (int dj = -1; dj <= 1; ++dj) {
      for (int dk = -1; dk <= 1; ++dk) {
        if (di == 0 && dj == 0 && dk == 0) {
          continue;
        }
        int ni = i + di;
        int nj = j + dj;
        int nk = k + dk;

        if (ni >= 0 && ni < dimensions[0] && nj >= 0 && nj < dimensions[1] &&
            nk >= 0 && nk < dimensions[2]) {
          int membrane_neighbour_index =
              ni + (nj * dimensions[0]) + nk * dimensions[0] * dimensions[1];
          int neighbour_id =
              particle_container.particles[membrane_neighbour_index].getId();

          // diagonal membrane member
          if (std::abs(di) + std::abs(dj) + std::abs(dk) > 1) {
            particle_container.cells_map[particle_id]
                ->diagonal_membrane_neighbours.push_back(
                    particle_container.cells_map[neighbour_id]);
          } else {
            // horizontal membrane member
            particle_container.cells_map[particle_id]
                ->membrane_neighbours.push_back(
                    particle_container.cells_map[neighbour_id]);
          }
        }
      }
    }
  }
}
