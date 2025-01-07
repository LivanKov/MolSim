/*
 * Created by sebastianpse on 11/9/24.
 */

#include "ParticleGenerator.h"
#include "../Simulation.h"
#include "container/LinkedCellContainer.h"
#include "io/input/cli/SimParams.h"
#include "utils/MaxwellBoltzmannDistribution.h"

#include <random>

#include "utils/logger/Logger.h"

ParticleGenerator::ParticleGenerator() = default;

void ParticleGenerator::insertCuboid(
    const std::array<double, 3> &lowerLeftFrontCorner,
    const std::array<size_t, 3> &dimensions, double h, double m,
    const std::array<double, 3> &initialVelocity,
    LinkedCellContainer &particles, double epsilon, double sigma) {
  for (size_t i = 0; i < dimensions[2]; ++i) {
    for (size_t j = 0; j < dimensions[1]; ++j) {
      for (size_t k = 0; k < dimensions[0]; ++k) {
        std::array position = {lowerLeftFrontCorner[0] + k * h,
                               lowerLeftFrontCorner[1] + j * h,
                               lowerLeftFrontCorner[2] + i * h};

        std::array<double, 3> velocity = initialVelocity;

        Particle particle(position, velocity, m, particles.particle_id, epsilon,
                          sigma);
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
  Logger::getInstance().info("New cuboid generated");
}

void ParticleGenerator::insertDisc(const std::array<double, 3> &center,
                                   const std::array<double, 3> &initialVelocity,
                                   size_t radius, double h, double mass,
                                   LinkedCellContainer &particles,
                                   double epsilon, double sigma) {

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
  Logger::getInstance().info("New single molecule generated.");
}
