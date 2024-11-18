/*
 * Created by sebastianpse on 11/9/24.
 */

#include "ParticleGenerator.h"
#include "ParticleContainer.h"
#include "utils/MaxwellBoltzmannDistribution.h"

#include <random>

#include "logger/Logger.h"

ParticleGenerator::ParticleGenerator() = default;

void ParticleGenerator::insertCuboid(
    const std::array<double, 3> &lowerLeftFrontCorner,
    const std::array<size_t, 3> &dimensions, double h, double m,
    const std::array<double, 3> &initialVelocity, double averageVelocity,
    ParticleContainer &particles) {
  for (size_t i = 0; i < dimensions[2]; ++i) {
    for (size_t j = 0; j < dimensions[1]; ++j) {
      for (size_t k = 0; k < dimensions[0]; ++k) {
        std::array position = {lowerLeftFrontCorner[0] + k * h,
                               lowerLeftFrontCorner[1] + j * h,
                               lowerLeftFrontCorner[2] + i * h};

        std::array<double, 3> velocity = initialVelocity;
        std::array<double, 3> randomVelocity =
            maxwellBoltzmannDistributedVelocity(averageVelocity, 2);
        for (size_t dim = 0; dim < 3; ++dim) {
          velocity[dim] += randomVelocity[dim];
        }

        Particle particle(position, velocity, m, 0);
        Logger::getInstance().trace("New Particle generated");
        particles.insert(particle);
        Logger::getInstance().trace("New Particle inserted into container");
      }
    }
  }
  Logger::getInstance().info("New cuboid generated");
}
