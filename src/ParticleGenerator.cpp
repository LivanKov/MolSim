//
// Created by sebastianpse on 11/9/24.
//

#include "ParticleGenerator.h"
#include "ParticleContainer.h"
#include "utils/MaxwellBoltzmannDistribution.h"


#include <random>



ParticleGenerator::ParticleGenerator() = default;


ParticleContainer ParticleGenerator::generateCuboid(
    const std::array<double, 3>& lowerLeftFrontCorner,
    const std::array<size_t, 3>& dimensions,
    double h,
    double m,
    const std::array<double, 3>& initialVelocity,
    double averageVelocity
) {
    ParticleContainer container;
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < dimensions[2]; ++k) {
                std::array position = {
                    lowerLeftFrontCorner[0] + i * h,
                    lowerLeftFrontCorner[1] + j * h,
                    lowerLeftFrontCorner[2] + k * h
                };


                std::array<double, 3> velocity = initialVelocity;
                std::array<double, 3> randomVelocity = maxwellBoltzmannDistributedVelocity(averageVelocity, 3);
                for (size_t dim = 0; dim < 3; ++dim) {
                    velocity[dim] += randomVelocity[dim];
                }


                Particle particle(position, velocity, m, 0);
                container.insert(particle);
            }
        }
    }
    return container;
}



