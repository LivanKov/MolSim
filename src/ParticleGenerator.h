//
// Created by sebastianpse on 11/9/24.
//

#pragma once

#include "ParticleContainer.h"
#include <array>

class ParticleGenerator {
public:
    ParticleGenerator();

    ParticleContainer generateCuboid(
        const std::array<double, 3>& lowerLeftFrontCorner,
        size_t N1, size_t N2, size_t N3,
        double h,
        double mass,
        const std::array<double, 3>& initialVelocity,
        double averageVelocity
    );

};
