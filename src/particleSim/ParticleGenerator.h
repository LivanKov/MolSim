/*
 * Created by sebastianpse on 11/9/24.
 */

#pragma once

#include "ParticleContainer.h"
#include <array>

class ParticleGenerator {
public:
  ParticleGenerator();

  static ParticleContainer
  generateCuboid(const std::array<double, 3> &lowerLeftFrontCorner,
                 const std::array<size_t, 3> &dimensions, double h, double mass,
                 const std::array<double, 3> &initialVelocity,
                 double averageVelocity);
};
