/*
 * Created by sebastianpse on 11/9/24.
 */

#pragma once

#include "ParticleContainer.h"
#include <array>

/**
 * @class ParticleGenerator.
 * @brief Generator that allows to create specially designed ParticleContainer
 * objects.
 */
class ParticleGenerator {
public:
  /**
   * @brief Constructor.
   */
  ParticleGenerator();

  /**
   * @brief Generate a cuboid consisting of particles, based on arguments.
   * @param lowerFrontCorner Coordinates of the lower left point of the cuboid
   * in 3d space.
   * @param dimensions Amount of particles in height, width and depth of the
   * cuboid.
   * @param h Distance between particles.
   * @param mass Mass of an individual particle.
   * @param initialVelocity contains an array with individual velocity in 3d
   * space. Applied to all particles in the container.
   * @param averageVelocity Mean value of the velocity of the Brownian Motion.
   * @return ParticleContainer in accordance to the arguments passed.
   */
  static void insertCuboid(const std::array<double, 3> &lowerLeftFrontCorner,
                           const std::array<size_t, 3> &dimensions, double h,
                           double mass,
                           const std::array<double, 3> &initialVelocity,
                           double averageVelocity,
                           ParticleContainer &particles);

  /**
   * @brief Generate a disc of particles. The disc gets plotted along the
   * XY-plane.
   * @param center Coordinates of the center.
   * @param initialVelocity Contains an array with individual velocity in 3d
   * space. Applied to all particles in the container.
   * @param radius Number of particles along the radius.
   * @param h Distance between particles.
   * @param mass Mass of an individual particle.
   * @param particles Container of particles that form the disc.
   * @return ParticleContainer in accordance to the arguments passed.
   */
  static void insertDisc(const std::array<double, 3> &center,
                         const std::array<double, 3> &initialVelocity,
                         size_t radius, double h, double mass,
                         ParticleContainer &particles);
};
