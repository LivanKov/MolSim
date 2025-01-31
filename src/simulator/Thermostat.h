//
// Created by sebastianpse on 12/8/24.
//
#pragma once

#include "Thermostat.h"
#include "particle/container/LinkedCellContainer.h"
#include "utils/MaxwellBoltzmannDistribution.h"

class Thermostat {

public:
  /**
   * @brief Constructor
   *
   * @param particles The particles of the system.
   * @param initial_temperature The initial temperature of the system.
   * @param target_temperature The temperature that the system should reach.
   * @param dimensions Number of dimensions in the system.
   * @param delta_temperature Maximum absolute temperature change per
   * application.
   * @param gradual Whether thermostat should apply gradually or directly
   * towards/to target temperature.
   * @param enable_brownian Whether particles should be initialized with
   * brownian motion.
   */
  Thermostat(LinkedCellContainer &particles, double initial_temperature,
             double target_temperature =
                 -1.0,              // when no target temperature got inserted
             size_t dimensions = 3, // by default 3 dimensions
             double delta_temperature =
                 std::numeric_limits<double>::infinity(), // Default to infinity
             bool gradual = true,                         // default is true
             bool enable_brownian =
                 true); // by default initialized with brownian motion

  /**
   * @brief This is the application method of the thermostat.
   *
   * It can heat, cool or hold the temperature of the systems with respect to
   * the parameters specified. This is done by scaling the particles' velocities
   * appropriate.
   */
  void apply();

  /**
   * @brief Calculates the current kinetic energy of the system.
   *
   * Here the formula for the kinetic energy of a systems gets applied.
   *
   * @return double The current kinetic energy.
   */
  double calculate_kinetic_energy() const;

  /**
   * @brief Calculates the current temperature of the system
   *
   * Uses the current kinetic energy to calculate the current temperature. We
   * use k_B = 1.
   *
   * @return The calculated current temperature.
   */
  void calculate_current_temperature();

  /**
   * @brief Calculates the velocity scaling factor.
   *
   * Applying the scaling factor formula here. Is used to calculate new velocity
   * of particles.
   *
   * @param new_temperature To this temperature should the particles' velocities
   * gets scaled.
   */
  void calculate_scaling_factor(double new_temperature);

  /**
   * @brief Initializes particles_ with brownian motion.
   *
   * This method ensures that the particles are getting initialized with
   * brownian motion when enable_brownian = true. The method is used in the
   * constructor.
   */
  void initialize_brownian();

  /**
   * @brief Gets the current temperature of the system.
   * @return double The current temperature.
   */
  double get_current_temperature() const;

  /**
   * @brief Gets the number of dimensions for the system.
   * @return size_t The number of dimensions.
   */
  size_t get_dimensions() const;

  /**
   * @brief Gets the target temperature for the thermostat.
   * @return double The target temperature.
   */
  double get_target_temperature() const;

  /**
   * @brief Gets the gradual setting of the thermostat.
   * @return bool `true` if gradual application is enabled, `false` otherwise.
   */
  bool get_gradual() const;

  /**
   * @brief Scales particles' velocities directly with respect to initial
   * temperature.
   *
   * Helper method primarily for testing purposes. Scales particles_ directly to
   * the velocities in regard to the initial temperature.
   */
  void initialize();

  // ---------------- modified methods: -----------------------

  /**
   * @brief This is the modified application method of the thermostat.
   *
   * It can heat, cool or hold the temperature of the systems with respect to
   * the parameters specified. This is done by scaling the particles' thermal
   * motions appropriate.
   */
  void apply_new();

  /**
   * @brief Calculates the average velocity of the particles in the system.
   *
   * @return std::array<double, 3> The average velocity
   */
  std::array<double, 3> calculate_average_velocity();

  /**
   *@brief Calculates the current thermal motion of each particle.
   *
   * When calling this method each particle's thermal motion gets updated
   * according to its current velocity and the average velocity of all
   * particles: thermal_motion = currentV - averageV
   *
   */
  void determine_thermal_motion();

  /**
   * @brief Calculates the current kinetic energy of the system.
   *
   * Here the formula for the kinetic energy of a systems gets applied.
   * The kinetic energy just depends, besides on the mass and
   * number of particles, on the particles' thermal motion.
   *
   * @return double The current kinetic energy.
   */
  double calculate_kinetic_energy_new();

  /**
   * @brief Gets the current temperature of the system.
   *
   * The current temperature gets calculated partly with the kinetic energy.
   * The kinetic energy only
   *
   * @return double The current temperature.
   */
  void calculate_current_temperature_new();

  /**
   * @brief Scales particles' velocities directly with respect to initial
   * temperature.
   *
   * Helper method primarily for testing purposes. Scales particles_ directly to
   * the velocities in regard to the initial temperature which is calculated
   * by thermal motion.
   */
  void initialize_new();

private:
  /** @brief Reference to the ParticleContainer holding all particles in the
   * system. */
  LinkedCellContainer &particles_;

  /** @brief The initial temperature of the system. */
  double initial_temperature_;

  /** @brief The target temperature of the system. */
  double target_temperature_;

  /** @brief The number of dimensions for the system (valid are 1, 2 or 3). */
  size_t dimensions_;

  /** @brief Maximum absolute temperature change per thermostat application. */
  double delta_temperature_;

  /** @brief Whether thermostat should apply gradually or directly towards/to
   * target temperature. */
  bool gradual_;

  /** @brief Scaling factor used to adjust particles' velocities based on
   * temperature. */
  double scaling_factor_;

  /** @brief The current temperature of the system. */
  double current_temperature_;
};
