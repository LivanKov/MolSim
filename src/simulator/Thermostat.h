//
// Created by sebastianpse on 12/8/24.
//
#pragma once

# include "Thermostat.h"
#include "particle/ParticleContainer.h"

class Thermostat {

public:


    Thermostat(double initial_temperature, size_t number_of_time_stamps, double target_temperature, double delta_temperature);

    /**
     * @brief Calculates the current kinetic energy of the system
     * @param particles The particles of the system
     * @return The current kinetic energy.
     */
    double calculate_kinetic_energy(ParticleContainer &particles) const;

    /**
     * @brief Calculates the current temperature of the system
     *
     * Uses the current kinetic energy to calculate the current temperature. We use k_B = 1.
     *
     * @param particles The particles of the system
     * @param dimension Dimension of system
     * @return The calculated current temperature.
     */
    double calculate_current_temperature(ParticleContainer &particles, size_t dimension);

    /**
     * @brief Calculates the temperature scaling factor. Is used to calculate new velocity of particles.
     */
    void calculate_scaling_factor();


private:
    double initial_temperature_;         /*! initial temperature of system */
    double target_temperature_;         /*! the final temperature that system should have */
    double current_temperature_;        /*! current temperature of system */
    double new_temperature_;            /*! target temperature for one increment */
    double scaling_factor_;             /*! scaling factor for velocities based on new and current temperature */
    double delta_temperature_;          /*! maximal absolute temperature change that is allowed per application */
    size_t number_of_time_stamps_;      /*! number of time stamps after thermostat gets applied */
    bool gradual_;                      /*! is true when the velocity scaling should be gradual */

};
