//
// Created by sebastianpse on 12/8/24.
//
#pragma once

# include "Thermostat.h"
#include "particle/ParticleContainer.h"
#include "utils/MaxwellBoltzmannDistribution.h"

class Thermostat {

public:


    Thermostat(
    ParticleContainer& particles,
    double initial_temperature,
    double target_temperature = -1.0,  // when no target temperature got inserted
    size_t dimensions = 3, // by default 3 dimensions
    double delta_temperature = std::numeric_limits<double>::infinity(),  // Default to infinity
    bool gradual = false, // default is false
    bool enable_brownian = true); // by default initialized with brownian motion


    /**
     * @brief Calculates the current kinetic energy of the system
     * @return The current kinetic energy.
     */
    double calculate_kinetic_energy() const;

    /**
     * @brief Calculates the current temperature of the system
     *
     * Uses the current kinetic energy to calculate the current temperature. We use k_B = 1.
     *
     * @return The calculated current temperature.
     */
    void calculate_current_temperature();

    /**
     * @brief Calculates the temperature scaling factor. Is used to calculate new velocity of particles.
     */
    void calculate_scaling_factor(double new_temperature);

    void initialize_brownian();

    double get_current_temperature();

    size_t get_dimensions() const;

    double get_target_temperature() const;

    bool get_gradual() const;

    void initialize();


    void apply();




private:
    ParticleContainer& particles_;
    double initial_temperature_;        /*! initial temperature of system */
    double target_temperature_;         /*! the final temperature that system should have */
    size_t dimensions_;
    double delta_temperature_;          /*! maximal absolute temperature change that is allowed per application */
    bool gradual_;
    double scaling_factor_;             /*! scaling factor for velocities based on new and current temperature */
    double current_temperature_;        /*! current temperature of system */



};
