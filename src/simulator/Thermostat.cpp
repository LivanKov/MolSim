//
// Created by sebastianpse on 12/8/24.
//

#include "Thermostat.h"
#include "particle/ParticleContainer.h"
#include "utils/logger/Logger.h"
#include <cmath>



Thermostat::Thermostat(ParticleContainer& particles,
    double initial_temperature,
    double target_temperature = -1.0,  // when no target temperature got inserted
    double delta_temperature = std::numeric_limits<double>::infinity(),  // Default to infinity
    bool gradual = false,
    size_t dimensions = 3) // by default 3
: particles_(particles),
initial_temperature_(initial_temperature),
delta_temperature_(delta_temperature),
current_temperature_(0),
target_temperature_(target_temperature),
scaling_factor_(1.0),
gradual_(gradual),
dimensions_(dimensions){
    // here we initialize the particles with a velocity with the respective initial temperature
    initialize();
    // No target temperature is provided -> T_target = T_init
    if (target_temperature == -1.0) {
        target_temperature_ = initial_temperature_;
    } else {
        target_temperature_ = target_temperature;
    }
    // Log initialization of thermostat
    Logger::getInstance().info("Thermostat created with initial temperature: " + std::to_string(initial_temperature_));
}


double Thermostat::calculate_kinetic_energy() const {
    double kinetic_energy = 0.0;
    size_t particle_count = particles_.size();
    for (size_t i = 0; i < particle_count; ++i) {
        const Particle& particle = particles_[i];
        auto mass = particle.getM();
        auto velocity = particle.getV();
        auto velocity_dot_product = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
        // using our formula (2) from WS 4, Task 1
        kinetic_energy += 0.5 * mass * velocity_dot_product;
    }
    Logger::getInstance().trace("Calculated kinetic energy of system");
    return kinetic_energy;
}

void Thermostat::calculate_current_temperature() {
    // handle case if no particles are in the system
    if(particles_.size() == 0) {
        Logger::getInstance().warn("There are no particles in the system!");
        return;
    }
    // also ensures that denominator is not 0, additionally checks for valid dimension parameters
    if(dimensions_ == 0 || dimensions_ > 3) {
        Logger::getInstance().warn("Invalid dimension param");
        return;
    }
    current_temperature_ = (calculate_kinetic_energy() * 2) / (dimensions_ *  particles_.size());
    Logger::getInstance().trace("Current temperature calculated");
}

void Thermostat::calculate_scaling_factor(double new_temperature) {
    scaling_factor_ = std::sqrt(new_temperature / current_temperature_);
    Logger::getInstance().debug("Scaling factor calculated: " + std::to_string(scaling_factor_));
}

void Thermostat::initialize_brownian( size_t dimensions) {
    for (auto &particle : particles_) {
        auto mass = particle.getM();
        if(mass <= 0) {
            Logger::getInstance().error("Mass of particle must be positive");
            return;
        }
        // factor for the Maxwell-Boltzmann distribution
        double averageVelocity = std::sqrt(initial_temperature_ / mass);

        // Generate random velocity for the particle
        std::array<double, 3> random_velocity = maxwellBoltzmannDistributedVelocity(averageVelocity, dimensions);

        particle.updateV(random_velocity);
    }
    Logger::getInstance().info("Particles initialized with Brownian motion.");
}

double Thermostat::get_current_temperature() {
    calculate_current_temperature();
    return current_temperature_;
}

void Thermostat::initialize() {
    // we check the temperature according to the particles velocity's
    calculate_current_temperature();
    // we now calculate the scaling factor. we want to have the initial temperature for the system
    // means the particles must have the corresponding velocity
    calculate_scaling_factor(initial_temperature_);
    // avoid loop if no scaling is applied
    if (scaling_factor_ == 1.0) {
        Logger::getInstance().debug("No scaling required, because scaling_factor is 1.0");
        return;
    }
    // in the very end we apply our scaling factor to the particles
    for (auto& p : particles_) {
        auto& current_velocity = p.getV();

        std::array<double, 3> new_velocity{};
        for (size_t i = 0; i < 3; ++i) {
            new_velocity[i] = current_velocity[i] * scaling_factor_;
        }
        p.updateV(new_velocity);
    }
}


void Thermostat::apply(bool gradual) {
    calculate_current_temperature();

    // calculate the temperature difference between current and target temperature
    double temp_diff = target_temperature_ - current_temperature_;

    if (std::abs(temp_diff) < 1e-5) {
        Logger::getInstance().debug("No change in temperature, since current and target are nearly the same");
        return;
    }

    // if we apply the thermostat gradual we ensure that the temperature change is maximum delta T
    if(gradual) {
        if (delta_temperature_ > 0) {
            // checks if temperature difference is greater than delta T
            if (std::abs(temp_diff) > delta_temperature_) {
                // if yes, we calculate our permitted temp_diff in one application of the thermostat
                temp_diff = (temp_diff > 0) ? delta_temperature_ : -delta_temperature_;
            }
        } else {
            Logger::getInstance().error("Invalid delta_temperature");
            return;
        }
    }

    // we now calculate the scaling factor with the permitted temperature difference for one application
    calculate_scaling_factor(current_temperature_ + temp_diff);

    // avoid loop if no scaling is applied
    if (scaling_factor_ == 1.0) {
        Logger::getInstance().debug("No scaling required, because scaling_factor is 1.0");
        return;
    }

    // in the very end we apply our scaling factor to the particles
    for (auto& p : particles_) {
        auto& current_velocity = p.getV();

        std::array<double, 3> new_velocity{};
        for (size_t i = 0; i < 3; ++i) {
            new_velocity[i] = current_velocity[i] * scaling_factor_;
        }
        p.updateV(new_velocity);
    }
}

