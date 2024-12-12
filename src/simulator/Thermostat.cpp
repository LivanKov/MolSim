//
// Created by sebastianpse on 12/8/24.
//

#include "Thermostat.h"
#include "particle/ParticleContainer.h"
#include "utils/logger/Logger.h"
#include <cmath>



Thermostat::Thermostat(
    ParticleContainer& particles,
    double initial_temperature,
    double target_temperature,
    size_t dimensions,
    double delta_temperature,
    bool gradual,
    bool enable_brownian)
    : particles_(particles),
    initial_temperature_(initial_temperature),
    target_temperature_(target_temperature),
    dimensions_(dimensions),
    delta_temperature_(delta_temperature),
    gradual_(gradual),
    scaling_factor_(1.0),
    current_temperature_(0) {

    // checks for correct dimensions param
    if (dimensions < 1 || dimensions > 3) {
        Logger::getInstance().error("Invalid parameter for dimensions!");
        throw std::invalid_argument("Dimensions must be between 1 and 3!");
    }

    // checks for valid initial temperature
    if (initial_temperature <= 0) {
        Logger::getInstance().error("Invalid initial temperature!");
        throw std::invalid_argument("Initial temperature must be positive!");
    }

    // checks if delta T is valid
    if (delta_temperature < 0) {
        Logger::getInstance().error("Invalid delta_temperature!");
        throw std::invalid_argument("delta_temperature must be non-negative!");
    }

    // here we initialize the particles with brownian motion if it's enabled
    if(enable_brownian)
        initialize_brownian();

    // No target temperature is provided -> T_target = T_init
    target_temperature_ = (target_temperature == -1.0) ? initial_temperature_ : target_temperature;

    Logger::getInstance().info(
        "Thermostat created with initial temperature: " +
        std::to_string(initial_temperature_) +
        ", target temperature: " + std::to_string(target_temperature_) +
        ", delta temperature: " + std::to_string(delta_temperature_) +
        ", dimensions: " + std::to_string(dimensions) +
        ", Brownian motion: " + (enable_brownian ? "enabled" : "disabled")
    );
}


double Thermostat::calculate_kinetic_energy() const {
    double kinetic_energy = 0.0;
    size_t particle_count = particles_.size();
    for (size_t i = 0; i < particle_count; ++i) {
        const Particle& particle = particles_[i];
        auto mass = particle.getM();
        auto velocity = particle.getV();
        double velocity_dot_product = 0.0;
        for (size_t j = 0; j < dimensions_; ++j) {
            velocity_dot_product += velocity[j] * velocity[j];
        }
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
    // check for dimensions param in constructor ensures that denominator is not 0
    current_temperature_ = (calculate_kinetic_energy() * 2) / (dimensions_ *  particles_.size());
    Logger::getInstance().trace("Current temperature calculated");
}


void Thermostat::calculate_scaling_factor(double new_temperature) {
    if (current_temperature_ < 1e-6) {
        Logger::getInstance().warn("Current temperature is near 0 or 0. Scaling factor may be too large!");
    }
    scaling_factor_ = std::sqrt(new_temperature / current_temperature_);
    Logger::getInstance().debug("Scaling factor calculated: " + std::to_string(scaling_factor_));
}


void Thermostat::initialize_brownian() {
    for (auto &particle : particles_) {
        auto mass = particle.getM();
        if(mass <= 0) {
            Logger::getInstance().error("Mass of particle must be positive");

            return;
        }
        // factor for the Maxwell-Boltzmann distribution
        double average_velocity = std::sqrt(initial_temperature_ / mass);

        // Generate random velocity for the particle
        std::array<double, 3> random_velocity = maxwellBoltzmannDistributedVelocity(average_velocity, dimensions_);

        auto current_velocity = particle.getV();

        for (size_t i = 0; i < dimensions_; ++i) {
            current_velocity[i] += random_velocity[i];
        }
        particle.updateV(current_velocity);
    }
    Logger::getInstance().info("Particles initialized with Brownian motion.");
}


double Thermostat::get_current_temperature() {
    calculate_current_temperature();
    return current_temperature_;
}




void Thermostat::apply() {
    calculate_current_temperature();

    // calculate the temperature difference between current and target temperature
    double temp_diff = target_temperature_ - current_temperature_;

    if (std::abs(temp_diff) < 1e-5) {
        Logger::getInstance().debug("No change in temperature, since current and target are nearly the same");
        return;
    }

    // if we apply the thermostat gradual we ensure that the temperature change is maximum delta T
    if(gradual_) {
        // checks if temperature difference is greater than delta T
        if (std::abs(temp_diff) > delta_temperature_) {
            // if yes, we calculate our permitted temp_diff in one application of the thermostat
            temp_diff = (temp_diff > 0) ? delta_temperature_ : -delta_temperature_;
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
        for (size_t i = 0; i < dimensions_; ++i) {
            new_velocity[i] = current_velocity[i] * scaling_factor_;
        }
        p.updateV(new_velocity);
    }
    Logger::getInstance().debug("Thermostat applied. The new temperature is now: " + std::to_string(get_current_temperature()));
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

size_t Thermostat::get_dimensions() const {
    return dimensions_;
}

double Thermostat::get_target_temperature() const {
    return target_temperature_;
}

bool Thermostat::get_gradual() const {
    return gradual_;
}


