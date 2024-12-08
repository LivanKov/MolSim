//
// Created by sebastianpse on 12/8/24.
//

#include "Thermostat.h"
#include "particle/ParticleContainer.h"
#include "utils/logger/Logger.h"
#include <cmath>


Thermostat::Thermostat(double initial_temperature,
    size_t number_of_time_stamps,
    double target_temperature,
    double delta_temperature)
        : initial_temperature_(initial_temperature),
number_of_time_stamps_(number_of_time_stamps),
target_temperature_(target_temperature),
delta_temperature_(delta_temperature),
current_temperature_(initial_temperature) {
    Logger::getInstance().info("Thermostat created");
}


double Thermostat::calculate_kinetic_energy(ParticleContainer &particles) const {
    double kinetic_energy = 0.0;
    size_t particle_count = particles.size();
    Logger::getInstance().info("" + particle_count);
    for (size_t i = 0; i < particle_count; ++i) {
        const Particle& particle = particles[i];
        auto mass = particle.getM();
        auto velocity = particle.getV();
        auto velocity_dot_product = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
        // using our formula (2) from WS 4, Task 1
        kinetic_energy += 0.5 * mass * velocity_dot_product;
    }
    Logger::getInstance().trace("Calculated kinetic energy of system");
    return kinetic_energy;
}

double Thermostat::calculate_current_temperature(ParticleContainer &particles, size_t dimension) {
    // handle case if no particles are in the system
    if(particles.size() == 0) {
        Logger::getInstance().warn("There are no particles in the system!");
        return current_temperature_;
    }
    // also ensures that denominator is not 0, additionally checks for valid dimension parameters
    if(dimension == 0 || dimension > 3) {
        Logger::getInstance().warn("Invalid dimension param");
        return current_temperature_;
    }
    current_temperature_ = (calculate_kinetic_energy(particles) * 2) / (dimension * particles.size());
    Logger::getInstance().trace("Current temperature calculated");
    return current_temperature_;
}

void Thermostat::calculate_scaling_factor() {
    scaling_factor_ = std::sqrt(new_temperature_ / current_temperature_);
    Logger::getInstance().debug("Scaling factor calculated: " + std::to_string(scaling_factor_));
}

