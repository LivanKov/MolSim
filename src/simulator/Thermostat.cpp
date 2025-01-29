//
// Created by sebastianpse on 12/8/24.
//

#include "Thermostat.h"
#include "io/input/cli/SimParams.h"
#include "particle/container/LinkedCellContainer.h"
#include "utils/logger/Logger.h"
#include <cmath>
#include <iostream>

Thermostat::Thermostat(LinkedCellContainer &particles,
                       double initial_temperature, double target_temperature,
                       size_t dimensions, double delta_temperature,
                       bool gradual, bool enable_brownian)
    : particles_(particles), initial_temperature_(initial_temperature),
      target_temperature_(target_temperature), dimensions_(dimensions),
      delta_temperature_(delta_temperature), gradual_(gradual),
      scaling_factor_(1.0), current_temperature_(0) {

  if (SimParams::enable_thermo) {
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
    if (enable_brownian)
      initialize_brownian();

    // No target temperature is provided -> T_target = T_init
    target_temperature_ = (target_temperature == -1.0) ? initial_temperature_
                                                       : target_temperature;

    Logger::getInstance().info(
        "Thermostat created with initial temperature: " +
        std::to_string(initial_temperature_) +
        ", target temperature: " + std::to_string(target_temperature_) +
        ", delta temperature: " + std::to_string(delta_temperature_) +
        ", dimensions: " + std::to_string(dimensions) +
        ", Brownian motion: " + (enable_brownian ? "enabled" : "disabled"));
  }
}

// application method of thermostat

void Thermostat::apply() {
  calculate_current_temperature();

  // calculate the temperature difference between current and target temperature
  double temp_diff = target_temperature_ - current_temperature_;

  if (std::abs(temp_diff) < 1e-5) {
    Logger::getInstance().debug("No change in temperature, since current and "
                                "target are nearly the same");
    return;
  }

  // if we apply the thermostat gradual we ensure that the temperature change is
  // maximum delta T
  if (gradual_) {
    // checks if temperature difference is greater than delta T
    if (std::abs(temp_diff) > delta_temperature_) {
      // if yes, we calculate our permitted temp_diff in one application of the
      // thermostat
      temp_diff = (temp_diff > 0) ? delta_temperature_ : -delta_temperature_;
    }
  }

  // we now calculate the scaling factor with the permitted temperature
  // difference for one application
  calculate_scaling_factor(current_temperature_ + temp_diff);

  // avoid loop if no scaling is applied
  if (scaling_factor_ == 1.0) {
    Logger::getInstance().debug(
        "No scaling required, because scaling_factor is 1.0");
    return;
  }

  // in the very end we apply our scaling factor to the particles
  for (size_t i = 0; i < particles_.size(); ++i) {
    auto &p = particles_[i];
    auto &current_velocity = p.getV();

    std::array<double, 3> new_velocity{};
    for (size_t j = 0; j < dimensions_; ++j) {
      new_velocity[j] = current_velocity[j] * scaling_factor_;
    }

    p.updateV(new_velocity);
  }

  Logger::getInstance().debug(
      "Thermostat applied. The new temperature is now: " +
      std::to_string(get_current_temperature()));
}

// ----------------- helper methods -----------------------------------

double Thermostat::calculate_kinetic_energy() const {
  double kinetic_energy = 0.0;
  size_t particle_count = particles_.size();
  for (size_t i = 0; i < particle_count; ++i) {
    const Particle &particle = particles_[i];
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
  if (particles_.size() == 0) {
    Logger::getInstance().warn("There are no particles in the system!");
    return;
  }
  // check for dimensions param in constructor ensures that denominator is not 0
  current_temperature_ =
      (calculate_kinetic_energy() * 2) / (dimensions_ * particles_.size());
  Logger::getInstance().trace("Current temperature calculated");
}

void Thermostat::calculate_scaling_factor(double new_temperature) {
  if (current_temperature_ < 1e-6) {
    Logger::getInstance().warn(
        "Current temperature is near 0 or 0. Scaling factor may be too large!");
  }
  scaling_factor_ = std::sqrt(new_temperature / current_temperature_);
  Logger::getInstance().debug("Scaling factor calculated: " +
                              std::to_string(scaling_factor_));
}

void Thermostat::initialize_brownian() {
  for (size_t i = 0; i < particles_.size(); ++i) {
    auto &particle = particles_[i];
    auto mass = particle.getM();
    if (mass <= 0) {
      Logger::getInstance().error("Mass of particle must be positive");
      return;
    }

    // factor for the Maxwell-Boltzmann distribution
    double average_velocity = std::sqrt(initial_temperature_ / mass);

    // Generate random velocity for the particle
    std::array<double, 3> random_velocity =
        maxwellBoltzmannDistributedVelocity(average_velocity, dimensions_);

    auto current_velocity = particle.getV();

    for (size_t j = 0; j < dimensions_; ++j) {
      current_velocity[j] += random_velocity[j];
    }

    particle.updateV(current_velocity);
  }
  Logger::getInstance().info("Particles initialized with Brownian motion.");
}

// -------- for testing purposes -------------------------------------

void Thermostat::initialize() {
  // we check the temperature according to the particles velocity's
  calculate_current_temperature();
  // we now calculate the scaling factor. we want to have the initial
  // temperature for the system means the particles must have the corresponding
  // velocity
  calculate_scaling_factor(initial_temperature_);
  // avoid loop if no scaling is applied
  if (scaling_factor_ == 1.0) {
    Logger::getInstance().debug(
        "No scaling required, because scaling_factor is 1.0");
    return;
  }
  // in the very end we apply our scaling factor to the particles
  for (size_t i = 0; i < particles_.size(); ++i) {
    auto &p = particles_[i];
    auto &current_velocity = p.getV();

    std::array<double, 3> new_velocity{};
    for (size_t j = 0; j < 3; ++j) {
      new_velocity[j] = current_velocity[j] * scaling_factor_;
    }
    p.updateV(new_velocity);
  }
}


// ----------------------- thermostat modification: -----------------------------------------

// new application method of thermostat (using thermal motion)

void Thermostat::apply_new() {
  // we use our new calculate current temperature method
  calculate_current_temperature_new();
  // we also store the current average velocity for later
  auto average_velocity = calculate_average_velocity();

  // calculate the temperature difference between current and target temperature
  double temp_diff = target_temperature_ - current_temperature_;

  if (std::abs(temp_diff) < 1e-5) {
    Logger::getInstance().debug("No change in temperature, since current and "
                                "target are nearly the same");
    return;
  }

  // if we apply the thermostat gradual we ensure that the temperature change is
  // maximum delta T
  if (gradual_) {
    // checks if temperature difference is greater than delta T
    if (std::abs(temp_diff) > delta_temperature_) {
      // if yes, we calculate our permitted temp_diff in one application of the
      // thermostat
      temp_diff = (temp_diff > 0) ? delta_temperature_ : -delta_temperature_;
    }
  }

  // we now calculate the scaling factor with the permitted temperature
  // difference for one application
  calculate_scaling_factor(current_temperature_ + temp_diff);

  // avoid loop if no scaling is applied
  if (scaling_factor_ == 1.0) {
    Logger::getInstance().debug(
        "No scaling required, because scaling_factor is 1.0");
    return;
  }

  // in the very end we apply our scaling factor to the particles
  // we only scale the thermal motion of the particle
  for (size_t i = 0; i < particles_.size(); ++i) {
    auto &particle = particles_[i];
    auto &current_thermal_motion= particle.getThermalMotion();

    // here we scale only the thermal motion part of the particle
    std::array<double, 3> scaled_thermal_motion{};
    for (size_t j = 0; j < dimensions_; ++j) {
      scaled_thermal_motion[j] = current_thermal_motion[j] * scaling_factor_;
    }

    std::array<double,3> new_velocity{};
    // then we add the average velocity to the scaled thermal motion
    // the result is the new velocity of the particle
    for (size_t l = 0; l < dimensions_; ++l) {
      new_velocity[l] = average_velocity[l] + scaled_thermal_motion[l];
    }
    // here we update the new velocity
    particle.updateV(new_velocity);
  }

  Logger::getInstance().debug(
      "Thermostat applied. The new temperature is now: " +
      std::to_string(get_current_temperature()));
}

// ------------ new helper methods --------------

std::array<double, 3> Thermostat::calculate_average_velocity() {
  const size_t amount_particles = particles_.size();
  std::array<double, 3> average_velocity = {0.,0.,0.};

  // also to avoid dividing by zero later
  if (amount_particles == 0) {
    Logger::getInstance().error("There are no particles in the system!");
    return average_velocity;
  }

  for (size_t i = 0; i < amount_particles; ++i) {
    const auto &velocity = particles_[i].getV();
    for (size_t j = 0; j < dimensions_; ++j) {
      average_velocity[j] += velocity[j];
    }
  }

  for (size_t i = 0; i < dimensions_; ++i) {
    average_velocity[i] /= static_cast<double>(amount_particles);
  }

  return average_velocity;
}

void Thermostat::determine_thermal_motion() {

  auto average_velocity = calculate_average_velocity();

  for (size_t i = 0; i < particles_.size(); ++i) {
    auto &particle = particles_[i];
    std::array<double, 3> thermal_motion = {0.,0.,0.};

    for (size_t l = 0; l < dimensions_; ++l) {
      thermal_motion[l] = particle.getV()[l] - average_velocity[l];
    }
    particle.updateThermalMotion(thermal_motion);
  }
}

double Thermostat::calculate_kinetic_energy_new() {
  // here, we need to determine the current thermal motion of each particle
  // in other words: to update it
  determine_thermal_motion();
  double kinetic_energy = 0.0;
  size_t particle_count = particles_.size();
  for (size_t i = 0; i < particle_count; ++i) {
    const Particle &particle = particles_[i];
    auto mass = particle.getM();
    std::array<double, 3> thermal_motion = particle.getThermalMotion();
    double velocity_dot_product = 0.0;
    for (size_t l = 0; l < dimensions_; ++l) {
      // only use thermal motion for calculating kinetic energy
      velocity_dot_product += thermal_motion[l] * thermal_motion[l];
    }
    // using our formula (2) from WS 4, Task 1
    kinetic_energy += 0.5 * mass * velocity_dot_product;
  }
  Logger::getInstance().trace("Calculated kinetic energy of system");
  return kinetic_energy;
}

void Thermostat::calculate_current_temperature_new() {
  // handle case if no particles are in the system
  if (particles_.size() == 0) {
    Logger::getInstance().warn("There are no particles in the system!");
    return;
  }
  // check for dimensions param in constructor ensures that denominator is not 0
  current_temperature_ =
      (calculate_kinetic_energy_new() * 2) / (dimensions_ * particles_.size());
  Logger::getInstance().trace("Current temperature calculated");
}

// -------- for testing purposes -------------------------------------

void Thermostat::initialize_new() {
  // we check the temperature according to the particles' velocities
  calculate_current_temperature_new();
  // and average velocity
  auto average_velocity = calculate_average_velocity();
  // we now calculate the scaling factor. we want to have the initial
  // temperature for the system means the particles must have the corresponding
  // velocity
  calculate_scaling_factor(initial_temperature_);
  // avoid loop if no scaling is applied
  if (scaling_factor_ == 1.0) {
    Logger::getInstance().debug(
        "No scaling required, because scaling_factor is 1.0");
    return;
  }
  // in the very end we apply our scaling factor to the particles
  // we only scale the thermal motion of the particle
  for (size_t i = 0; i < particles_.size(); ++i) {
    auto &particle = particles_[i];
    auto &current_thermal_motion= particle.getThermalMotion();

    // here we scale only the thermal motion part of the particle
    std::array<double, 3> scaled_thermal_motion{};
    for (size_t j = 0; j < dimensions_; ++j) {
      scaled_thermal_motion[j] = current_thermal_motion[j] * scaling_factor_;
    }

    std::array<double,3> new_velocity{};
    // then we add the average velocity to the scaled thermal motion
    // the result is the new velocity of the particle
    for (size_t l = 0; l < dimensions_; ++l) {
      new_velocity[l] = average_velocity[l] + scaled_thermal_motion[l];
    }
    // here we update the new velocity
    particle.updateV(new_velocity);
  }
}


// ------------- getters  --------------------------------------------

double Thermostat::get_current_temperature() const { return current_temperature_; }

size_t Thermostat::get_dimensions() const { return dimensions_; }

double Thermostat::get_target_temperature() const { return target_temperature_; }

bool Thermostat::get_gradual() const { return gradual_; }