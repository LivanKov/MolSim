#include "Simulation.h"
#include "Thermostat.h"
#include "io/input/CheckpointReader.h"
#include "io/input/FileReader.h"
#include "io/input/XMLReader.h"
#include "io/output/CheckpointWriter.h"
#include "io/output/FileWriter.h"
#include "io/output/VTKWriter.h"
#include "io/output/XYZWriter.h"
#include "particle/ParticleGenerator.h"
#include "particle/container/ParticleContainer.h"
#include "simulator/calculations/BoundaryConditions.h"
#include "simulator/calculations/Calculation.h"
#include "simulator/calculations/Force.h"
#include "simulator/calculations/Position.h"
#include "simulator/calculations/Velocity.h"
#include "utils/logger/Logger.h"
#include <chrono>
#include <iostream>
#include <memory>
#include <spdlog/spdlog.h>
#include <utility>

std::unique_ptr<Simulation> Simulation::generate_simulation(SimParams &params) {
  std::unique_ptr<Simulation> ptr = std::make_unique<Simulation>(params);
  return ptr;
}

Simulation::Simulation(SimParams &params) : params_(params) {}

LinkedCellContainer Simulation::readFile(SimParams &params) {
  LinkedCellContainer particles{};
  XMLReader::readXMLFile(particles, params);
  return particles;
}

void Simulation::run(LinkedCellContainer &particles) {

  Logger &logger = Logger::getInstance(params_.log_level);

  logger.info("Starting a simulation with:");
  logger.info("\tEnd time: " + std::to_string(params_.end_time));
  logger.info("\tDelta: " + std::to_string(params_.time_delta));

  int iteration{0};
  double current_time{0};
  size_t total_molecule_updates = 0;

  ForceType FORCE_TYPE = params_.calculate_grav_force
                             ? ForceType::GRAVITATIONAL
                             : (params_.is_membrane ? ForceType::MEMBRANE
                                                    : ForceType::LENNARD_JONES);
  std::unique_ptr<output::FileWriter> writer;
  if (params_.xyz_output) {
    writer = std::make_unique<output::XYZWriter>(particles);
  } else {
    writer = std::make_unique<output::VTKWriter>(particles);
  }

  OPTIONS option =
      params_.linked_cells ? OPTIONS::LINKED_CELLS : OPTIONS::DIRECT_SUM;

  // Initialize Thermostat
  Thermostat thermostat(particles, params_.initial_temp, params_.target_temp,
                        params_.dimensions, params_.delta_temp,
                        params_.is_gradual, params_.enable_brownian);

  if (params_.checkpoint_only) {
    while (current_time < params_.end_time) {
      Calculation<Position>::run(particles, params_.time_delta, option);
      Calculation<BoundaryConditions>::run(particles);
      Calculation<Force>::run(particles, FORCE_TYPE, option);
      Calculation<Velocity>::run(particles, params_.time_delta);
      current_time += params_.time_delta;
    }

    CheckpointWriter::writeCheckpoint(particles, "../output/checkpoint.chk",
                                      params_.time_delta, params_.end_time);
    logger.info("Equilibration completed.");
    return;
  }

  if (params_.resume_from_checkpoint) {
    CheckpointReader::readCheckpoint(particles, params_.time_delta,
                                     params_.resume_start_time);
    logger.info("Resumed from checkpoint. Adding additional input...");
    current_time = params_.resume_start_time;
  }

  // Start measuring time for the main simulation loop
  auto start_time = std::chrono::high_resolution_clock::now();

  while (current_time < params_.end_time) {

    if(current_time >= SimParams::additional_force_time_limit){
      SimParams::apply_fzup = false;
    }

    size_t molecules_this_iteration = particles.size();

    Calculation<Position>::run(particles, params_.time_delta, option);
    Calculation<BoundaryConditions>::run(particles);
    Calculation<Force>::run(particles, FORCE_TYPE, option);
    Calculation<Velocity>::run(particles, params_.time_delta);

    total_molecule_updates += molecules_this_iteration;

    // Apply the thermostat periodically
    if (SimParams::enable_thermo && iteration % params_.n_thermostats == 0) {
      thermostat.apply();
      logger.info("Thermostat applied at iteration: " +
                  std::to_string(iteration));
    }

    iteration++;
    if (iteration % params_.write_frequency == 0 && !params_.disable_output) {
      writer->plot_particles(params_.output_path, iteration);
    }

    logger.info("Iteration " + std::to_string(iteration) + " finished.");
    current_time += params_.time_delta;
  }

  // End measuring time for the main simulation loop
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> runtime = end_time - start_time;

  // Calculate updates per second
  double updates_per_second = total_molecule_updates / runtime.count();

  std::cout << "Total runtime: " << std::to_string(runtime.count())
            << " seconds" << std::endl;
  std::cout << "Molecules updated per second: "
            << std::to_string(updates_per_second) << std::endl;

  logger.info("output written. Terminating...");

  logger.info("Number of particles: " + std::to_string(particles.size()));

  logger.info("Particles left the domain: " +
              std::to_string(particles.particles_left_domain));

  logger.info("Amount of halo particles:" +
              std::to_string(particles.halo_count));

  logger.info("Simulation finished.");
};