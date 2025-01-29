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
#include "particle/container/DirectSumContainer.h"
#include "simulator/calculations/BoundaryConditions.h"
#include "simulator/calculations/Calculation.h"
#include "simulator/calculations/Force.h"
#include "simulator/calculations/Position.h"
#include "simulator/calculations/Velocity.h"
#include "utils/logger/Logger.h"
#include <chrono>
#include <iostream>
#include <memory>
#include <omp.h>
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
  // Set OpenMP Threads
  if (SimParams::enable_omp) {
    omp_set_num_threads(omp_get_num_procs());
    std::cout << "Using " << omp_get_max_threads()
              << " threads for the simulation." << std::endl;
    std::cout << "Using " << to_string(SimParams::ompstrategy)
              << " Strategy for the simulation." << std::endl;
  }

  // Initialize Logger
  Logger &logger = Logger::getInstance(params_.log_level);
  logger.info("Starting a simulation with:");
  logger.info("\tEnd time: " + std::to_string(params_.end_time));
  logger.info("\tDelta: " + std::to_string(params_.time_delta));

  int iteration{0};
  double current_time{0};
  size_t total_molecule_updates = 0;
  ForceType FORCE_TYPE = params_.calculate_grav_force
                             ? ForceType::GRAVITATIONAL
                             : ForceType::LENNARD_JONES;

  std::unique_ptr<output::FileWriter> writer = createFileWriter(particles);

  OPTIONS option =
      params_.linked_cells ? OPTIONS::LINKED_CELLS : OPTIONS::DIRECT_SUM;

  // Initialize Thermostat
  Thermostat thermostat(particles, params_.initial_temp, params_.target_temp,
                        params_.dimensions, params_.delta_temp,
                        params_.is_gradual, params_.enable_brownian);
  // Checkout-only mode
  if (params_.checkpoint_only) {
    checkpointMode(particles, current_time, option, FORCE_TYPE);
    return;
  }

  // Resume from checkpoint if enabled
  if (params_.resume_from_checkpoint) {
    CheckpointReader::readCheckpoint(particles, params_.time_delta,
                                     params_.resume_start_time);
    logger.info("Resumed from checkpoint. Adding additional input...");
    current_time = params_.resume_start_time;
  }

  // Main simulation loop
  simulate(particles, current_time, iteration, total_molecule_updates, writer,
           thermostat, option, FORCE_TYPE);

  logger.info("Simulation finished.");
}

// ----------------- Helper Functions -----------------------------------

std::unique_ptr<output::FileWriter>
Simulation::createFileWriter(LinkedCellContainer &particles) const {
  if (params_.xyz_output) {
    return std::make_unique<output::XYZWriter>(particles);
  } else {
    return std::make_unique<output::VTKWriter>(particles);
  }
}

void Simulation::checkpointMode(LinkedCellContainer &particles,
                                double &current_time, OPTIONS option,
                                ForceType force_type) {
  while (current_time < params_.end_time) {
    Calculation<Position>::run(particles, params_.time_delta, option);
    Calculation<BoundaryConditions>::run(particles);
    Calculation<Force>::run(particles, force_type, option);
    Calculation<Velocity>::run(particles, params_.time_delta);
    current_time += params_.time_delta;
  }

  CheckpointWriter::writeCheckpoint(particles, "../output/checkpoint.chk",
                                    params_.time_delta, params_.end_time);
  Logger::getInstance().info("Equilibration completed.");
}

void Simulation::simulate(LinkedCellContainer &particles, double &current_time,
                          int &iteration, size_t &total_molecule_updates,
                          std::unique_ptr<output::FileWriter> &writer,
                          Thermostat &thermostat, OPTIONS option,
                          ForceType force_type) {
  auto start_time = std::chrono::high_resolution_clock::now();

  while (current_time < params_.end_time) {
    size_t molecules_this_iteration = particles.size();

    Calculation<Position>::run(particles, params_.time_delta, option);
    Calculation<BoundaryConditions>::run(particles);
    Calculation<Force>::run(particles, force_type, option);
    Calculation<Velocity>::run(particles, params_.time_delta);

    total_molecule_updates += molecules_this_iteration;

    if (SimParams::enable_thermo && iteration % params_.n_thermostats == 0) {
      thermostat.apply_new();
      Logger::getInstance().warn("Thermostat applied at iteration: " +
                                 std::to_string(iteration));
    }

    iteration++;
    if (iteration % params_.write_frequency == 0 && !params_.disable_output) {
      writer->plot_particles(params_.output_path, iteration);
    }

    Logger::getInstance().info("Iteration " + std::to_string(iteration) +
                               " finished.");
    current_time += params_.time_delta;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> runtime = end_time - start_time;

  logPerformance(runtime, total_molecule_updates);
}

void Simulation::logPerformance(const std::chrono::duration<double> &runtime,
                                size_t total_molecule_updates) const {
  double updates_per_second = total_molecule_updates / runtime.count();

  Logger::getInstance().info("output written. Terminating...");
  std::cout << "Total runtime: " << runtime.count() << " seconds" << std::endl;
  std::cout << "Molecules updated per second: " << updates_per_second
            << std::endl;
}
