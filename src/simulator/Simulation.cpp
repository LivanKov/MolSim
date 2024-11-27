#include "Simulation.h"
#include "io/input/FileReader.h"
#include "io/output/FileWriter.h"
#include "io/output/VTKWriter.h"
#include "io/output/XYZWriter.h"
#include "particle/ParticleContainer.h"
#include "simulator/calculations/Calculation.h"
#include "simulator/calculations/Force.h"
#include "simulator/calculations/Position.h"
#include "simulator/calculations/Velocity.h"
#include "utils/logger/Logger.h"
#include <memory>
#include <spdlog/spdlog.h>
#include <utility>

std::unique_ptr<Simulation> Simulation::generate_simulation(SimParams &params) {
  std::unique_ptr<Simulation> ptr = std::make_unique<Simulation>(params);
  return ptr;
}

Simulation::Simulation(SimParams &params) : params_(params) {}

void Simulation::run() {

  Logger &logger = Logger::getInstance(params_.log_level);

  logger.info("Starting a simulation with:");
  logger.info("\tEnd time: " + std::to_string(params_.end_time));
  logger.info("\tDelta: " + std::to_string(params_.time_delta));

  ParticleContainer particles{};

  FileReader::readFile(particles, params_.input_path.data());

  int iteration{0};
  double current_time{0};

  ForceType FORCE_TYPE =
      params_.calculate_lj_force ? ForceType::LENNARD_JONES : ForceType::VERLET;

  std::unique_ptr<output::FileWriter> writer;
  if (params_.xyz_output) {
    writer = std::make_unique<output::XYZWriter>(particles);
  } else {
    writer = std::make_unique<output::VTKWriter>(particles);
  }

  while (current_time < params_.end_time) {

    Calculation<Position>::run(particles, params_.time_delta);
    Calculation<Force>::run(particles, FORCE_TYPE);
    Calculation<Velocity>::run(particles, params_.time_delta);

    if (iteration % 10 == 0 && params_.enable_output)
      writer->plot_particles(params_.output_path, iteration++);
    logger.trace("Iteration " + std::to_string(iteration) + " finished.");
    current_time += params_.time_delta;
  }
  logger.info("output written. Terminating...");

  logger.info("Number of particles: " + std::to_string(particles.size()));

  logger.info("Simulation finished.");
}