#include "Simulation.h"
#include "io/input/FileReader.h"
#include "io/input/XMLReader.h"
#include "io/output/FileWriter.h"
#include "io/output/VTKWriter.h"
#include "io/output/XYZWriter.h"
#include "particle/ParticleGenerator.h"
#include "particle/container/DirectSumContainer.h"
#include "simulator/calculations/Calculation.h"
#include "simulator/calculations/Force.h"
#include "simulator/calculations/Position.h"
#include "simulator/calculations/Velocity.h"
#include "utils/logger/Logger.h"
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
  LinkedCellContainer particles{{180.0, 90.0}, 3.0};
  XMLReader::readXMLFile(particles, params);
  if (params.reflective) {
    particles.reflective_flag = true;
    particles.boundary_conditions_ = DomainBoundaryConditions{
        BoundaryCondition::Reflecting, BoundaryCondition::Reflecting,
        BoundaryCondition::Reflecting, BoundaryCondition::Reflecting,
        BoundaryCondition::Reflecting, BoundaryCondition::Reflecting};
  }
  if (params.periodic) {
    particles.reflective_flag = false;
    particles.periodic_flag = true;
    particles.boundary_conditions_ = DomainBoundaryConditions{
        BoundaryCondition::Periodic, BoundaryCondition::Periodic,
        BoundaryCondition::Periodic, BoundaryCondition::Periodic,
        BoundaryCondition::Periodic, BoundaryCondition::Periodic};
  }
  return particles;
}

void Simulation::run(LinkedCellContainer &particles) {

  Logger &logger = Logger::getInstance(params_.log_level);

  logger.warn("Starting a simulation with:");
  logger.info("\tEnd time: " + std::to_string(params_.end_time));
  logger.info("\tDelta: " + std::to_string(params_.time_delta));

  int iteration{0};
  double current_time{0};

  ForceType FORCE_TYPE = params_.calculate_grav_force
                             ? ForceType::GRAVITATIONAL
                             : ForceType::LENNARD_JONES;
  std::unique_ptr<output::FileWriter> writer;
  if (params_.xyz_output) {
    writer = std::make_unique<output::XYZWriter>(particles);
  } else {
    writer = std::make_unique<output::VTKWriter>(particles);
  }

  OPTIONS option =
      params_.linked_cells ? OPTIONS::LINKED_CELLS : OPTIONS::DIRECT_SUM;

  while (current_time < params_.end_time) {

    Calculation<Position>::run(particles, params_.time_delta, option);
    Calculation<Force>::run(particles, FORCE_TYPE, option);
    Calculation<Velocity>::run(particles, params_.time_delta);

    iteration++;
    if (iteration % params_.write_frequency == 0 && !params_.disable_output) {
      writer->plot_particles(params_.output_path, iteration);
    }

    logger.info("Iteration " + std::to_string(iteration) + " finished.");
    current_time += params_.time_delta;
  }
  logger.info("output written. Terminating...");

  logger.info("Number of particles: " + std::to_string(particles.size()));

  logger.info("Particles left the domain: " +
              std::to_string(particles.particles_left_domain));

  logger.info("Amount of halo particles:" +
              std::to_string(particles.halo_count));

  logger.warn("Simulation finished.");
};