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
#include <iostream>
#include <memory>
#include <spdlog/spdlog.h>
#include <utility>
#include "io/input/XMLReader.h"
#include "particle/ParticleGenerator.h"

std::unique_ptr<Simulation> Simulation::generate_simulation(SimParams &params) {
  std::unique_ptr<Simulation> ptr = std::make_unique<Simulation>(params);
  return ptr;
}

Simulation::Simulation(SimParams &params) : params_(params) {}

ParticleContainer Simulation::readFile(char *argsv1, SimParams &params) {
  ParticleContainer particles{};
  // FileReader::readFile(particles, params_.input_path.data());
  XMLReader::readXMLFile(particles, params, argsv1);

  return particles;
}

void Simulation::run(ParticleContainer &particles) {

  Logger &logger = Logger::getInstance(params_.log_level);

  logger.warn("Starting a simulation with:");
  logger.info("\tEnd time: " + std::to_string(params_.end_time));
  logger.info("\tDelta: " + std::to_string(params_.time_delta));

  int iteration{0};
  double current_time{0};

  ForceType FORCE_TYPE = params_.calculate_grav_force
                             ? ForceType::VERLET
                             : ForceType::LENNARD_JONES;
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

    iteration++;
    if (iteration % params_.write_frequency == 0 && !params_.disable_output){
      writer->plot_particles(params_.output_path, iteration);
    }

    logger.info("Iteration " + std::to_string(iteration) + " finished.");
    current_time += params_.time_delta;
  }
  logger.info("output written. Terminating...");

  logger.info("Number of particles: " + std::to_string(particles.size()));

  logger.warn("Simulation finished.");
}

/*void Simulation::runDisc() {

  const std::array<double, 3> &center{10,10,0};
  const std::array<double, 3> &initialVelocity{0,0,0};
  size_t radius{15};
  double h{2^(1/6)};
  double mass{1.0};
  ParticleContainer particles{};

  ParticleGenerator::insertDisc(center,initialVelocity,radius,h,mass,particles);

  const std::array<double, 3> &center2{0,0,1};
  const std::array<double, 3> &initialVelocity2{0,0,0};
  size_t radius2{3};
  double h2{2^(1/6)};
  double mass2{1.0};

  ParticleGenerator::insertDisc(center2,initialVelocity2,radius2,h2,mass2,particles);


  std::unique_ptr<output::FileWriter> writer;
  writer = std::make_unique<output::VTKWriter>(particles);
  writer->plot_particles(params_.output_path,1);



}*/;