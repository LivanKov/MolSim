#include "io/input/FileReader.h"
#include "io/output/VTKWriter.h"
#include "io/output/XYZWriter.h"
#include "io/output/FileWriter.h"
#include "simulator/particle/ParticleContainer.h"
#include "utils/ArrayUtils.h"
#include "simulator/particle/ParticleGenerator.h"
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <list>
#include <unordered_map>
#include <variant>
#include <memory>
#include "logger/Logger.h"
#include "utils/SimParams.h"
#include "utils/CommandParser.h"
#include "simulator/calculations/Calculation.h"
#include "simulator/calculations/Force.h"
#include "simulator/calculations/Position.h"
#include "simulator/calculations/Velocity.h"



int main(int argc, char *argsv[]) {

  SimParams parameters = CommandParser::parse(argc,argsv);
  ParticleContainer particles{};

  Logger &logger = Logger::getInstance(parameters.log_level);

  FileReader fileReader;
  fileReader.readFile(particles, parameters.input_path.data());

  int iteration = 0;
  double current_time = parameters.start_time;

  logger.warn("Starting a simulation with:");
  logger.info("\tStart time: " + std::to_string(parameters.start_time));
  logger.info("\tEnd time: " + std::to_string(parameters.end_time));
  logger.info("\tDelta: " + std::to_string(parameters.time_delta));

  output::VTKWriter writer(particles);

  while (current_time < parameters.end_time) {
    
    Calculation<Position>::run(particles, parameters.time_delta);
    Calculation<Force>::run(particles,Type::VERLET);
    Calculation<Velocity>::run(particles, parameters.time_delta);

    iteration++;
    if (parameters.sparse_output && iteration % 10 == 0 && parameters.enable_output)
      writer.plot_particles(parameters.output_path ,iteration);
    else if (!parameters.sparse_output && parameters.enable_output)
      writer.plot_particles(parameters.output_path ,iteration);
    logger.trace("Iteration " + std::to_string(iteration) + " finished.");
    current_time += parameters.time_delta;
  }

  logger.info("output written. Terminating...");

  logger.debug("Number of particles: " + std::to_string(particles.size()));

  logger.warn("Simulation finished.");

  return 0;
}


