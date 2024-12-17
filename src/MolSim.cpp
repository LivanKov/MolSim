#include "io/input/cli/CommandParser.h"
#include "io/input/cli/SimParams.h"
#include "simulator/Simulation.h"
#include "simulator/particle/container/LinkedCellContainer.h"
#include "utils/logger/Logger.h"
#include <iostream>

int main(int argc, char *argsv[]) {
  SimParams parameters{};
  CommandParser::parse(argc, argsv, parameters);
  Logger &logger = Logger::getInstance(parameters.log_level);
  logger.info("Log level is set." + parameters.log_level);
  LinkedCellContainer particles = Simulation::readFile(parameters);
  auto simulation = Simulation::generate_simulation(parameters);
  simulation->run(particles);
  return 0;
}
