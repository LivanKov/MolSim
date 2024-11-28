#include "io/input/cli/CommandParser.h"
#include "io/input/cli/SimParams.h"
#include "simulator/Simulation.h"
#include <iostream>

int main(int argc, char *argsv[]) {
  SimParams parameters{};
  ParticleContainer particles = Simulation::readFile(argsv[1], parameters);
  SimParams overridedParams = CommandParser::parse(argc, argsv, parameters);
  auto simulation = Simulation::generate_simulation(overridedParams);
  simulation->run(particles);
  return 0;
}
