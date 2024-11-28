#include "io/input/cli/CommandParser.h"
#include "io/input/cli/SimParams.h"
#include "simulator/Simulation.h"
#include <iostream>

int main(int argc, char *argsv[]) {
  SimParams parameters{};
  ParticleContainer particles = Simulation::readFile(argsv[1], parameters);
  parameters = CommandParser::parse(argc, argsv);
  auto simulation = Simulation::generate_simulation(parameters);
  simulation->run(particles);
  return 0;
}
