#include <iostream>
#include "simulator/Simulation.h"
#include "io/input/cli/SimParams.h"
#include "io/input/cli/CommandParser.h"



int main(int argc, char *argsv[]) {
  SimParams parameters = CommandParser::parse(argc,argsv);
  auto simulation = Simulation::generate_simulation(parameters);
  simulation->run();
  return 0;
}


