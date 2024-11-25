#include <iostream>
#include "simulator/Simulation.h"
#include "utils/SimParams.h"
#include "utils/CommandParser.h"



int main(int argc, char *argsv[]) {
  SimParams parameters = CommandParser::parse(argc,argsv);
  auto simulation = Simulation::generate_simulation(parameters);
  simulation->run();
  return 0;
}


