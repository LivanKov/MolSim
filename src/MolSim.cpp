#include "io/input/FileReader.h"
#include "io/output/VTKWriter.h"
#include "io/output/XYZWriter.h"
#include "io/output/FileWriter.h"
#include "simulator/particle/ParticleContainer.h"
#include "utils/ArrayUtils.h"
#include "simulator/particle/ParticleGenerator.h"
#include <spdlog/spdlog.h>
#include <iostream>
#include "utils/logger/Logger.h"
#include "utils/SimParams.h"
#include "utils/CommandParser.h"
#include "simulator/calculations/Calculation.h"
#include "simulator/calculations/Force.h"
#include "simulator/calculations/Position.h"
#include "simulator/calculations/Velocity.h"
#include "simulator/Simulation.h"



int main(int argc, char *argsv[]) {
  SimParams parameters = CommandParser::parse(argc,argsv);
  auto simulation = Simulation::generate_simulation(parameters);
  simulation->run();
  return 0;
}


