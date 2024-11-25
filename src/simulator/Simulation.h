#include "io/input/cli/SimParams.h"
#include <memory>

#pragma once

// Simulation.h
class Simulation {
private:
  SimParams &params_;

public:
  Simulation(SimParams &params);
  static std::unique_ptr<Simulation> generate_simulation(SimParams &params);
  void run();
};
