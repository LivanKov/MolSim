#include "io/input/cli/SimParams.h"
#include <memory>

#pragma once

/**
 * @class Simulation
 * @brief Main simulation class, provides an interface to run the simulation.
 */
class Simulation {
private:
  SimParams &params_;

public:
  /** 
   * @brief Constructor
   * @param params SimParams reference, pass the simulatin parameters.
   */
  Simulation(SimParams &params);
  /** 
   * @brief Generate a simulation object.
   * @param params SimParams reference, pass the simulatin parameters.
   * @return std::unique_ptr<Simulation> pointer to the simulation object.
   */
  static std::unique_ptr<Simulation> generate_simulation(SimParams &params);
  /** 
   * @brief Run the simulation, contains all the necessary logic for the simulation.
   */
  void run();
};
