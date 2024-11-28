#include "io/input/cli/SimParams.h"
#include <memory>
#include "particle/ParticleContainer.h"

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
  void run(ParticleContainer &particles);
  /**
   * @brief read file from XML input, initialize particles.
   * @return ParticleContainer initialized particle container.
   */
  static ParticleContainer readFile(char *argv1, SimParams &SimParams);
};
