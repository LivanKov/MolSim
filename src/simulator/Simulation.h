#include "io/input/cli/SimParams.h"
#include "particle/ParticleContainer.h"
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
   * @param params SimParams reference, pass the overrided simulation parameters
   * by Command line arguments.
   * @return std::unique_ptr<Simulation> pointer to the simulation object.
   */
  static std::unique_ptr<Simulation> generate_simulation(SimParams &params);
  /**
   * @brief Run the simulation, contains all the necessary logic for the
   * simulation.
   * @param particles ParticleContainer reference, pass the particle container
   * initialized in readFile.
   */
  void run(ParticleContainer &particles);
  /**
   * @brief read file from XML input, initialize particles.
   * @param simParams SimParams reference, pass the initial simulation
   * parameters by XML input.
   * @return ParticleContainer initialized particle container.
   */
  static ParticleContainer readFile(SimParams &simParams);
};
