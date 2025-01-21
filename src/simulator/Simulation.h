#include "Thermostat.h"
#include "io/input/cli/SimParams.h"
#include "io/output/FileWriter.h"
#include "particle/container/DirectSumContainer.h"
#include "particle/container/LinkedCellContainer.h"
#include "simulator/calculations/Calculation.h"
#include "simulator/calculations/Force.h"
#include <chrono>
#include <memory>

#pragma once

/**
 * @brief Global variables for particle tracking.
 * @particles_left_domain: Number of particles left in the domain.
 * @particle_id: Particle ID.
 */

/**
 * @class Simulation
 * @brief Main simulation class, provides an interface to run the simulation.
 */

class Simulation {
private:

  /**
   * @brief Creates the appropriate file writer (e.g., XYZWriter or VTKWriter)
   * based on the parameters.
   * @param particles A reference to the particle container.
   * @return A unique pointer to the created FileWriter instance.
   */
  std::unique_ptr<output::FileWriter>
  createFileWriter(LinkedCellContainer &particles) const;

  /**
   * @brief Handles simulation execution in checkpoint-only mode.
   * @param particles A reference to the particle container.
   * @param current_time The current simulation time.
   * @param option The simulation method option (e.g., DIRECT_SUM or
   * LINKED_CELLS).
   * @param force_type The type of force calculation (e.g., GRAVITATIONAL or
   * LENARD_JONES).
   */
  void checkpointMode(LinkedCellContainer &particles, double &current_time,
                      OPTIONS option, ForceType force_type);

  /**
   * @brief Handles the main simulation loop.
   * @param particles A reference to the particle container.
   * @param current_time The current simulation time.
   * @param iteration The current iteration number.
   * @param total_molecule_updates The total number of molecule updates during
   * the simulation.
   * @param writer A unique pointer to the file writer for output.
   * @param thermostat A reference to the Thermostat instance.
   * @param option The simulation method option (e.g., DIRECT_SUM or
   * LINKED_CELLS).
   * @param force_type The type of force calculation (e.g., GRAVITATIONAL or
   * LENARD_JONES).
   */
  void simulate(LinkedCellContainer &particles, double &current_time,
                int &iteration, size_t &total_molecule_updates,
                std::unique_ptr<output::FileWriter> &writer,
                Thermostat &thermostat, OPTIONS option, ForceType force_type);

  /**
   * @brief Logs simulation performance metrics such as runtime and updates per
   * second.
   * @param runtime The total runtime of the simulation.
   * @param total_molecule_updates The total number of molecule updates during
   * the simulation.
   */
  void logPerformance(const std::chrono::duration<double> &runtime,
                      size_t total_molecule_updates) const;

  /// Simulation parameters
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
   * @param particles DirectSumContainer reference, pass the particle container
   * initialized in readFile.
   */
  void run(LinkedCellContainer &particles);

  /**
   * @brief read file from XML input, initialize particles.
   * @param simParams SimParams reference, pass the initial simulation
   * parameters by XML input.
   * @return DirectSumContainer initialized particle container.
   */
  static LinkedCellContainer readFile(SimParams &simParams);
};
