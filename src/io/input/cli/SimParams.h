#include "simulator/calculations/Force.h"
#include "simulator/particle/container/LinkedCellContainer.h"
#include <array>
#include <string>
#include <unordered_set>
#pragma once

/**
 * @struct SimParams
 * @brief Stores simulation parameters.
 */
struct SimParams {

  /**
   * @brief The path to the output directory.
   */
  std::string output_path;
  /**
   * @brief The path to the file, from which the input is read
   */
  std::string input_path;
  /**
   * @brief The end time of the simulation
   */
  double end_time;
  /**
   * @brief The time step of the simulation
   */
  double time_delta;
  /**
   * @brief The cutoff radius for the cell container
   */
  double r_cutoff_radius;
  /**
   * @brief The frequency at which the output is written
   */
  unsigned int write_frequency;
  /**
   * @brief Flag to indicate whether to output XYZ files
   */
  bool xyz_output;
  /**
   * @brief Flag to indicate whether to calculate the gravitational force
   */
  bool calculate_grav_force;
  /**
   * @brief Flag to indicate whether to disable output
   */
  bool disable_output;
  /**
   * @brief Control the options of the logger
   */
  std::string log_level;
  /*
   * @brief The size of the domain
   */
  std::array<double, 3> domain_size;
  /*
   * @brief Flag to indicate whether to use linked cells
   */
  bool linked_cells;
  /*
   * @brief Control the boundary conditions of the domain
   */
  DomainBoundaryConditions boundaryConditions;
  /*
   * @brief Flag to indicate whether to use reflective boundary conditions
   */
  bool reflective;
  /*
   * @brief Flag to indicate whether to use periodic boundary conditions
   */
  bool periodic;

  /*
   * @brief The gravitational constant. Global variable
   */
  static double gravity;

  /*
   * @brief Flag to indicate whether to calculate the gravitational force.
   * Global variable
   */
  static bool enable_gravity;

  /*
   * @brief The gravitational constant on the z-axis. Global variable
   */
  static double z_gravity;

  /*
   * @brief Flag to indicate whether to calculate the gravitational force on the
   * z-axis. Global variable
   */
  static bool enable_z_gravity;

  /*
   * @brief The lower left corner of the domain. Global variable
   */
  static std::array<double, 3> lower_left_corner;

  /*
   * @brief Flag to indicate whether the domain is fixed. Global variable
   */
  static bool fixed_Domain;

  // Thermostat variables
  /*
   * @brief Flag to indicate whether to use thermostats. Global variable
   */
  static bool enable_thermo;
  /*
   * @brief Initial temperature of the system.
   */
  double initial_temp;
  /*
   * @brief Target temperature of the system.
   */
  double target_temp;
  /*
   * @brief Temperature change per time step.
   */
  double delta_temp;
  /*
   * @brief Flag to indicate whether to apply gradual temperature change.
   */
  bool is_gradual;
  /*
   * @brief Number of applications of the thermostat.
   */
  unsigned int n_thermostats;
  /*
   * @brief Number of dimensions.
   */
  size_t dimensions;
  /*
   * @brief Flag to indicate whether to enable Brownian motion.
   */
  bool enable_brownian;

  // Checkpoint variables
  /*
   * @brief Flag to indicate whether to resume simulation from checkpoint.
   */
  bool resume_from_checkpoint;
  /*
   * @brief Flag to indicate whether to only execute from checkpoint.
   */
  bool checkpoint_only;
  /*
   * @brief Time at which to resume the simulation.
   */
  double resume_start_time;

  // openmp variables
  /**
   * @brief Flag to indicate whether to use OpenMP parallelization.Global
   * Variable
   */
  static bool enable_omp;
  /**
   * @brief Choose a strategy for OpenMP parallelization.Global Variable
   */
  static OMPSTRATEGY ompstrategy;

  // additional force (FZup)
  /**
   * @brief Flag to indicate whether to apply an additional force.Global
   * Variable
   */
  static bool apply_fzup;
  /**
   * @brief The strength of the additional force.Global Variable
   */
  static double additional_force_zup;
  /**
   * @brief The time limit for the additional force.Global Variable
   */
  static double additional_force_time_limit;

  // membrane parameters
  /**
   * @brief The stiffness of the membrane.Global Variable
   */
  static double membrane_stiffness;
  /**
   * @brief The bond length of the membrane.Global Variable
   */
  static double membrane_bond_length;
  /**
   * @brief Flag to indicate whether the simulation is a membrane
   * simulation.Global Variable
   */
  static bool is_membrane;
};
