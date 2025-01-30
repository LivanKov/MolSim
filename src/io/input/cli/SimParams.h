#include "simulator/particle/container/LinkedCellContainer.h"
#include <array>
#include <string>
#include <unordered_set>
#pragma once

struct SimParams {
  std::string output_path;
  std::string input_path;
  double end_time;
  double time_delta;
  double r_cutoff_radius;
  unsigned int write_frequency;
  bool xyz_output;
  bool calculate_grav_force;
  bool disable_output;
  std::string log_level;
  std::array<double, 3> domain_size;
  bool linked_cells;
  DomainBoundaryConditions boundaryConditions;
  bool reflective;
  bool periodic;
  static double gravity;
  static bool enable_gravity;
  static double z_gravity;
  static bool enable_z_gravity;
  static std::array<double, 3> lower_left_corner;
  static bool fixed_Domain;

  // Thermostats
  static bool enable_thermo;
  double initial_temp;
  double target_temp;
  double delta_temp;
  bool is_gradual;
  unsigned int n_thermostats;
  size_t dimensions;
  bool enable_brownian;

  // Checkpoint
  bool resume_from_checkpoint;
  bool checkpoint_only;
  double resume_start_time;

  // velocity threshold
  static bool enable_v_threshold;
  static double v_threshold;

  // additional force (FZup)
  static double additional_force_zup;
  static double additional_force_time_limit;
  static bool apply_fzup;

  // membrane parameters
  static double membrane_stiffness;
  static double membrane_bond_length;

  //precompute parameters
  static bool precompute_sigma;
  static bool precompute_epsilon;

  static bool is_membrane;
};
