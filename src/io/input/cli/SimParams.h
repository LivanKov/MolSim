#include <string>
#include <array>
#pragma once

struct SimParams {
  std::string output_path;
  double end_time;
  double time_delta;
  double r_cutoff_radius;
  unsigned int write_frequency;
  bool xyz_output;
  bool calculate_grav_force;
  bool disable_output;
  std::string log_level;
  std::array<double, 3> domain_size;
};