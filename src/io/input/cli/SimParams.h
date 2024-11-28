#include <string>

#pragma once

struct SimParams {
  std::string input_path;
  std::string output_path;
  double end_time;
  double time_delta;
  unsigned int write_frequency;
  bool xyz_output;
  bool calculate_lj_force;
  bool disable_output;
  std::string log_level;
};