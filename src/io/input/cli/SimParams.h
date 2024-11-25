#include <string>

#pragma once

struct SimParams {
  std::string input_path;
  std::string output_path;
  double end_time;
  double time_delta;
  bool xyz_output;
  bool calculate_lj_force;
  bool enable_output;
  std::string log_level;
};