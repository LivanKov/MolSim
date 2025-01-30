#include "SimParams.h"
#include "simulator/calculations/Calculation.h"
#include <iostream>
#include <unistd.h>

#pragma once

namespace CommandParser {

void print_help();

SimParams &parse(int argc, char **argv, SimParams &parameters);

  while ((opt = getopt(argc, argv, "e:d:i:t:o:hxl:fnurcv:pa")) != -1) {
    switch (opt) {
    case 'e':
      parameters.end_time = atof(optarg);
      break;
    case 'd':
      parameters.time_delta = atof(optarg);
      break;
    case 'i':
      parameters.input_path = std::string(optarg);
      break;
    case 't':
      parameters.write_frequency = atoi(optarg);
      break;
    case 'o':
      parameters.output_path = std::string(optarg);
      break;
    case 'h':
      print_help();
      break;
    case 'x':
      parameters.xyz_output = true;
      break;
    case 'l':
      parameters.log_level = std::string(optarg);
      break;
    case 'f':
      parameters.calculate_grav_force = true;
      break;
    case 'n':
      parameters.disable_output = true;
      break;
    case 'u':
      parameters.linked_cells = true;
      break;
    case 'r':
      parameters.resume_from_checkpoint = true;
      break;
    case 'c':
      parameters.checkpoint_only = true;
      break;
    case 'v':
      SimParams::enable_v_threshold = true;
      SimParams::v_threshold = atof(optarg);
      break;
    case 'p':
      SimParams::enable_omp = true;
      SimParams::ompstrategy = OMPSTRATEGY::FORK_JOIN;
      break;
    case 'a':
      SimParams::enable_omp = true;
      SimParams::ompstrategy = OMPSTRATEGY::TASKING;
      break;
    default:
      fprintf(stderr, "Usage: %s [-h] help\n", argv[0]);
    }
  }
  return parameters;
}
} // namespace CommandParser