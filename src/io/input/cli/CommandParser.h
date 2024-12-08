#include "SimParams.h"
#include <iostream>
#include <unistd.h>

namespace CommandParser {

void print_help() {
  std::cout << "Usage: MolSim [options]\n";
  std::cout << "Options:\n";
  std::cout << "  -h                  Show this help message\n";
  std::cout << "  -o   <file_path>    Specify output file path (including "
               "filename)\n";
  std::cout << "  -i   <file_path>   Specify input file path\n";
  std::cout
      << "  -e   <end_time>    Specify how long the simulation should run\n";
  std::cout << "  -d   <time_delta>   Specify time increments\n";
  std::cout << "  -t   <write_freq>   Write a file for every -t iteration\n";

  std::cout << "  -x                  Output .xyz files instead of .vpu\n";
  std::cout << "  -l   <log_level>    Option to choose the logging level "
               "[trace, debug, info, warn, error, off]\n";
  std::cout
      << "  -f                   Calculate Gravitational Force instead of "
         "Lennard-Jones Force\n";
  std::cout << "  -n                  Disable all file output for the sake of "
               "performance\n";
}

SimParams parse(int argc, char **argv, SimParams &parameters) {
  if (argc < 2) {
    print_help();
  }

  int opt;

  while ((opt = getopt(argc, argv, "e:d:i:t:o:hxl:fn")) != -1) {
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
    default:
      fprintf(stderr, "Usage: %s [-h] help\n", argv[0]);
    }
  }
  return parameters;
}
} // namespace CommandParser