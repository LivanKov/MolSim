#include <unistd.h>

namespace CommandParser{
    void parse(int argc, char**argv, SimParams params){
        if (argc < 2) {
    print_help();
    return 1;
  }

  int opt;

  while ((opt = getopt(argc, argsv, "e:d:i:o:thxl:fn")) != -1) {
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
    case 'o':
      parameters.output_path = std::string(optarg);
      break;
    case 't':
      parameters.sparse_output = false;
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
      parameters.calculate_lj_force = false;
      break;
    case 'n':
      parameters.enable_output = false;
      break;
    default:
      fprintf(stderr, "Usage: %s [-h] help\n", argsv[0]);
      return 1;
    }
  }
    }

    void print_help() {
  std::cout << "Usage: MolSim [options]\n";
  std::cout << "Options:\n";
  std::cout << "  -h                 Show this help message\n";
  std::cout << "  -o   <file_path>   Specify output file path\n";
  std::cout << "  -i   <file_path>   Specify input file path\n";
  std::cout
      << "  -e   <end_time>    Specify how long the simulation should run\n";
  std::cout << "  -d   <time_delta>  Specify time increments\n";
  std::cout << "  -t                 Enable testing mode (Writes a file for "
               "each iteration)\n";
  std::cout << "  -x                 Output .xyz files instead of .vpu\n";
  std::cout << "  -l   <log_level>   Option to choose the logging level "
               "[trace, debug, info, warn, error, off]\n";
  std::cout << "  -f                 Calculate Gravitational Force instead of "
               "Lennard-Jones Force\n";
  std::cout << "  -n   Disable all file output for the sake of performance\n";
    }
}