namespace CommandParser{
    void parse(int argc, char**argv, SimParams params){

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