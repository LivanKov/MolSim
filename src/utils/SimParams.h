#include <string>


struct SimParams{
    std::string input_path;
    std::string output_path;

    double start_time;
    double end_time;
    double time_delta;

    bool sparse_output;
    bool xyz_output;
    bool calculate_lj_force;
    bool enable_output;
    std::string log_level;
};