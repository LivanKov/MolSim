#include "SimParams.h"

// Define and initialize static members
std::array<double, 3> SimParams::lower_left_corner = {0.0, 0.0, 0.0};
bool SimParams::fixed_Domain = false;
bool SimParams::enable_gravity = false;
bool SimParams::enable_thermo = false;
double SimParams::gravity = 0.0;
bool SimParams::enable_v_threshold = false;
double SimParams::v_threshold = 0.0;
std::unordered_set<int> SimParams::additional_force_particle_ids = {};
double SimParams::additional_force_z_gravity = 0.0;
double SimParams::additional_force_time_limit = 0.0;
bool SimParams::enable_additional_force = false;
double SimParams::membrane_stiffness = 0.0;
double SimParams::membrane_bond_length = 0.0;