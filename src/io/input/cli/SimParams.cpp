#include "SimParams.h"
#include "simulator/calculations/Force.h"

// Define and initialize static members
std::array<double, 3> SimParams::lower_left_corner = {0.0, 0.0, 0.0};
bool SimParams::fixed_Domain = false;
bool SimParams::enable_gravity = false;
bool SimParams::enable_thermo = false;
double SimParams::gravity = 0.0;
double SimParams::z_gravity = 0.0;
bool SimParams::enable_z_gravity = false;
bool SimParams::enable_v_threshold = false;
double SimParams::v_threshold = 0.0;
bool SimParams::enable_omp = false;
OMPSTRATEGY SimParams::ompstrategy = OMPSTRATEGY::FORK_JOIN;
bool SimParams::apply_fzup = false;
double SimParams::additional_force_zup = 0.0;
double SimParams::additional_force_time_limit = 0.0;
double SimParams::membrane_stiffness = 0.0;
double SimParams::membrane_bond_length = 0.0;
bool SimParams::precompute_epsilon = false;
bool SimParams::precompute_sigma = false;
bool SimParams::is_membrane = false;

