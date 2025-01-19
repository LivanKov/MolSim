#include "SimParams.h"

// Define and initialize static members
std::array<double, 3> SimParams::lower_left_corner = {0.0, 0.0, 0.0};
bool SimParams::fixed_Domain = false;
bool SimParams::enable_gravity = false;
bool SimParams::enable_thermo = false;
double SimParams::gravity = 0.0;
bool SimParams::enable_v_threshold = false;
double SimParams::v_threshold = 0.0;
