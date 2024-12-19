#include "SimParams.h"

// Define and initialize static members
std::array<double, 3> SimParams::lower_left_corner = {0.0, 0.0, 0.0};
bool SimParams::fixed_Domain = false;
bool SimParams::gravity_applied = false;
bool SimParams::enable_thermo = false;
