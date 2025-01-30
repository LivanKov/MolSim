#include "SimParams.h"
#include "simulator/calculations/Calculation.h"
#include <iostream>
#include <unistd.h>

#pragma once

namespace CommandParser {

void print_help();

SimParams &parse(int argc, char **argv, SimParams &parameters);
} // namespace CommandParser