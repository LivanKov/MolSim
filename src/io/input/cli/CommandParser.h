#include "SimParams.h"
#include "simulator/calculations/Calculation.h"
#include <iostream>
#include <unistd.h>

#pragma once

namespace CommandParser {

/*
 * @brief Print input help message.
 */
void print_help();

/*
 * @brief Parse input arguments.
 */
SimParams &parse(int argc, char **argv, SimParams &parameters);
} // namespace CommandParser