#include "SimParams.h"

#pragma once

namespace CommandParser {

void print_help();

SimParams& parse(int argc, char **argv, SimParams &parameters);

} // namespace CommandParser