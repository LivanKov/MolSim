#pragma once

#include "simulator/particle/container/LinkedCellContainer.h"

class CheckpointWriter {

public:
  CheckpointWriter();
  ~CheckpointWriter();

  static void writeCheckpoint(LinkedCellContainer &particles,
                              const std::string &filename, double delta_t,
                              double t_end);

private:
};
