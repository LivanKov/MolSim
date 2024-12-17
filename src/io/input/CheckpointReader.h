#pragma once

#include "simulator/particle/container/LinkedCellContainer.h"

class CheckpointReader {

public:
  CheckpointReader();
  ~CheckpointReader();

  static void readCheckpoint(LinkedCellContainer &particles, double &delta_t,
                             double &t_end);

private:
};
