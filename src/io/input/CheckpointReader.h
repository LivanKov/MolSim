#pragma once

#include "simulator/particle/container/LinkedCellContainer.h"

namespace input {

/**
 * @class CheckpointReader
 * @brief Class for reading checkpoint files.
 */

class CheckpointReader {

public:
  /**
   * @brief Default constructor for CheckpointReader.
   */
  CheckpointReader();
  ~CheckpointReader();

  /**
   * @brief Reads a checkpoint file and populates the particle container with
   * particles.
   * @param particles A reference to the particle container to populate.
   * @param delta_t A reference to the delta_t variable to populate.
   * @param t_end A reference to the t_end variable to populate.
   */

  static void readCheckpoint(LinkedCellContainer &particles, double &delta_t,
                             double &t_end);

private:
};

} // namespace input
