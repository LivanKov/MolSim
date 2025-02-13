#pragma once

#include "simulator/particle/container/LinkedCellContainer.h"

namespace output {
/**
 * @class CheckpointWriter
 * @brief Class for writing checkpoint files.
 */

class CheckpointWriter {

public:
  /**
   * @brief Default constructor for CheckpointWriter.
   */
  CheckpointWriter();
  ~CheckpointWriter();

  /**
   * @brief Writes a checkpoint file from the given particle container.
   * @param particles A reference to the particle container to write from.
   * @param filename The name of the file to write to.
   * @param delta_t The delta_t value to write.
   * @param t_end The t_end value to write.
   */
  static void writeCheckpoint(LinkedCellContainer &particles,
                              const std::string &filename, double delta_t,
                              double t_end);

private:
};
} // namespace output