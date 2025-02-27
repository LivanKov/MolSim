#include "simulator/particle/container/LinkedCellContainer.h"
#include <memory>
#include <string>

#pragma once

namespace output {

/**
 * @class FileWriter
 * @brief Abstract file writer, used to declare writers with different
 * properties.
 */
class FileWriter {
public:
  /**
   * @brief Constructor, should be overriden by the derived classes.
   * @param particles ParticleContainer reference.
   */
  FileWriter(LinkedCellContainer &particles);

  /**
   * @brief Virtual function to plot particles, should be overriden by the
   * derived classes.
   * @param filepath std::string reference, path to the file.
   * @param iteration int, iteration number.
   */
  virtual void plot_particles(const std::string &filepath, int iteration);

  /**
   * @brief ParticleContainer reference
   */
  LinkedCellContainer &particles;
};
} // namespace output