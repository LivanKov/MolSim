#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <memory>
#include "simulator/particle/container/ParticleContainer.h"
#include <string>

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
  FileWriter(ParticleContainer &particles);

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
  ParticleContainer &particles;
};
} // namespace output
#endif