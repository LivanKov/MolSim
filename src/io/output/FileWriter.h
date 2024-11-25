#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <memory>
#include <simulator/particle/ParticleContainer.h>
#include <string>

namespace output {

class FileWriter {
public:
  FileWriter(ParticleContainer &particles);

  virtual void plot_particles(const std::string &filepath, int iteration);

  ParticleContainer &particles;
};
} // namespace output
#endif