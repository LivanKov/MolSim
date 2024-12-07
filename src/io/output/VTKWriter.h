/*
 * VTKWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include "io/output/FileWriter.h"
#include "outputUtils/vtk-unstructured.h"
#include "simulator/particle/Particle.h"

#include <list>

namespace output {

/**
 * This class implements the functionality to generate vtk output from
 * particles.
 */
class VTKWriter : public FileWriter {

public:
  /**
   * set up internal data structures and prepare to plot a particle.
   */
  VTKWriter(DirectSumContainer &particles);

  /**
   * writes the final output file.
   *
   * @param filename the base name of the file to be written.
   * @param iteration the number of the current iteration,
   *        which is used to generate an unique filename
   */

  void plot_particles(const std::string &filename, int iteration) override;

private:
  /**
   * plot type, mass, position, velocity and force of a particle.
   *
   * @note: initializeOutput() must have been called before.
   */

  void plotParticle(Particle &p);

  VTKFile_t *vtkFile;

  std::string out_name{"MD_vtk"};
};

} // namespace output
