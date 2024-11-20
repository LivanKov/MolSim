/*
 * VTKWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include "outputUtils/vtk-unstructured.h"
#include "particleSim/Particle.h"
#include "io/output/FileWriter.h"

#include <list>

namespace output {

/**
 * This class implements the functionality to generate vtk output from
 * particles.
 */
class VTKWriter : public FileWriter{

public:
  /**
   * set up internal data structures and prepare to plot a particle.
   */
  VTKWriter(std::shared_ptr<ParticleContainer>& particles);

  /**
   * plot type, mass, position, velocity and force of a particle.
   *
   * @note: initializeOutput() must have been called before.
   */
  void plotParticle(Particle &p);

  /**
   * writes the final output file.
   *
   * @param filename the base name of the file to be written.
   * @param iteration the number of the current iteration,
   *        which is used to generate an unique filename
   */
  void write_file(const std::string &data, int iteration) override;

  void plot_particles() override;

private:
  
  VTKFile_t *vtkFile;
};

} // namespace outputWriter
