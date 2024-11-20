/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include "particleSim/ParticleContainer.h"
#include "FileWriter.h"

#include <fstream>
#include <list>

namespace output {

class XYZWriter : FileWriter {

public:

  void plotParticles(ParticleContainer &particles, const std::string &filename,
                     int iteration);

  void write_file(const std::string &data, int iteration) override;

  void plot_particles() override;
};

} // namespace outputWriter
