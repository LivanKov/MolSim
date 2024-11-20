/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include "simulator/ParticleContainer.h"
#include "FileWriter.h"

#include <fstream>
#include <list>

namespace output {

class XYZWriter : FileWriter {

public:

  XYZWriter(std::shared_ptr<ParticleContainer>& particles);

  void plot_particles(const std::string &filename,
                     int iteration) override;
};

} // namespace outputWriter
