/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include "FileWriter.h"
#include "simulator/particle/ParticleContainer.h"

#include <fstream>
#include <list>

namespace output {

class XYZWriter : public FileWriter {


public:
  XYZWriter(ParticleContainer &particles);

  void plot_particles(const std::string &filename, int iteration) override;

  std::string out_name{"MD_xyz"};
};

} // namespace output
