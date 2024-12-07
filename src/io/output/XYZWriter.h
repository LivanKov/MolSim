/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include "FileWriter.h"
#include "simulator/particle/container/DirectSumContainer.h"

#include <fstream>
#include <list>

namespace output {

class XYZWriter : public FileWriter {

public:
  XYZWriter(DirectSumContainer &particles);

  void plot_particles(const std::string &filename, int iteration) override;

  std::string out_name{"MD_xyz"};
};

} // namespace output
