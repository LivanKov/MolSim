/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include "FileWriter.h"
#include "simulator/particle/container/LinkedCellContainer.h"

#include <fstream>
#include <list>

namespace output {

class XYZWriter : public FileWriter {

public:
  /**
   * @brief Constructor.
   */

  XYZWriter(LinkedCellContainer &particles);

  /**
   * @brief Write the current state of the particles to a xyz-file.
   * @param filename Name of the file to write to.
   * @param iteration Current iteration of the simulation.
   */
  void plot_particles(const std::string &filename, int iteration) override;

  /**
   * @brief Default prefix for the output files
   */
  std::string out_name{"MD_xyz"};
};

} // namespace output
