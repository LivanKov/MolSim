/*
 * XYZWriter.cpp
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#include "XYZWriter.h"
#include <iomanip>
#include <sstream>

namespace output {

XYZWriter::XYZWriter(LinkedCellContainer &particles) : FileWriter(particles) {}

void XYZWriter::plot_particles(const std::string &filename, int iteration) {

  std::ofstream file;
  std::stringstream strstr;
  strstr << filename << "/" << out_name << "_" << std::setfill('0')
         << std::setw(4) << iteration << ".xyz";

  file.open(strstr.str().c_str());
  file << particles.size() << std::endl;
  file << "Generated by MolSim. See http://openbabel.org/wiki/XYZ_(format) for "
          "file format doku."
       << std::endl;

  for (auto &p : particles.particles) {
    std::array<double, 3> x = p.getX();
    file << "Ar ";
    file.setf(std::ios_base::showpoint);

    for (auto &xi : x) {
      file << xi << " ";
    }

    file << std::endl;
  }

  file.close();
}
// namespace outputWriter
} // namespace output