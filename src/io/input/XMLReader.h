#pragma once

#include "cli/SimParams.h"
#include "simulator/particle/container/LinkedCellContainer.h"
#include <iostream>
#include <sstream>

/**
 * @class XMLReader
 * @brief A class for reading XML files and parsing simulation parameters
 *
 * This class provides functionality to read XML configuration files
 * and parse the contents into simulation parameters and particle data.
 */

class XMLReader {
private:
public:
  /**
   * @brief Default constructor for XMLReader
   */
  XMLReader();
  ~XMLReader();

  static void readXMLFile(LinkedCellContainer &particles, SimParams &simParams);

  static BoundaryCondition parseBoundaryCondition(const std::string &value);
};
