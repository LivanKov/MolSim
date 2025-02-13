#pragma once

#include "cli/SimParams.h"
#include "simulator/particle/container/LinkedCellContainer.h"
#include <iostream>
#include <sstream>

namespace input {
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

  /**
   * @brief Destructor for XMLReader
   */
  ~XMLReader();

  /**
   * @brief Reads an XML file and populates the simulation parameters and
   * particle data
   * @param particles A reference to the particle container to populate
   * @param simParams A reference to the simulation parameters to populate
   */
  static void readXMLFile(LinkedCellContainer &particles, SimParams &simParams);

  /**
   * @brief Parses a boundary condition string into an enum value
   * @param value The string to parse
   * @return The parsed boundary condition
   */
  static BoundaryCondition parseBoundaryCondition(const std::string &value);
};
} // namespace input