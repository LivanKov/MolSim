/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "simulator/particle/Particle.h"
#include "simulator/particle/container/LinkedCellContainer.h"

#include <list>

#include "utils/logger/Logger.h"

namespace input {

/**
 * @class FileReader
 * @brief Reads particle data from a file.
 */

class FileReader {

public:
  /**
   * @brief Default constructor for FileReader.
   */
  FileReader();
  virtual ~FileReader() = 0;

  /**
   * @brief Reads particle data from a file.
   * @param particles A reference to the particle container to populate.
   * @param filename The name of the file to read from.
   */
  static void readFile(LinkedCellContainer &particles, char *filename);
};
} // namespace input