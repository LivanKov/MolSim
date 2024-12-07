/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "simulator/particle/Particle.h"
#include "simulator/particle/container/ParticleContainer.h"

#include <list>

#include "utils/logger/Logger.h"

class FileReader {

public:
  FileReader();
  virtual ~FileReader() = 0;

  static void readFile(ParticleContainer &particles, char *filename);
};
