/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "simulator/particle/Particle.h"
#include "simulator/particle/ParticleContainer.h"

#include <list>

#include "logger/Logger.h"

class FileReader {

public:
  FileReader();
  virtual ~FileReader();

  void readFile(ParticleContainer &particles, char *filename);
};