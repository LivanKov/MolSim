#pragma once

#include "cli/SimParams.h"
#include "simulator/particle/ParticleContainer.h"
#include <iostream>
#include <sstream>

class XMLReader {
private:
public:
  XMLReader();
  ~XMLReader();

  static void readXMLFile(ParticleContainer &particles, SimParams &simParams);
};
