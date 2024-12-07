#pragma once

#include "cli/SimParams.h"
#include "simulator/particle/container/ParticleContainer.h"
#include <iostream>
#include <sstream>

class XMLReader {
private:
public:
  XMLReader();
  ~XMLReader();

  static void readXMLFile(ParticleContainer &particles, SimParams &simParams);

  static BoundaryCondition parseBoundaryCondition(const std::string &value);
};
