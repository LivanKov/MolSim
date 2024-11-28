#pragma once

#include "simulator/particle/ParticleContainer.h"
#include <iostream>
#include <sstream>
#include "cli/SimParams.h"

class XMLReader {
private:
public:
  XMLReader();
  ~XMLReader();

  static void readXMLFile(ParticleContainer &particles,
                   SimParams &simParams,
                   const std::string &filename);
};
