#pragma once

#include "cli/SimParams.h"
#include "simulator/particle/container/DirectSumContainer.h"
#include <iostream>
#include <sstream>

class XMLReader {
private:
public:
  XMLReader();
  ~XMLReader();

  static void readXMLFile(DirectSumContainer &particles, SimParams &simParams);

  static BoundaryCondition parseBoundaryCondition(const std::string &value);
};
