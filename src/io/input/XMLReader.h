#pragma once

#include "simulator/particle/ParticleContainer.h"
#include <iostream>
#include <sstream>

struct SimulationParameters {
  double t_end;
  double delta_t;
  std::string output_basename;
  unsigned int write_frequency;
};

class XMLReader {
private:
public:
  XMLReader();
  ~XMLReader();

  void readXMLFile(ParticleContainer &particles,
                   SimulationParameters &simParams,
                   const std::string &filename);

  template <typename Container>
  std::string containerToString(const Container &container) {
    std::ostringstream oss;
    oss << "{ ";
    for (auto it = container.begin(); it != container.end(); ++it) {
      oss << *it;
      if (std::next(it) != container.end()) {
        oss << ", ";
      }
    }
    oss << " }";
    return oss.str();

    oss << " }";
    return oss.str();
  }
};
