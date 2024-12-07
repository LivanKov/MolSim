/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"
#include "simulator/particle/container/DirectSumContainer.h"
#include "simulator/particle/ParticleGenerator.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

#include "utils/logger/Logger.h"

FileReader::FileReader() = default;

FileReader::~FileReader() = default;

auto containerToString = [](const auto &container) {
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
};

void FileReader::readFile(DirectSumContainer &particles, char *filename) {
  std::array<double, 3> x;
  std::array<size_t, 3> mesh;
  double d;
  double m;
  std::array<double, 3> v;
  size_t num_objects = 0;

  Logger &logger = Logger::getInstance();

  std::ifstream input_file(filename);
  std::string tmp_string;

  if (input_file.is_open()) {

    getline(input_file, tmp_string);
    logger.debug("Read line: " + tmp_string);

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      logger.debug("Read line: " + tmp_string);
    }
    std::istringstream numstream(tmp_string);

    if (numstream.str().size() == 1) {

      numstream >> num_objects;

      logger.info("    Reading " + num_objects);
      getline(input_file, tmp_string);
      logger.debug("Read line: " + tmp_string);

      for (size_t i = 0; i < num_objects; i++) {
        std::istringstream datastream(tmp_string);

        for (auto &xj : x) {
          datastream >> xj;
        }
        for (auto &vj : v) {
          datastream >> vj;
        }
        if (datastream.eof()) {
          logger.error("Error reading file: eof reached unexpectedly reading "
                       "from line " +
                       std::to_string(i));

          exit(-1);
        }
        datastream >> m;
        particles.emplace_back(x, v, m);

        getline(input_file, tmp_string);
        logger.info("Read line: " + tmp_string);
      }
    } else {
      logger.info("    Creating a cuboid ");

      getline(input_file, tmp_string);
      std::istringstream numstream(tmp_string);

      numstream >> num_objects;

      logger.info("    Reading " + num_objects);

      getline(input_file, tmp_string);

      for (size_t i = 0; i < num_objects; ++i) {
        std::istringstream datastream(tmp_string);

        for (auto &xj : x) {
          datastream >> xj;
        }

        for (auto &mj : mesh) {
          datastream >> mj;
        }

        datastream >> d;

        datastream >> m;

        for (auto &vj : v) {
          datastream >> vj;
        }

        logger.info("Created a cuboid\n");
        logger.info("Amount of particles: " +
                    std::to_string(std::accumulate(mesh.begin(), mesh.end(), 1,
                                                   std::multiplies<size_t>())));
        logger.info("Dimensions: " + containerToString(mesh));
        logger.info("Distance: " + std::to_string(d));
        logger.info("Mass: " + std::to_string(m));
        logger.info("Initial velocity: " + containerToString(v));

        ParticleGenerator::insertCuboid(x, mesh, d, m, v, 0.1, particles);

        getline(input_file, tmp_string);
      }
    }
  } else {
    logger.warn("Error: could not open file " + std::string(filename));
    exit(-1);
  }
}
