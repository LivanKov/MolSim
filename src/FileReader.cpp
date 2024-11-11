/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"
#include "ParticleContainer.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

#include "Logger.h"

FileReader::FileReader() = default;

FileReader::~FileReader() = default;

void FileReader::readFile(ParticleContainer &particles, char *filename) {
  std::array<double, 3> x;
  std::array<double, 3> v;
  double m;
  int num_particles = 0;

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

      numstream >> num_particles;

      logger.info("    Reading " + num_particles);
      getline(input_file, tmp_string);
      logger.debug("Read line: " + tmp_string);

      for (int i = 0; i < num_particles; i++) {
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
      logger.debug("Read line: " + tmp_string);
      std::array<double, 3> xyz;
      std::istringstream xyz_datastream(tmp_string);
      for (auto &p : xyz) {
        xyz_datastream >> p;
      }

      getline(input_file, tmp_string);
      logger.debug("Read line: " + tmp_string);
      std::array<size_t, 3> cube_dim;
      std::istringstream dim_datastream(tmp_string);
      for (auto &p : cube_dim) {
        dim_datastream >> p;
      }

      getline(input_file, tmp_string);
      logger.debug("Read line: " + tmp_string);
      std::istringstream dist_datastream(tmp_string);
      double distance;
      dist_datastream >> distance;

      getline(input_file, tmp_string);
      logger.debug("Read line: " + tmp_string);
      std::istringstream mass_datastream(tmp_string);
      double mass;
      mass_datastream >> mass;

      getline(input_file, tmp_string);
      logger.debug("Read line: " + tmp_string);
      std::istringstream velocity_datastream(tmp_string);
      std::array<double, 3> velocity;
      for (auto &p : velocity) {
        velocity_datastream >> p;
      }

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
      std::cout << containerToString(xyz) << std::endl;
      std::cout << containerToString(cube_dim) << std::endl;
 
      logger.info("Created a cuboid\n");
      logger.info("Amount of particles: " +
                  std::to_string(std::accumulate(cube_dim.begin(), cube_dim.end(), 1,
                                  std::multiplies<size_t>())));
      logger.info("Dimensions: " + containerToString(cube_dim));
      logger.info("Distance: " + std::to_string(distance));
      logger.info("Mass: " + std::to_string(mass));
      logger.info("Initial velocity: " + containerToString(velocity));
    }
  } else {
    std::cout << "Error: could not open file " << filename << std::endl;
    // logger.warn("Error: could not open file "  filename);
    exit(-1);
  }
}
