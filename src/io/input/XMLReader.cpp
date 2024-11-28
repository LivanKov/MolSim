#include "XMLReader.h"
#include "MolSim.hxx"
#include "cli/SimParams.h"
#include "simulator/particle/ParticleContainer.h"
#include "simulator/particle/ParticleGenerator.h"
#include "utils/logger/Logger.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

XMLReader::XMLReader() = default;

XMLReader::~XMLReader() = default;

auto containerToStrings = [](const auto &container) {
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

void XMLReader::readXMLFile(ParticleContainer &particles,
                            SimParams &simParameters,
                            const std::string &filename) {
  Logger &logger = Logger::getInstance();
  try {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
      throw std::runtime_error("Could not open file: " + filename);
    }

    std::unique_ptr<MolSim> doc =
        MolSim_(filename, xml_schema::flags::dont_validate);

    logger.info("Parsed XML file successfully");

    // Extract simulation parameters
    auto simParams = doc->simulation_parameters();
    simParameters.end_time = simParams.end_time();
    simParameters.time_delta = simParams.delta_time();
    simParameters.output_path = simParams.output_basename();
    simParameters.write_frequency = simParams.write_frequency();
    simParameters.r_cutoff_radius = simParams.r_cutoff_radius();
    simParameters.domain_size = {simParams.domain_size().x(),
                                 simParams.domain_size().y(),
                                 simParams.domain_size().z()};

    logger.info("Simulation parameters loaded:");
    logger.info("End Time: " + std::to_string(simParameters.end_time));
    logger.info("Delta Time: " + std::to_string(simParameters.time_delta));
    logger.info("Output Base Name: " + simParameters.output_path);
    logger.info("Write Frequency: " +
                std::to_string(simParameters.write_frequency));
    logger.info("Cutoff Radius: " + std::to_string(simParameters.r_cutoff_radius));
    logger.info("Domain Size: " + containerToStrings(simParameters.domain_size));

    // auto cuboids = doc->cuboids();

    // Extract cuboid specification
    logger.info("Number of cuboids found: " +
                std::to_string(doc->cuboids().cuboid().size()));

    for (const auto &cuboid : doc->cuboids().cuboid()) {
      // const auto &cuboid = cuboids_instance.cuboid();
      std::array<double, 3> position = {cuboid.coordinate().x(),
                                        cuboid.coordinate().y(),
                                        cuboid.coordinate().z()};

      std::array<size_t, 3> dimensions = {cuboid.dimensions().x(),
                                          cuboid.dimensions().y(),
                                          cuboid.dimensions().z()};
      double mesh_width = cuboid.mesh_width();
      double mass = cuboid.mass();
      std::array<double, 3> initial_velocity = {cuboid.initial_velocity().x(),
                                                cuboid.initial_velocity().y(),
                                                cuboid.initial_velocity().z()};
      double avg_velocity = cuboid.average_velocity();

      logger.info("Creating cuboid: \n");
      logger.info(
          "Amount of particles: " +
          std::to_string(std::accumulate(dimensions.begin(), dimensions.end(),
                                         1, std::multiplies<size_t>())));
      logger.info("Position: " + containerToStrings(position));
      logger.info("Dimensions: " + containerToStrings(dimensions));
      logger.info("Mesh Width: " + std::to_string(mesh_width));
      logger.info("Mass: " + std::to_string(mass));
      logger.info("Initial Velocity: " + containerToStrings(initial_velocity));
      logger.info("Average Velocity: " + std::to_string(avg_velocity));

      ParticleGenerator::insertCuboid(position, dimensions, mesh_width, mass,
                                      initial_velocity, avg_velocity,
                                      particles);
    }

  } catch (const std::exception &e) {
    logger.warn("Error while reading XML file: " + std::string(filename));
    exit(-1);
  }
}