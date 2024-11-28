#include "XMLReader.h"
#include "MolSim.hxx"
#include "simulator/particle/ParticleContainer.h"
#include "simulator/particle/ParticleGenerator.h"
#include "utils/logger/Logger.h"
#include <fstream>
#include <iostream>
#include <sstream>

XMLReader::XMLReader() = default;

XMLReader::~XMLReader() = default;

void XMLReader::readXMLFile(ParticleContainer &particles, SimulationParameters &simParameters,
                            const std::string &filename) {
  try {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
      throw std::runtime_error("Could not open file: " + filename);
    }

    std::unique_ptr<MolSim> doc =
        MolSim_(filename, xml_schema::flags::dont_validate);

    Logger::getInstance().info("Parsed XML file successfully");

    // Extract simulation parameters
    auto simParams = doc->simulation_parameters();
    simParameters.t_end = simParams.end_time();
    simParameters.delta_t = simParams.delta_time();
    simParameters.output_basename = simParams.output_basename();
    simParameters.write_frequency = simParams.write_frequency();

    Logger::getInstance().info("Simulation parameters loaded:");
    Logger::getInstance().info("End Time: " +
                               std::to_string(simParameters.t_end));
    Logger::getInstance().info("Delta Time: " +
                                std::to_string(simParameters.delta_t));
    Logger::getInstance().info("Output Base Name: " +
                               simParameters.output_basename);
    Logger::getInstance().info("Write Frequency: " +
                               std::to_string(simParameters.write_frequency));

    // auto cuboids = doc->cuboids();


    // Extract cuboid specification
    Logger::getInstance().info("Number of cuboids found: " +
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

      Logger::getInstance().info("Creating cuboid:");
      Logger::getInstance().info("Position: " + containerToString(position));
      Logger::getInstance().info("Dimensions: " +
                                      containerToString(dimensions));
      Logger::getInstance().info("Mesh Width: " +
                                      std::to_string(mesh_width));
      Logger::getInstance().info("Mass: " + std::to_string(mass));
      Logger::getInstance().info("Initial Velocity: " +
                                      containerToString(initial_velocity));
      Logger::getInstance().info("Average Velocity: " +
                                      std::to_string(avg_velocity));

      ParticleGenerator::insertCuboid(position, dimensions, mesh_width, mass,
                                      initial_velocity, avg_velocity,
                                      particles);
    } 

  } catch (const std::exception &e) {
    std::cerr << "Error while reading XML file: " << e.what() << std::endl;
    exit(-1);
  }
}