#include "XMLReader.h"
#include "MolSim.hxx"
#include "cli/SimParams.h"
#include "simulator/particle/ParticleGenerator.h"
#include "simulator/particle/container/LinkedCellContainer.h"
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

void XMLReader::readXMLFile(LinkedCellContainer &particles,
                            SimParams &simParameters) {
  Logger &logger = Logger::getInstance(simParameters.log_level);
  std::string filename = simParameters.input_path;
  try {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
      throw std::runtime_error("Could not open file: " + filename);
    }

    std::unique_ptr<MolSim> doc =
        MolSim_(filename, xml_schema::flags::dont_validate);

    logger.info("Parsed XML file successfully");

    // Extract simulation parameters
    auto xmlParams = doc->simulation_parameters();
    simParameters.end_time = (xmlParams.end_time() != simParameters.end_time &&
                              simParameters.end_time == 0.0)
                                 ? xmlParams.end_time()
                                 : simParameters.end_time;
    simParameters.time_delta =
        (xmlParams.delta_time() != simParameters.time_delta &&
         simParameters.time_delta == 0.0)
            ? xmlParams.delta_time()
            : simParameters.time_delta;
    simParameters.output_path =
        (xmlParams.output_basename() != simParameters.output_path &&
         simParameters.output_path == "")
            ? xmlParams.output_basename()
            : simParameters.output_path;
    simParameters.write_frequency =
        (xmlParams.write_frequency() != simParameters.write_frequency &&
         simParameters.write_frequency == 0)
            ? xmlParams.write_frequency()
            : simParameters.write_frequency;
    simParameters.r_cutoff_radius = xmlParams.r_cutoff_radius();
    simParameters.enable_brownian = xmlParams.enable_brownian();

    logger.info("Simulation parameters loaded:");
    logger.info("End Time: " + std::to_string(simParameters.end_time));
    logger.info("Delta Time: " + std::to_string(simParameters.time_delta));
    logger.info("Output Base Name: " + simParameters.output_path);
    logger.info("Write Frequency: " +
                std::to_string(simParameters.write_frequency));
    logger.info("Cutoff Radius: " +
                std::to_string(simParameters.r_cutoff_radius));

    // Extract domain size, if no domain is passed, then we use original
    // particle container.
    if (xmlParams.domain_size().present()) {
      simParameters.linked_cells = true;
      simParameters.domain_size = {xmlParams.domain_size().get().x(),
                                   xmlParams.domain_size().get().y(),
                                   xmlParams.domain_size().get().z()};
      logger.info("Domain Size: " +
                  containerToStrings(simParameters.domain_size));
    }

    if (xmlParams.gravity().present()) {
      simParameters.gravity = xmlParams.gravity().get();
      logger.info("g_gravity: " + std::to_string(simParameters.gravity));
    }

    // Read Thermostats
    if (doc->thermostats().present()) {
      const auto &xmlThermostats = doc->thermostats().get();
      simParameters.initial_temp = xmlThermostats.initial_temp();
      simParameters.target_temp = xmlThermostats.target_temp();
      simParameters.delta_temp = xmlThermostats.delta_temp();
      simParameters.is_gradual = xmlThermostats.is_gradual();
      simParameters.n_thermostats = xmlThermostats.n_thermostats();
      logger.info("\n Thermostats loaded:");
      logger.info("Initial temperature: " + std::to_string(simParameters.initial_temp));
      logger.info("target temperature: " + std::to_string(simParameters.target_temp));
      logger.info("Temperature difference: " + std::to_string(simParameters.delta_temp));
      logger.info("The number of time steps: " + std::to_string(simParameters.n_thermostats));
    }

    // Read boundary conditions
    if (doc->boundary_conditions().present()) {
      const auto &xmlBoundaryConditions = doc->boundary_conditions().get();
      simParameters.boundaryConditions.left =
          parseBoundaryCondition(xmlBoundaryConditions.left());
      simParameters.boundaryConditions.right =
          parseBoundaryCondition(xmlBoundaryConditions.right());
      simParameters.boundaryConditions.top =
          parseBoundaryCondition(xmlBoundaryConditions.top());
      simParameters.boundaryConditions.bottom =
          parseBoundaryCondition(xmlBoundaryConditions.bottom());

      if (xmlBoundaryConditions.front().present()) {
        simParameters.boundaryConditions.front =
            parseBoundaryCondition(xmlBoundaryConditions.front().get());
      }

      if (xmlBoundaryConditions.back().present()) {
        simParameters.boundaryConditions.back =
            parseBoundaryCondition(xmlBoundaryConditions.back().get());
      }

      logger.info("Boundary conditions loaded successfully");
    }

    std::initializer_list<double> init_list;

    if (simParameters.domain_size[2] == 0) {
      init_list = {simParameters.domain_size[0], simParameters.domain_size[1]};
      simParameters.dimensions = 2;
    } else {
      init_list = {simParameters.domain_size[0], simParameters.domain_size[1],
                   simParameters.domain_size[2]};
      simParameters.dimensions = 3;
    }

    // Initialize the linkedcellcontainer with received parameters from XMLinput
    particles.initialize(init_list, simParameters.r_cutoff_radius,
                         simParameters.boundaryConditions);

    // Extract cuboid specification
    if (doc->cuboids().present()) {
      logger.info("Number of cuboids found: " +
                  std::to_string(doc->cuboids().get().cuboid().size()));

      for (const auto &cuboid : doc->cuboids().get().cuboid()) {
        std::array<double, 3> position = {cuboid.coordinate().x(),
                                          cuboid.coordinate().y(),
                                          cuboid.coordinate().z()};

        std::array<size_t, 3> dimensions = {cuboid.dimensions().x(),
                                            cuboid.dimensions().y(),
                                            cuboid.dimensions().z()};
        double mesh_width = cuboid.mesh_width();
        double mass = cuboid.mass();

        // Specify epsilon and sigma of the cuboid
        double epsilon = cuboid.epsilon();
        double sigma = cuboid.sigma();

        std::array<double, 3> initial_velocity = {
            cuboid.initial_velocity().x(), cuboid.initial_velocity().y(),
            cuboid.initial_velocity().z()};
        double avg_velocity = cuboid.average_velocity();

        logger.info("\n Creating cuboid: ");
        logger.info(
            "Amount of particles: " +
            std::to_string(std::accumulate(dimensions.begin(), dimensions.end(),
                                           1, std::multiplies<size_t>())));
        logger.info("Position: " + containerToStrings(position));
        logger.info("Dimensions: " + containerToStrings(dimensions));
        logger.info("Mesh Width: " + std::to_string(mesh_width));
        logger.info("Mass: " + std::to_string(mass));
        logger.info("Initial Velocity: " +
                    containerToStrings(initial_velocity));
        logger.info("Average Velocity: " + std::to_string(avg_velocity));

        ParticleGenerator::insertCuboid(position, dimensions, mesh_width, mass,
                                        initial_velocity, avg_velocity,
                                        particles);
        logger.info("Particles check: " + std::to_string(particles.size()));
        logger.info("Particles' cell check: " +
                    std::to_string(particles.cells.size()));
      }
    }

    // Extract disc specification
    if (doc->discs().present()) {
      for (const auto &disc : doc->discs().get().disc()) {
        std::array<double, 3> center = {disc.center().x(), disc.center().y(),
                                        disc.center().z()};

        std::array<double, 3> initial_velocity = {disc.initial_velocity().x(),
                                                  disc.initial_velocity().y(),
                                                  disc.initial_velocity().z()};

        size_t radius = disc.radius();
        double mesh_width = disc.mesh_width();
        double mass = disc.mass();

        // Specify epsilon and sigma of the disc
        double epsilon = disc.epsilon();
        double sigma = disc.sigma();
        double avg_velocity = disc.average_velocity();

        logger.info("\n Creating disc: ");
        logger.info("Center: " + containerToStrings(center));
        logger.info("Radius: " + std::to_string(radius));
        logger.info("Mesh Width: " + std::to_string(mesh_width));
        logger.info("Mass: " + std::to_string(mass));
        logger.info("Initial Velocity: " +
                    containerToStrings(initial_velocity));
        logger.info("Average Velocity: " + std::to_string(avg_velocity));

        ParticleGenerator::insertDisc(center, initial_velocity, radius,
                                      mesh_width, mass, particles);
      }
    }

  } catch (const std::exception &e) {
    logger.warn("Error while reading XML file: " + std::string(filename));
    exit(-1);
  }
}

BoundaryCondition XMLReader::parseBoundaryCondition(const std::string &value) {
  if (value == "Outflow") {
    return BoundaryCondition::Outflow;
  } else if (value == "Reflecting") {
    return BoundaryCondition::Reflecting;
  } else {
    throw std::runtime_error("Invalid boundary condition: " + value);
  }
}
