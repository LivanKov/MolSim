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
    logger.info("Parsing XML file: " + filename);

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
    simParameters.enable_brownian = xmlParams.enable_brownian();

    logger.info("Simulation parameters loaded:");
    logger.info("End Time: " + std::to_string(simParameters.end_time));
    logger.info("Delta Time: " + std::to_string(simParameters.time_delta));
    logger.info("Output Base Name: " + simParameters.output_path);
    logger.info("Write Frequency: " +
                std::to_string(simParameters.write_frequency));

    if (xmlParams.gravity().present()) {
      SimParams::enable_gravity = true;
      simParameters.gravity = xmlParams.gravity().get();
      logger.info("g_gravity: " + std::to_string(simParameters.gravity));
    }

    // Read Thermostats
    if (doc->thermostats().present()) {
      SimParams::enable_thermo = true;
      const auto &xmlThermostats = doc->thermostats().get();
      simParameters.initial_temp = xmlThermostats.initial_temp();
      simParameters.target_temp = xmlThermostats.target_temp();
      simParameters.delta_temp = xmlThermostats.delta_temp();
      simParameters.is_gradual = xmlThermostats.is_gradual();
      simParameters.n_thermostats = xmlThermostats.n_thermostats();
      logger.info(" Thermostats loaded:\n");
      logger.info("Initial temperature: " +
                  std::to_string(simParameters.initial_temp));
      logger.info("target temperature: " +
                  std::to_string(simParameters.target_temp));
      logger.info("Temperature difference: " +
                  std::to_string(simParameters.delta_temp));
      logger.info("The number of time steps: " +
                  std::to_string(simParameters.n_thermostats));
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

    // Extract domain size, if no domain is passed, then we use original
    // particle container.
    if (xmlParams.domain_size().present()) {
      if (!doc->boundary_conditions().present()) {
        DomainBoundaryConditions boundary_conditions{
            BoundaryCondition::Outflow, BoundaryCondition::Outflow,
            BoundaryCondition::Outflow, BoundaryCondition::Outflow,
            BoundaryCondition::Outflow, BoundaryCondition::Outflow};
        simParameters.boundaryConditions = boundary_conditions;
      }

      simParameters.linked_cells = true;
      simParameters.domain_size = {xmlParams.domain_size().get().x(),
                                   xmlParams.domain_size().get().y(),
                                   xmlParams.domain_size().get().z()};
      logger.info("Domain Size: " +
                  containerToStrings(simParameters.domain_size));
      simParameters.r_cutoff_radius = xmlParams.r_cutoff_radius();
      logger.info("Cutoff Radius: " +
                  std::to_string(simParameters.r_cutoff_radius));

      if (xmlParams.domain_size().get().lower_left_corner().present()) {
        SimParams::fixed_Domain = true;
        auto lowerLeftCorner =
            xmlParams.domain_size().get().lower_left_corner().get();
        SimParams::lower_left_corner = {
            lowerLeftCorner.x(), lowerLeftCorner.y(), lowerLeftCorner.z()};
        logger.info("Domain is fixed to: " +
                    containerToStrings(simParameters.lower_left_corner));
      }

      if (simParameters.domain_size[2] == 0) {
        particles.initialize(
            {simParameters.domain_size[0], simParameters.domain_size[1]},
            simParameters.r_cutoff_radius, simParameters.boundaryConditions);
        simParameters.dimensions = 2;
      } else {
        particles.initialize(
            {simParameters.domain_size[0], simParameters.domain_size[1],
             simParameters.domain_size[2]},
            simParameters.r_cutoff_radius, simParameters.boundaryConditions);
        simParameters.dimensions = 3;
      }
    }

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

        if (cuboid.additional_force().present()) {
          SimParams::enable_additional_force = true;
          SimParams::additional_force_z_gravity =
              cuboid.additional_force().get().z_grav();
          SimParams::additional_force_time_limit =
              cuboid.additional_force().get().time_limit();

          std::vector<std::array<double, 3>> additional_force_coordinates{};

          for (const auto &coordinate :
               cuboid.additional_force().get().particle_coordinates()) {
            additional_force_coordinates.push_back(std::array<double, 3>{
                coordinate.x(), coordinate.y(), coordinate.z()});
          }
        }

        bool membrane = false;

        if (cuboid.membrane().present()) {
          SimParams::membrane_bond_length = cuboid.membrane().get().r_0();
          SimParams::membrane_stiffness = cuboid.membrane().get().k();
          SimParams::is_membrane = true;
          membrane = true;
        }

        std::array<double, 3> initial_velocity = {
            cuboid.initial_velocity().x(), cuboid.initial_velocity().y(),
            cuboid.initial_velocity().z()};

        logger.info("Creating cuboid: \n");
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

        ParticleGenerator::insertCuboid(position, dimensions, mesh_width, mass,
                                        initial_velocity, particles, epsilon,
                                        sigma, membrane);

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

        logger.info("Creating disc: \n ");
        logger.info("Center: " + containerToStrings(center));
        logger.info("Radius: " + std::to_string(radius));
        logger.info("Mesh Width: " + std::to_string(mesh_width));
        logger.info("Mass: " + std::to_string(mass));
        logger.info("Initial Velocity: " +
                    containerToStrings(initial_velocity));

        ParticleGenerator::insertDisc(center, initial_velocity, radius,
                                      mesh_width, mass, particles, epsilon,
                                      sigma);
      }
    }

    // Extract single particle
    if (doc->particles().present()) {
      for (const auto &p : doc->particles().get().particle()) {
        std::array<double, 3> position = {p.position().x(), p.position().y(),
                                          p.position().z()};

        std::array<double, 3> initial_velocity = {
            p.velocity().x(), p.velocity().y(), p.velocity().z()};

        double mass = p.mass();
        logger.info("Creating single particle: \n");
        logger.info("Position: " + containerToStrings(position));
        logger.info("Initial Velocity: " +
                    containerToStrings(initial_velocity));
        logger.info("Mass: " + std::to_string(mass));

        ParticleGenerator::insertSingleMolecule(position, initial_velocity,
                                                mass, particles);
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
  } else if (value == "Periodic") {
    return BoundaryCondition::Periodic;
  } else {
    throw std::runtime_error("Invalid boundary condition: " + value);
  }
}
