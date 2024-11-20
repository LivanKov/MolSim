#include "io/input/FileReader.h"
#include "io/output/VTKWriter.h"
#include "io/output/XYZWriter.h"
#include "particleSim/ParticleContainer.h"
#include "utils/ArrayUtils.h"

#include "particleSim/ParticleGenerator.h"

#include <spdlog/spdlog.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <list>
#include <unordered_map>
#include <variant>

#include "logger/Logger.h"
#include "utils/SimParams.h"
#include "utils/CommandParser.h"

/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 * calculate the position for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);


ParticleContainer particles{};
SimParams parameters{};
std::string out_name("MD_vtk");
outputWriter::XYZWriter writer;
outputWriter::VTKWriter v_writer;


int main(int argc, char *argsv[]) {

  CommandParser::parse(argc,argsv,parameters);

  Logger &logger = Logger::getInstance(parameters.log_level);

  FileReader fileReader;
  fileReader.readFile(particles, parameters.input_path.data());

  int iteration = 0;
  double current_time = parameters.start_time;

  logger.warn("Starting a simulation with:");
  logger.info("\tStart time: " + std::to_string(parameters.start_time));
  logger.info("\tEnd time: " + std::to_string(parameters.end_time));
  logger.info("\tDelta: " + std::to_string(parameters.time_delta));
/*
@brief
*/
  // for this loop, we assume: current x, current f and current v are known
  while (current_time < parameters.end_time) {
    calculateX();
    calculateF();
    calculateV();

    iteration++;
    if (parameters.sparse_output && iteration % 10 == 0 && parameters.enable_output)
      plotParticles(iteration);
    else if (!parameters.sparse_output && parameters.enable_output)
      plotParticles(iteration);
    logger.trace("Iteration " + std::to_string(iteration) + " finished.");
    current_time += parameters.time_delta;
  }

  logger.info("output written. Terminating...");

  logger.debug("Number of particles: " + std::to_string(particles.size()));

  logger.warn("Simulation finished.");

  return 0;
}


void calculateF() {
  // store the current force as the old force and reset current to 0
  for (auto &p : particles) {
    p.updateOldF(p.getF());
    p.updateF(0, 0, 0);
  }

  // Iterate each pair
  for (auto it = particles.pair_begin(); it != particles.pair_end(); ++it) {
    ParticlePair &pair = *it;
    Particle &p1 = *(pair.first);
    Particle &p2 = *(pair.second);
    auto r12 = p2.getX() - p1.getX();
    // distance ||x_i - x_j ||
    double distance = ArrayUtils::L2Norm(r12);

    // avoid extermely small distance
    if (distance > 1e-5) {
      // switch Lennard-Jones/ Simple force
      double totalForce;
      if (parameters.calculate_lj_force) {
        // Lennard-Jones parameters
        const double epsilon = 5.0;
        const double sigma = 1.0;
        // Lennard-Jones Force Formula (3)
        double term = sigma / distance;
        double term6 = pow(term, 6);
        double term12 = pow(term, 12);
        totalForce = 24 * epsilon * (term6 - 2 * term12) / distance;
      } else {
        // Simple Force Calculation Formula (14)
        totalForce = p1.getM() * p2.getM() / pow(distance, 2);
      }
      auto force = (totalForce / distance) * r12;

      p1.updateF(p1.getF() + force);
      // Newton's third law
      p2.updateF(p2.getF() - force);
    }
  }
}

void calculateX() {
  for (auto &p : particles) {
    auto x = p.getX();
    auto v = p.getV();
    auto f = p.getF();
    double m = p.getM();

    // Velocity-Störmer-Verlet formula (8)
    x = x + parameters.time_delta * v + pow(parameters.time_delta, 2) * f / (2 * m);

    p.updateX(x);
  }
}

void calculateV() {
  for (auto &p : particles) {
    auto v = p.getV();
    auto old_f = p.getOldF();
    auto new_f = p.getF();
    double m = p.getM();

    // Velocity-Störmer-Verlet formula (9)
    v = v + parameters.time_delta * (old_f + new_f) / (2 * m);

    p.updateV(v);
  }
}

void plotParticles(int iteration) {
  if (parameters.xyz_output) {
    writer.plotParticles(particles, parameters.output_path + "/" + out_name, iteration);
  } else {
    v_writer.initializeOutput(particles.size());
    for (auto p : particles)
      v_writer.plotParticle(p);
    v_writer.writeFile(parameters.output_path + "/" + out_name, iteration);
  }
}
