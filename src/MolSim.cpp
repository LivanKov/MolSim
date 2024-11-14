
#include "io/FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "particleSim/ParticleContainer.h"
#include "utils/ArrayUtils.h"

#include "particleSim/ParticleGenerator.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <list>
#include <unistd.h>
#include <unordered_map>
#include <variant>

#include "logger/Logger.h"

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

/*
 * print help flag
 */
void print_help();

// e : time-end
// d : delta
// i : input path
// o : output  path
// t : testing flag -> writes a file for each iteration
// h: help

// TODO: what data structure to pick?
ParticleContainer particles;
double start_time = 0;
double end_time, delta_t;
std::string input_path, output_path;
bool sparse_output = true;
bool xyz_output = false;
bool calculateLJForce = true;
bool enableOutput = true;
std::string log_level;

std::string out_name("MD_vtk");
outputWriter::XYZWriter writer;
outputWriter::VTKWriter v_writer;

int main(int argc, char *argsv[]) {

  if (argc < 2) {
    print_help();
    return 1;
  }

  int opt;

  while ((opt = getopt(argc, argsv, "e:d:i:o:thxl:fn")) != -1) {
    switch (opt) {
    case 'e':
      end_time = atof(optarg);
      break;
    case 'd':
      delta_t = atof(optarg);
      break;
    case 'i':
      input_path = std::string(optarg);
      break;
    case 'o':
      output_path = std::string(optarg);
      break;
    case 't':
      sparse_output = false;
      break;
    case 'h':
      print_help();
      break;
    case 'x':
      xyz_output = true;
      break;
    case 'l':
      log_level = std::string(optarg);
      break;
    case 'f':
      calculateLJForce = false;
      break;
    case 'n':
      enableOutput = false;
      break;
    default:
      fprintf(stderr, "Usage: %s [-h] help\n", argsv[0]);
      return 1;
    }
  }

  Logger &logger = Logger::getInstance(log_level);

  FileReader fileReader;
  fileReader.readFile(particles, input_path.data());

  int iteration = 0;
  double current_time = start_time;

  logger.info("Starting a simulation with:");
  logger.info("\tStart time: " + std::to_string(start_time));
  logger.info("\tEnd time: " + std::to_string(end_time));
  logger.info("\tDelta: " + std::to_string(delta_t));

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    calculateX();
    calculateF();
    calculateV();

    iteration++;
    if (sparse_output && iteration % 10 == 0 && enableOutput)
      plotParticles(iteration);
    else if (!sparse_output && enableOutput)
      plotParticles(iteration);
    logger.trace("Iteration " + std::to_string(iteration) + " finished.");
    current_time += delta_t;
  }

  logger.info("output written. Terminating...");

  logger.debug("Number of particles: " + std::to_string(particles.size()));

  logger.info("Simulation finished.");

  return 0;
}

// ---------------------------------------------------------------------------------------------------------------------

void print_help() {
  std::cout << "Usage: MolSim [options]\n";
  std::cout << "Options:\n";
  std::cout << "  -h                 Show this help message\n";
  std::cout << "  -o   <file_path>   Specify output file path\n";
  std::cout << "  -i   <file_path>   Specify input file path\n";
  std::cout
      << "  -e   <end_time>    Specify how long the simulation should run\n";
  std::cout << "  -d   <time_delta>  Specify time increments\n";
  std::cout << "  -t                 Enable testing mode (Writes a file for "
               "each iteration)\n";
  std::cout << "  -x                 Output .xyz files instead of .vpu\n";
  std::cout << "  -l   <log_level>   Option to choose the logging level "
               "[trace, debug, info, warn, error, off]\n";
  std::cout << "  -f                 Calculate Gravitational Force instead of "
               "Lennard-Jones Force\n";
  std::cout << "  -n   Disable all file output for the sake of performance\n";
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
      if (calculateLJForce) {
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
    x = x + delta_t * v + pow(delta_t, 2) * f / (2 * m);

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
    v = v + delta_t * (old_f + new_f) / (2 * m);

    p.updateV(v);
  }
}

void plotParticles(int iteration) {
  if (xyz_output) {
    writer.plotParticles(particles, output_path + "/" + out_name, iteration);
  } else {
    v_writer.initializeOutput(particles.size());
    for (auto p : particles)
      v_writer.plotParticle(p);
    v_writer.writeFile(output_path + "/" + out_name, iteration);
  }
}
