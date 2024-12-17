#include "CheckpointReader.h"
#include "utils/logger/Logger.h"
#include <fstream>
#include <iomanip>

CheckpointReader::CheckpointReader() = default;
CheckpointReader::~CheckpointReader() = default;

void CheckpointReader::readCheckpoint(LinkedCellContainer &particles,
                                      double &delta_t, double &t_end) {
  const std::string &filename = "../output/checkpoint.chk";
  std::ifstream file(filename);
  if (!file.is_open()) {
    Logger::getInstance().error("Failed to open checkpoint file for reading: " +
                                filename);
    throw std::runtime_error("Cannot open checkpoint file");
  }

  Logger::getInstance().info("Reading checkpoint file: " + filename);

  std::string line;
  size_t particle_count = 0;

  // Read header
  while (std::getline(file, line)) {
    if (line[0] == '#') {
      // Parse metadata
      if (line.find("delta_t:") != std::string::npos) {
        std::istringstream(line.substr(line.find(":") + 1)) >> delta_t;
      } else if (line.find("t_end:") != std::string::npos) {
        std::istringstream(line.substr(line.find(":") + 1)) >> t_end;
      }
    } else {
      // First non-comment line is the particle count
      particle_count = std::stoul(line);
      break;
    }
  }

  Logger::getInstance().info("Number of particles: " +
                             std::to_string(particle_count));

  // Read particle data
  for (size_t i = 0; i < particle_count; ++i) {
    double x, y, z, vx, vy, vz, mass, oldFx, oldFy, oldFz;
    int type;

    file >> x >> y >> z >> vx >> vy >> vz >> mass >> oldFx >> oldFy >> oldFz >>
        type;

    if (file.fail()) {
      Logger::getInstance().error("Error reading particle data at line " +
                                  std::to_string(i + 1));
      throw std::runtime_error("Error reading checkpoint file");
    }

    // Create particle and insert into the container
    Particle particle({x, y, z}, {vx, vy, vz}, mass, particles.particle_id);
    particles.particle_id++;
    particles.insert(particle, true);
    particles.readjust();
  }

  file.close();
  Logger::getInstance().info("Checkpoint successfully read from " + filename);
}