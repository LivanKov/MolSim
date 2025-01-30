#include "CheckpointWriter.h"
#include "utils/logger/Logger.h"
#include <fstream>
#include <iomanip>

CheckpointWriter::CheckpointWriter() = default;
CheckpointWriter::~CheckpointWriter() = default;

void CheckpointWriter::writeCheckpoint(LinkedCellContainer &particles,
                                       const std::string &filename,
                                       double delta_t, double t_end) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    Logger::getInstance().error("Failed to open checkpoint file for writing: " +
                                filename);
    throw std::runtime_error("Cannot open checkpoint file");
  }

  // Write header information
  file << "# Checkpoint file for simulation restart\n";
  file << "# delta_t: " << delta_t << "\n";
  file << "# t_end: " << t_end << "\n";
  file << particles.size() << "\n";

  // Write particle data
  for (const auto &p : particles) {
    file << std::fixed << std::setprecision(6) << p.getX()[0] << " "
         << p.getX()[1] << " " << p.getX()[2] << " " << p.getV()[0] << " "
         << p.getV()[1] << " " << p.getV()[2] << " " << p.getM() << " "
         << p.getId() << "\n";
  }

  file.close();
  Logger::getInstance().info("Checkpoint written to " + filename);
}