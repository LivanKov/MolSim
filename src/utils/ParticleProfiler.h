//
// Created by sebastianpse on 1/27/25.
//

#pragma once
#include "../simulator/particle/container/LinkedCellContainer.h"
#include <fstream>
#include <iostream>

/**
 *
 * This class offers a component which divides the system into bins along the
 * x-axis. It calculates the density and average velociy of the particles for
 * each bin.
 */
class ParticleProfiler {
public:
  /**
   * @brief Constructor
   *
   * @param particles: The particles of the system.
   * @param number_of_bins: The number of bins in which the system should be
   * divided.
   * @param x_min: The x-point where the system of the bins should begin.
   * @param x_max: The x-point where the system of the bins should end.
   * @param output_file: The output file that gets written with the values
   *
   */
  ParticleProfiler(LinkedCellContainer &particles, size_t number_of_bins,
                   double x_min, double x_max, const std::string &output_file);

  /**
   * @brief Writes the current state of the system to the output file.
   *
   */
  void apply_profiler();

private:
  LinkedCellContainer &particles_;
  size_t number_of_bins_;
  double x_min_;
  double x_max_;
  double bin_width_;
  std::string output_file_;
  /**
   * @brief Checks if the header is written to ensure that the header is only
   * written once.
   */
  bool header_written_;
};
