//
// Created by sebastianpse on 1/27/25.
//

#include "ParticleProfiler.h"

ParticleProfiler::ParticleProfiler(LinkedCellContainer &particles,
                                   size_t number_of_bins, double x_min,
                                   double x_max, const std::string &output_file)
    : particles_(particles), number_of_bins_(number_of_bins), x_min_(x_min),
      x_max_(x_max), output_file_(output_file), header_written_(false) {
  // check for invalid parameters
  if (number_of_bins == 0) {
    Logger::getInstance().error("Number of bins must be greater than 0");
    throw std::invalid_argument("Number of bins must be greater than 0");
  }
  // here again
  if (x_min_ > x_max_) {
    Logger::getInstance().error("Invalid arguments for x_min and x_max!");
    throw std::invalid_argument("Not allowed: x_min > x_max");
  }
  // calculates the bin width
  bin_width_ = (x_max_ - x_min_) / static_cast<double>(number_of_bins_);
}

void ParticleProfiler::apply_profiler() {

  std::ofstream file(output_file_, std::ios::app);
  // when file can not be opened
  if (!file.is_open()) {
    Logger::getInstance().error("Error: Unable to open file for writing" +
                                output_file_);
    return;
  }

  if (!header_written_) {
    // write the header for the CSV (only once)
    file << "Bin index,middle X-Position of Bin,Density (Particle per "
            "Binwidth),Average Velocity-x,Average Velocity-y,Average "
            "Velocity-z\n";
    header_written_ = true;
  }

  std::vector<size_t> particle_count(number_of_bins_, 0);
  std::vector<std::array<double, 3>> total_velocity(number_of_bins_);

  // iterate over each particle here
  for (size_t i = 0; i < particles_.size(); ++i) {
    auto &particle = particles_[i];

    auto x_position = particle.getX()[0];

    // check if particle is in correct scope
    if (x_position >= x_min_ && x_position <= x_max_) {
      // calculate the correct bin with respect to position and x_min_
      auto bin_index = static_cast<size_t>((x_position - x_min_) / bin_width_);

      // when you have an x at x_max, bin_index would be 1 too large -> so
      // reduce by 1
      if (bin_index == number_of_bins_) {
        bin_index--; // Ensure it goes into the last bin
      }

      // additional safety check
      if (bin_index < number_of_bins_) {
        // increase particle count, since particle is in this bin
        particle_count[bin_index]++;

        // also add the particle's velocity to the bin
        total_velocity[bin_index][0] += particle.getV()[0]; // x-value
        total_velocity[bin_index][1] += particle.getV()[1]; // y-value
        total_velocity[bin_index][2] += particle.getV()[2]; // z-value
      }
    }
  }

  Logger::getInstance().info("Finished with iterating!");

  // now write the correct values for each bin to the csv file
  for (size_t i = 0; i < number_of_bins_; ++i) {
    size_t n_particles_bin = particle_count[i];
    if (n_particles_bin > 0) {
      // compute the density
      double density = static_cast<double>(n_particles_bin) / bin_width_;

      // Compute the average velocity in each direction
      double avg_velo_x =
          total_velocity[i][0] / static_cast<double>(n_particles_bin);
      double avg_velo_y =
          total_velocity[i][1] / static_cast<double>(n_particles_bin);
      double avg_velo_z =
          total_velocity[i][2] / static_cast<double>(n_particles_bin);

      // Get the x-position of the center of the bin
      double bin_center_x = x_min_ + (i + 0.5) * bin_width_;

      // Write the data to the CSV file
      file << i << "," << bin_center_x << "," << density << "," << avg_velo_x
           << "," << avg_velo_y << "," << avg_velo_z << "\n";
    }
  }
}
