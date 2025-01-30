//
// Created by sebastianpse on 1/28/25.
//

#include <gtest/gtest.h>

#include "simulator/particle/ParticleGenerator.h"
#include "simulator/particle/container/LinkedCellContainer.h"
#include "utils/ParticleProfiler.h"

#include <filesystem>


class ParticleProfilerTest : public ::testing::Test {
protected:
    LinkedCellContainer particles{};

};

TEST(ParticleProfilerTest, BinAssignment) {
    LinkedCellContainer particles{};

    // Insert particles with known positions
    Particle particle1{{5.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{15.0, 0.0, 0.0}, {0.5, 0.0, 0.0}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{25.0, 0.0, 0.0}, {0.0, 0.5, 0.0}, 1.0, 0};
    particles.insert(particle3);

    // init the profiler with correct bin configuration
    size_t bins = 10;
    double x_min = 0.0;
    double x_max = 50.0;
    std::string output_file = "test_output.csv";

    ParticleProfiler profiler(particles, bins, x_min, x_max, output_file);

    profiler.apply_profiler();

    // open output file to read results
    std::ifstream file(output_file);
    ASSERT_TRUE(file.is_open());

    // skip header line
    std::string line;
    std::getline(file, line);

    // init vector to track bin counts
    std::vector<size_t> bin_counts(bins, 0);

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        size_t bin_index;
        ss >> bin_index;
        bin_counts[bin_index]++;
    }
    file.close();

    // check that the particles are assigned to the correct bins
    ASSERT_EQ(bin_counts[1], 1);  // Particle at x=5.0 should be in bin 1
    ASSERT_EQ(bin_counts[3], 1);  // Particle at x=15.0 should be in bin 3
    ASSERT_EQ(bin_counts[5], 1);  // Particle at x=25.0 should be in bin 5

    // remove output file after test
    std::filesystem::remove(output_file);
}


