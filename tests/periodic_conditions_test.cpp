#include "simulator/particle/ParticleGenerator.h"
#include "simulator/particle/container/LinkedCellContainer.h"
#include "io/input/cli/SimParams.h"
#include "utils/logger/Logger.h"
#include "simulator/calculations/BoundaryConditions.h"
#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include "simulator/calculations/Calculation.h"
#include "simulator/calculations/Position.h"

class PeriodicBoundaryTest : public ::testing::Test {
protected:
    LinkedCellContainer container;  

    LinkedCellContainer b_container;

    LinkedCellContainer d_container;

    PeriodicBoundaryTest() : container{}, b_container{}, d_container{} {
        container.initialize({10.0, 10.0}, 2.5,
            {BoundaryCondition::Periodic, // Left
             BoundaryCondition::Periodic, // Right
             BoundaryCondition::Periodic, // Top
             BoundaryCondition::Periodic  // Bottom
            });

        b_container.initialize({10.0, 10.0}, 2.5,
            {BoundaryCondition::Periodic,   // Left
             BoundaryCondition::Reflecting, // Right
             BoundaryCondition::Periodic,   // Top
             BoundaryCondition::Periodic    // Bottom
            });
        d_container.initialize({4.0, 4.0, 6.0}, 2.0,
            {BoundaryCondition::Periodic, // Left
             BoundaryCondition::Periodic, // Right
             BoundaryCondition::Periodic, // Top
             BoundaryCondition::Periodic, // Bottom
             BoundaryCondition::Periodic, // Front
             BoundaryCondition::Periodic  // Back
            });
    }

};


TEST_F(PeriodicBoundaryTest, BasicNeighbourTest) {
    SimParams::fixed_Domain = false;
    ParticleGenerator::insertCuboid(
      std::array<double, 3>{0.0, 0.0, 0.0}, std::array<size_t, 3>{4, 4, 1}, 2.5,
      1.0, std::array<double, 3>{0.0, 0.0, 0.0}, container);

    Calculation<BoundaryConditions>::run(container);


    // Test basic container setup
    ASSERT_EQ(container.size(), 16) << "Container should have 16 particles";
    ASSERT_EQ(container.halo_cell_indices.size(), 12) << "Should have 12 halo cells";

    // Verify cell contents
    for(size_t i = 0; i < container.cells.size(); i++){
        auto& cell = container.cells[i];
        ASSERT_EQ(cell.size(), 1) << "Cell " << i << " should contain exactly 1 particle";
        ASSERT_TRUE(std::find(cell.particle_ids.begin(), cell.particle_ids.end(), i) != cell.particle_ids.end()) << "Cell " << i << " should contain particle " << i;
    }

    // Test corner cell properties
    auto& corner_cell = container.cells[0];
    ASSERT_TRUE(corner_cell.is_halo) << "Corner cell should be a halo cell";
    ASSERT_EQ(corner_cell.placement, Placement::BOTTOM_LEFT_CORNER) << "Corner cell should be bottom-left";

    // Get and verify neighbors for corner particle
    auto neighbours = container.get_neighbours(0);

    auto additional_neighbour_indices = container.get_periodic_neighbours(0);
    
    // Print debug info
    std::cout << "Number of neighbors found: " << neighbours.size() << std::endl;
    std::cout << "Neighbor particle IDs: ";
    for(const auto& n : neighbours) {
        std::cout << n->getId() << " ";
    }
    std::cout << std::endl;

    // TODO: Uncomment and adjust expected neighbor count once verified
    EXPECT_EQ(neighbours.size() + additional_neighbour_indices.size(), 8) << "Corner particle should have 8 neighbors with periodic conditions";

    std::vector<int> expected_ids = {1, 4, 5, 12, 13, 3, 7, 15};
    std::vector<int> actual_ids;
    
    for (const auto& n : neighbours) {
        actual_ids.push_back(n->getId());
    }

    for (const auto& n : additional_neighbour_indices) {
        actual_ids.push_back(n.id);
    }

    std::sort(actual_ids.begin(), actual_ids.end());
    std::sort(expected_ids.begin(), expected_ids.end());
    
    EXPECT_EQ(actual_ids, expected_ids) << "Neighbor IDs do not match expected values";


    // Test bottom right corner
    auto& bottom_right_cell = container.cells[3];
    ASSERT_TRUE(bottom_right_cell.is_halo) << "Bottom right cell should be a halo cell";
    ASSERT_EQ(bottom_right_cell.placement, Placement::BOTTOM_RIGHT_CORNER) << "Cell should be bottom-right";

    // Get and verify neighbors for bottom right particle
    auto bottom_right_neighbours = container.get_neighbours(3);

    additional_neighbour_indices = container.get_periodic_neighbours(3);
    
    std::cout << "Number of neighbors for bottom right particle: " << bottom_right_neighbours.size() << std::endl;
    std::cout << "Bottom right neighbor particle IDs: ";
    for(const auto& n : bottom_right_neighbours) {
        std::cout << n->getId() << " ";
    }
    std::cout << std::endl;

    EXPECT_EQ(bottom_right_neighbours.size() + additional_neighbour_indices.size(), 8) << "Bottom right particle should have 8 neighbors with periodic conditions";

    std::vector<int> bottom_right_expected_ids = {2, 3, 6, 7, 0, 4, 14, 15, 12};
    std::vector<int> bottom_right_actual_ids;
    for (const auto& n : bottom_right_neighbours) {
        bottom_right_actual_ids.push_back(n->getId());
    }

    for (const auto& n : additional_neighbour_indices) {
        bottom_right_actual_ids.push_back(n.id);
    }

    std::sort(bottom_right_actual_ids.begin(), bottom_right_actual_ids.end());
    std::sort(bottom_right_expected_ids.begin(), bottom_right_expected_ids.end());
    
    EXPECT_EQ(bottom_right_actual_ids, bottom_right_expected_ids) << "Bottom right neighbor IDs do not match expected values";


    // Test upper left corner
    auto& upper_left_cell = container.cells[12];
    ASSERT_TRUE(upper_left_cell.is_halo) << "Upper left cell should be a halo cell";
    ASSERT_EQ(upper_left_cell.placement, Placement::TOP_LEFT_CORNER) << "Cell should be top-left";

    // Get and verify neighbors for upper left particle
    auto upper_left_neighbours = container.get_neighbours(12);

    additional_neighbour_indices = container.get_periodic_neighbours(12);
    
    std::cout << "Number of neighbors for upper left particle: " << upper_left_neighbours.size() << std::endl;
    std::cout << "Upper left neighbor particle IDs: ";
    for(const auto& n : upper_left_neighbours) {
        std::cout << n->getId() << " ";
    }
    std::cout << std::endl;

    EXPECT_EQ(upper_left_neighbours.size() + additional_neighbour_indices.size(), 8) << "Upper left particle should have 8 neighbors with periodic conditions";

    std::vector<int> upper_left_expected_ids = {8, 9, 12, 13, 0, 1, 11, 15, 3};
    std::vector<int> upper_left_actual_ids;
    for (const auto& n : upper_left_neighbours) {
        upper_left_actual_ids.push_back(n->getId());
    }

    for (const auto& n : additional_neighbour_indices) {
        upper_left_actual_ids.push_back(n.id);
    }

    std::sort(upper_left_actual_ids.begin(), upper_left_actual_ids.end());
    std::sort(upper_left_expected_ids.begin(), upper_left_expected_ids.end());
    
    EXPECT_EQ(upper_left_actual_ids, upper_left_expected_ids) << "Upper left neighbor IDs do not match expected values";

    // Test upper right corner
    auto& upper_right_cell = container.cells[15];
    ASSERT_TRUE(upper_right_cell.is_halo) << "Upper right cell should be a halo cell";
    ASSERT_EQ(upper_right_cell.placement, Placement::TOP_RIGHT_CORNER) << "Cell should be top-right";

    // Get and verify neighbors for upper right particle
    auto upper_right_neighbours = container.get_neighbours(15);

    additional_neighbour_indices = container.get_periodic_neighbours(15);
    
    std::cout << "Number of neighbors for upper right particle: " << upper_right_neighbours.size() << std::endl;
    std::cout << "Upper right neighbor particle IDs: ";
    for(const auto& n : upper_right_neighbours) {
        std::cout << n->getId() << " ";
    }
    std::cout << std::endl;

    EXPECT_EQ(upper_right_neighbours.size() + additional_neighbour_indices.size(), 8) << "Upper right particle should have 8 neighbors with periodic conditions";

    std::vector<int> upper_right_expected_ids = {10, 11, 14, 15, 2, 3, 8, 12, 0};
    std::vector<int> upper_right_actual_ids;
    for (const auto& n : upper_right_neighbours) {
        upper_right_actual_ids.push_back(n->getId());
    }

    for (const auto& n : additional_neighbour_indices) {
        upper_right_actual_ids.push_back(n.id);
    }

    std::sort(upper_right_actual_ids.begin(), upper_right_actual_ids.end());
    std::sort(upper_right_expected_ids.begin(), upper_right_expected_ids.end());
    
    EXPECT_EQ(upper_right_actual_ids, upper_right_expected_ids) << "Upper right neighbor IDs do not match expected values";
}


TEST_F(PeriodicBoundaryTest, PeriodicTransitionTest) {
    // Create a particle at the left edge of the domain
    Particle left_edge_particle({0.1, 5.0, 0.0}, {-2.0, 0.0, 0.0}, 1.0, 0);
    container.insert(left_edge_particle,true);

    EXPECT_TRUE(container.size() == 1);
    EXPECT_EQ(container.cells[8].size(), 1);

    Calculation<Position>::run(container, 1, OPTIONS::LINKED_CELLS);

    // Check that all cells are empty
    for (const auto& cell : container.cells) {
        EXPECT_EQ(cell.size(), 0) << "Cell should be empty after particle transition";
    }

    EXPECT_EQ(container.particles_outbound.size(), 1);

    Calculation<BoundaryConditions>::run(container);

    EXPECT_EQ(container.particles_outbound.size(), 0);

    EXPECT_EQ(container.cells[11].size(), 1);

    EXPECT_EQ(container.size(), 1);

    //EXPECT_TRUE(container[0].is_periodic_copy);

    EXPECT_TRUE(!container[0].left_domain);

    EXPECT_TRUE(!container[0].outbound);

    Calculation<Position>::run(container, 1, OPTIONS::LINKED_CELLS);

    EXPECT_EQ(container.cells[11].size(), 0);

    EXPECT_EQ(container.cells[10].size(), 1);

    //EXPECT_TRUE(!container[0].is_periodic_copy);

    EXPECT_TRUE(!container[0].left_domain);

    EXPECT_TRUE(!container[0].outbound);
    // Check that all cells are empty    
}


TEST_F(PeriodicBoundaryTest, ParticleMovesFromCornerToCorner) {
    // Place a particle in the bottom-left corner
    Particle p({0.1, 0.1, 0.0}, {-1.0, -1.0, 0.0}, 1.0, 0);
    container.insert(p, true);

    EXPECT_EQ(container.cells[0].size(), 1);

    Calculation<Position>::run(container, 1, OPTIONS::LINKED_CELLS);

    // Check that the particle has moved to the top-right corner
    EXPECT_EQ(container.cells[0].size(), 0);

    Calculation<BoundaryConditions>::run(container);

    EXPECT_EQ(container.cells[container.cells.size() - 1].size(), 1);
}


TEST_F(PeriodicBoundaryTest, ReflectingTransitionTest) {
    // Create a particle at the left edge of the domain
    Particle left_edge_particle({0.1, 5.0, 0.0}, {-2.0, 0.0, 0.0}, 1.0, 0);
    b_container.insert(left_edge_particle,true);

    EXPECT_TRUE(b_container.size() == 1);
    EXPECT_EQ(b_container.cells[8].size(), 1);

    Calculation<Position>::run(b_container, 1, OPTIONS::LINKED_CELLS);

    // Check that all cells are empty
    for (const auto& cell : b_container.cells) {
        EXPECT_EQ(cell.size(), 0) << "Cell should be empty after particle transition";
    }

    EXPECT_EQ(b_container.particles_outbound.size(), 1);

    Calculation<BoundaryConditions>::run(b_container);

    EXPECT_EQ(b_container.particles_outbound.size(), 0);

    EXPECT_EQ(b_container.cells[11].size(), 1);

    EXPECT_EQ(b_container.size(), 1);

    EXPECT_TRUE(!b_container[0].left_domain);

    EXPECT_TRUE(!b_container[0].outbound);

    Calculation<Position>::run(b_container, 1, OPTIONS::LINKED_CELLS);

    EXPECT_EQ(b_container.cells[11].size(), 0);

    //EXPECT_EQ(b_container.cells[10].size(), 1);

    EXPECT_TRUE(!b_container[0].left_domain);

    EXPECT_TRUE(!b_container[0].outbound);
    // Check that all cells are empty    
}


TEST_F(PeriodicBoundaryTest, GhostParticlesTest) {
    // Insert a particle in the corner
    Particle corner_particle({0.1, 0.1, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0);
    container.insert(corner_particle, true);

    EXPECT_EQ(container.size(), 1) << "Container should have 1 particle initially";

    // Run boundary conditions to generate ghost particles
    Calculation<BoundaryConditions>::run(container);

    // Check for ghost particles
    size_t total_particles = 0;
    for (const auto& ghost : container.cell_ghost_particles_map) {
        total_particles++;
    }

    EXPECT_EQ(total_particles, 5) << "Container should have 8 ghost particles";
    /*
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(7)) << "Particle should be marked as left domain";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(3)) << "Particle should be marked as left domain";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(12)) << "Particle should be marked as left domain";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(13)) << "Particle should be marked as left domain";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(15)) << "Particle should be marked as left domain";
    */

    // Check ghost particle locations
    EXPECT_TRUE(container.cell_ghost_particles_map[3][0].position[0] == 10.1) << "Ghost particle in cell 7 should be to the right of original";
    EXPECT_TRUE(container.cell_ghost_particles_map[3][0].position[1] == 0.1) << "Ghost particle in cell 3 should be above original";

    EXPECT_TRUE(container.cell_ghost_particles_map[7][0].position[0] == 10.1) << "Ghost particle in cell 7 should be to the left of original";
    EXPECT_TRUE(container.cell_ghost_particles_map[7][0].position[1] == 0.1) << "Ghost particle in cell 7 should be below original";
    
    
    
    EXPECT_TRUE(container.cell_ghost_particles_map[12][0].position[0] == 0.1) << "Ghost particle in cell 12 should be to the left of original";
    EXPECT_TRUE(container.cell_ghost_particles_map[12][0].position[1] == 10.1) << "Ghost particle in cell 12 should be below original";
    
    EXPECT_TRUE(container.cell_ghost_particles_map[13][0].position[0] == 0.1) << "Ghost particle in cell 13 should be to the left and above original";
    EXPECT_TRUE(container.cell_ghost_particles_map[13][0].position[1] == 10.1) << "Ghost particle in cell 13 should be above original";


    EXPECT_TRUE(container.cell_ghost_particles_map[15][0].position[0] == 10.1) << "Ghost particle in cell 15 should be to the right and above original";
    EXPECT_TRUE(container.cell_ghost_particles_map[15][0].position[1] == 10.1) << "Ghost particle in cell 15 should be above original";
}

TEST_F(PeriodicBoundaryTest, LowerRightCornerGhostParticlesTest) {
    // Insert a particle in the lower right corner
    Particle corner_particle({9.9, 0.1, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0);
    container.insert(corner_particle, true);

    EXPECT_EQ(container.size(), 1) << "Container should have 1 particle initially";

    // Run boundary conditions to generate ghost particles
    Calculation<BoundaryConditions>::run(container);

    // Check for ghost particles
    size_t total_particles = 0;
    for (const auto& ghost : container.cell_ghost_particles_map) {
        total_particles++;
    }

    EXPECT_EQ(total_particles, 5) << "Container should have 5 ghost particles";
    
    /*
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(12)) << "Ghost particle should be in cell 0";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(0)) << "Ghost particle should be in cell 1";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(4)) << "Ghost particle should be in cell 4";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(14)) << "Ghost particle should be in cell 12";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(15)) << "Ghost particle should be in cell 13";
    */


    // Check for correct ghost particle locations
    EXPECT_NEAR(container.cell_ghost_particles_map[12][0].position[0], -0.1, 0.0001) << "Ghost particle in cell 12 should be at the same x-position as original";
    EXPECT_EQ(container.cell_ghost_particles_map[12][0].position[1], 10.1) << "Ghost particle in cell 12 should be above original";

    EXPECT_NEAR(container.cell_ghost_particles_map[0][0].position[0], -0.1, 0.0001) << "Ghost particle in cell 0 should be to the left of original";
    EXPECT_EQ(container.cell_ghost_particles_map[0][0].position[1], 0.1) << "Ghost particle in cell 0 should be at the same y-position as original";

    EXPECT_NEAR(container.cell_ghost_particles_map[4][0].position[0], -0.1, 0.0001) << "Ghost particle in cell 4 should be to the left of original";
    EXPECT_EQ(container.cell_ghost_particles_map[4][0].position[1], 0.1) << "Ghost particle in cell 4 should be above original";

    EXPECT_EQ(container.cell_ghost_particles_map[14][0].position[0], 9.9) << "Ghost particle in cell 14 should be at the same x-position as original";
    EXPECT_EQ(container.cell_ghost_particles_map[14][0].position[1], 10.1) << "Ghost particle in cell 14 should be below original";

    EXPECT_TRUE(container.cell_ghost_particles_map[15][0].position[0] == 9.9) << "Ghost particle in cell 15 should be to the left of original";
    EXPECT_TRUE(container.cell_ghost_particles_map[15][0].position[1] == 10.1) << "Ghost particle in cell 15 should be below original";

}

TEST_F(PeriodicBoundaryTest, UpperRightCornerGhostParticlesTest) {
    // Insert a particle in the upper right corner
    Particle corner_particle({9.9, 9.9, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0);
    container.insert(corner_particle, true);

    EXPECT_EQ(container.size(), 1) << "Container should have 1 particle initially";

    // Run boundary conditions to generate ghost particles
    Calculation<BoundaryConditions>::run(container);

    // Check for ghost particles
    size_t total_particles = 0;
    for (const auto& ghost : container.cell_ghost_particles_map) {
        total_particles++;
    }

    EXPECT_EQ(total_particles, 5) << "Container should have 5 ghost particles";

    /*EXPECT_TRUE(container.cell_ghost_particles_map.contains(0)) << "Ghost particle should be in cell 0";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(2)) << "Ghost particle should be in cell 1";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(3)) << "Ghost particle should be in cell 3";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(8)) << "Ghost particle should be in cell 4";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(12)) << "Ghost particle should be in cell 12";
    */

     // Check for correct ghost particle locations
    EXPECT_NEAR(container.cell_ghost_particles_map[12][0].position[0], -0.1, 0.0001) << "Ghost particle in cell 12 should be at the same x-position as original";
    EXPECT_NEAR(container.cell_ghost_particles_map[12][0].position[1], 9.9, 0.0001) << "Ghost particle in cell 12 should be above original";

    EXPECT_NEAR(container.cell_ghost_particles_map[0][0].position[0], -0.1, 0.0001) << "Ghost particle in cell 0 should be to the left of original";
    EXPECT_NEAR(container.cell_ghost_particles_map[0][0].position[1], -0.1, 0.0001) << "Ghost particle in cell 0 should be at the same y-position as original";

    EXPECT_NEAR(container.cell_ghost_particles_map[8][0].position[0], -0.1, 0.0001) << "Ghost particle in cell 4 should be to the left of original";
    EXPECT_NEAR(container.cell_ghost_particles_map[8][0].position[1], 9.9, 0.0001) << "Ghost particle in cell 4 should be above original";

    EXPECT_NEAR(container.cell_ghost_particles_map[2][0].position[0], 9.9, 0.0001) << "Ghost particle in cell 14 should be at the same x-position as original";
    EXPECT_NEAR(container.cell_ghost_particles_map[2][0].position[1], -0.1, 0.0001) << "Ghost particle in cell 14 should be below original";

    EXPECT_NEAR(container.cell_ghost_particles_map[3][0].position[0], 9.9, 0.0001) << "Ghost particle in cell 15 should be to the left of original";
    EXPECT_NEAR(container.cell_ghost_particles_map[3][0].position[1], -0.1, 0.0001) << "Ghost particle in cell 15 should be below original";
}


TEST_F(PeriodicBoundaryTest, UpperLeftCornerGhostParticlesTest) {
    // Insert a particle in the upper left corner
    Particle corner_particle({0.1, 9.9, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0);
    container.insert(corner_particle, true);

    EXPECT_EQ(container.size(), 1) << "Container should have 1 particle initially";

    // Run boundary conditions to generate ghost particles
    Calculation<BoundaryConditions>::run(container);

    // Check for ghost particles
    size_t total_particles = 0;
    for (const auto& ghost : container.cell_ghost_particles_map) {
        total_particles++;
    }

    EXPECT_EQ(total_particles, 5) << "Container should have 5 ghost particles";

    /*
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(15)) << "Ghost particle should be in cell 3";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(11)) << "Ghost particle should be in cell 7";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(3)) << "Ghost particle should be in cell 11";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(1)) << "Ghost particle should be in cell 15";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(0)) << "Ghost particle should be in cell 0";
    */

    EXPECT_NEAR(container.cell_ghost_particles_map[0][0].position[0], 0.1, 0.0001) << "Ghost particle in cell 12 should be at the same x-position as original";
    EXPECT_NEAR(container.cell_ghost_particles_map[0][0].position[1], -0.1, 0.0001) << "Ghost particle in cell 12 should be above original";

    EXPECT_NEAR(container.cell_ghost_particles_map[1][0].position[0], 0.1, 0.0001) << "Ghost particle in cell 0 should be to the left of original";
    EXPECT_NEAR(container.cell_ghost_particles_map[1][0].position[1], -0.1, 0.0001) << "Ghost particle in cell 0 should be at the same y-position as original";

    EXPECT_NEAR(container.cell_ghost_particles_map[15][0].position[0], 10.1, 0.0001) << "Ghost particle in cell 4 should be to the left of original";
    EXPECT_NEAR(container.cell_ghost_particles_map[15][0].position[1], 9.9, 0.0001) << "Ghost particle in cell 4 should be above original";

    EXPECT_NEAR(container.cell_ghost_particles_map[11][0].position[0], 10.1, 0.0001) << "Ghost particle in cell 14 should be at the same x-position as original";
    EXPECT_NEAR(container.cell_ghost_particles_map[11][0].position[1], 9.9, 0.0001) << "Ghost particle in cell 14 should be below original";

    EXPECT_NEAR(container.cell_ghost_particles_map[3][0].position[0], 10.1, 0.0001) << "Ghost particle in cell 15 should be to the left of original";
    EXPECT_NEAR(container.cell_ghost_particles_map[3][0].position[1], -0.1, 0.0001) << "Ghost particle in cell 15 should be below original";
}

TEST_F(PeriodicBoundaryTest, LeftSideGhostParticlesTest) {
    // Insert a particle on the left side, but not in the corner
    Particle left_side_particle({0.1, 5.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0);
    container.insert(left_side_particle, true);

    EXPECT_EQ(container.size(), 1) << "Container should have 1 particle initially";

    // Run boundary conditions to generate ghost particles
    Calculation<BoundaryConditions>::run(container);

    // Check for ghost particles
    size_t total_particles = 0;
    for (const auto& ghost : container.cell_ghost_particles_map) {
        total_particles += ghost.second.size();
    }

    EXPECT_EQ(total_particles, 3) << "Container should have 3 ghost particles";
    /*
    // Check if ghost particles are in the correct cells on the opposite side
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(15)) << "Ghost particle should be in right-side cell";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(7)) << "Ghost particle should be in right-side cell";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(11)) << "Ghost particle should be in right-side cell";
    */

    // Verify the positions of ghost particles
    for (const auto& [cell_index, ghosts] : container.cell_ghost_particles_map) {
        for (const auto& ghost : ghosts) {
            EXPECT_NEAR(ghost.position[0], 10.1, 1e-6) << "Ghost particle X position should be on the right side";
            EXPECT_NEAR(ghost.position[1], 5.0, 1e-6) << "Ghost particle Y position should match the original";
            EXPECT_NEAR(ghost.position[2], 0.0, 1e-6) << "Ghost particle Z position should match the original";
        }
    }
}

TEST_F(PeriodicBoundaryTest, RightSideGhostParticlesTest) {
    Particle right_side_particle({9.9, 5.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0);
    container.insert(right_side_particle, true);

    EXPECT_EQ(container.size(), 1) << "Container should have 1 particle initially";

    Calculation<BoundaryConditions>::run(container);

    size_t total_particles = 0;
    for (const auto& ghost : container.cell_ghost_particles_map) {
        total_particles += ghost.second.size();
    }

    EXPECT_TRUE(container.cells[11].size() == 1);

    EXPECT_EQ(total_particles, 3) << "Container should have 3 ghost particles";

    /*
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(12)) << "Ghost particle should be in left-side cell";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(4)) << "Ghost particle should be in left-side cell";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(8)) << "Ghost particle should be in left-side cell";
    */

    for (const auto& [cell_index, ghosts] : container.cell_ghost_particles_map) {
        for (const auto& ghost : ghosts) {
            EXPECT_NEAR(ghost.position[0], -0.1, 1e-6) << "Ghost particle X position should be on the left side";
            EXPECT_NEAR(ghost.position[1], 5.0, 1e-6) << "Ghost particle Y position should match the original";
            EXPECT_NEAR(ghost.position[2], 0.0, 1e-6) << "Ghost particle Z position should match the original";
        }
    }
}

TEST_F(PeriodicBoundaryTest, UpperSideGhostParticlesTest) {
    Particle upper_side_particle({5.0, 9.9, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0);
    container.insert(upper_side_particle, true);

    EXPECT_EQ(container.size(), 1) << "Container should have 1 particle initially";

    Calculation<BoundaryConditions>::run(container);

    size_t total_particles = 0;
    for (const auto& ghost : container.cell_ghost_particles_map) {
        total_particles += ghost.second.size();
    }

    EXPECT_EQ(total_particles, 3) << "Container should have 3 ghost particles";

    /*
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(1)) << "Ghost particle should be in bottom-side cell";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(2)) << "Ghost particle should be in bottom-side cell";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(3)) << "Ghost particle should be in bottom-side cell";
    */

    for (const auto& [cell_index, ghosts] : container.cell_ghost_particles_map) {
        for (const auto& ghost : ghosts) {
            EXPECT_NEAR(ghost.position[0], 5.0, 1e-6) << "Ghost particle X position should match the original";
            EXPECT_NEAR(ghost.position[1], -0.1, 1e-6) << "Ghost particle Y position should be on the bottom side";
            EXPECT_NEAR(ghost.position[2], 0.0, 1e-6) << "Ghost particle Z position should match the original";
        }
    }
}

TEST_F(PeriodicBoundaryTest, BottomSideGhostParticlesTest) {
    Particle bottom_side_particle({5.0, 0.1, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0);
    container.insert(bottom_side_particle, true);

    EXPECT_EQ(container.size(), 1) << "Container should have 1 particle initially";

    Calculation<BoundaryConditions>::run(container);

    size_t total_particles = 0;
    for (const auto& ghost : container.cell_ghost_particles_map) {
        total_particles += ghost.second.size();
    }

    EXPECT_EQ(total_particles, 3) << "Container should have 3 ghost particles";

    /*
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(13)) << "Ghost particle should be in upper-side cell";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(14)) << "Ghost particle should be in upper-side cell";
    EXPECT_TRUE(container.cell_ghost_particles_map.contains(15)) << "Ghost particle should be in upper-side cell";
    */

    for (const auto& [cell_index, ghosts] : container.cell_ghost_particles_map) {
        for (const auto& ghost : ghosts) {
            EXPECT_NEAR(ghost.position[0], 5.0, 1e-6) << "Ghost particle X position should match the original";
            EXPECT_NEAR(ghost.position[1], 10.1, 1e-6) << "Ghost particle Y position should be on the upper side";
            EXPECT_NEAR(ghost.position[2], 0.0, 1e-6) << "Ghost particle Z position should match the original";
        }
    }
}

TEST_F(PeriodicBoundaryTest, PreComputeNeighbouringCellsTest) {

    std::vector<size_t> expected_indices{4,5,1};
    auto& indices = container.cells[0].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {0,4,5,6,2};
    indices = container.cells[1].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {1,3,5,6,7};
    indices = container.cells[2].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {2,6,7};
    indices = container.cells[3].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {0,1,5,8,9};
    indices = container.cells[4].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);


    expected_indices = {0,1,2,4,6,8,9,10};
    indices = container.cells[5].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);


    expected_indices = {1,2,3,5,7,9,10,11};
    indices = container.cells[6].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {2,3,6,10,11};
    indices = container.cells[7].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {4,5,9,12,13};
    indices = container.cells[8].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {4,5,6,8,10,12,13,14};
    indices = container.cells[9].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {5,6,7,9,11,13,14,15};
    indices = container.cells[10].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {6,7,10,14,15};
    indices = container.cells[11].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {8,9,13};
    indices = container.cells[12].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {8,9,10,12,14};
    indices = container.cells[13].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {9,10,11,13,15};
    indices = container.cells[14].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {10,11,14};
    indices = container.cells[15].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);
}

TEST_F(PeriodicBoundaryTest, PreComputingThreeDimensionsTest) {

    ASSERT_EQ(d_container.x, 2);
    ASSERT_EQ(d_container.y, 2);
    ASSERT_EQ(d_container.z, 3);

    std::vector<size_t>expected_indices = {1,2,3,4,5,6,7};
    auto& indices = d_container.cells[0].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {0,2,3,4,5,6,7};
    indices = d_container.cells[1].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {0,1,3,4,5,6,7};
    indices = d_container.cells[2].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);

    expected_indices = {0,1,2,4,5,6,7};
    indices = d_container.cells[3].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);


    expected_indices = {0,1,2,4,5,6,7};
    indices = d_container.cells[3].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);


    expected_indices = {0,1,2,3,5,6,7,8,9,10,11};
    indices = d_container.cells[4].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);


    expected_indices = {0,1,2,3,4,6,7,8,9,10,11};
    indices = d_container.cells[5].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);


    expected_indices = {0,1,2,3,4,5,7,8,9,10,11};
    indices = d_container.cells[6].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);


    expected_indices = {0,1,2,3,4,5,6,8,9,10,11};
    indices = d_container.cells[7].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);


    expected_indices = {4,5,6,7,9,10,11};
    indices = d_container.cells[8].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);


    expected_indices = {4,5,6,7,8,10,11};
    indices = d_container.cells[9].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);


    expected_indices = {4,5,6,7,8,9,11};
    indices = d_container.cells[10].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);


    expected_indices = {4,5,6,7,8,9,10};
    indices = d_container.cells[11].neighbour_indices;
    std::sort(indices.begin(), indices.end());
    std::sort(expected_indices.begin(), expected_indices.end());
    ASSERT_EQ(expected_indices, indices);
}




