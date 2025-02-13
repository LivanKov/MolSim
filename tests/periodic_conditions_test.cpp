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


//Main test fixture, initialize multiple containers with different boundary conditions
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

//Check that the ghost particles are correctly identified as such
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

    std::vector<int> bottom_right_expected_ids = {2, 6, 7, 0, 4, 14, 15, 12};
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

    std::vector<int> upper_left_expected_ids = {8, 9, 13, 0, 1, 11, 15, 3};
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

    std::vector<int> upper_right_expected_ids = {10, 11, 14, 2, 3, 8, 12, 0};
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

//Check that the periodic transition from one side of the domain to the other is correctly handled
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

//Check that a particle that moves from one corner to the other stays in the correct cell. Edge case
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

//Check the behaviour of a domain that incorporates reflecting and periodic boundary conditions
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

//Check for single ghost particles being created
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

//Check for single ghost particles being created in the lower right corner
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

//Check for single ghost particles being created in the upper right corner
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

//Check that ghost particles are created correctly when a particle is inserted on the left side of the domain
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

//Check that ghost particles are created correctly when a particle is inserted on the right side of the domain
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

    for (const auto& [cell_index, ghosts] : container.cell_ghost_particles_map) {
        for (const auto& ghost : ghosts) {
            EXPECT_NEAR(ghost.position[0], -0.1, 1e-6) << "Ghost particle X position should be on the left side";
            EXPECT_NEAR(ghost.position[1], 5.0, 1e-6) << "Ghost particle Y position should match the original";
            EXPECT_NEAR(ghost.position[2], 0.0, 1e-6) << "Ghost particle Z position should match the original";
        }
    }
}

//Check that ghost particles are created correctly when a particle is inserted on the upper side of the domain
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

    for (const auto& [cell_index, ghosts] : container.cell_ghost_particles_map) {
        for (const auto& ghost : ghosts) {
            EXPECT_NEAR(ghost.position[0], 5.0, 1e-6) << "Ghost particle X position should match the original";
            EXPECT_NEAR(ghost.position[1], -0.1, 1e-6) << "Ghost particle Y position should be on the bottom side";
            EXPECT_NEAR(ghost.position[2], 0.0, 1e-6) << "Ghost particle Z position should match the original";
        }
    }
}

//Check that ghost particles are created correctly when a particle is inserted on the bottom side of the domain
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

// Verifies that a particle initially at (8.9, 8.9, 9.3) (Cell 27) in a 9×9×9
// periodic boundary container with a cutoff radius of 3.0 is correctly placed
// into Cell 8 after applying periodic boundary conditions. (potentially out of
// bound)
TEST_F(PeriodicBoundaryTest, PeriodicCellWrapping) {
  std::initializer_list<double> domain_siz = {9.0, 9.0, 9.0};
  double cutoff_radius = 3.0;
  DomainBoundaryConditions boundary_condition{
      BoundaryCondition::Periodic, BoundaryCondition::Periodic,
      BoundaryCondition::Periodic, BoundaryCondition::Periodic,
      BoundaryCondition::Periodic, BoundaryCondition::Periodic};

  LinkedCellContainer container_3d{};
  container_3d.initialize(domain_siz, cutoff_radius, boundary_condition);

  Particle particle({8.9, 8.9, 8.9}, {0.0, 0.0, 2.0}, 1.0, 0);
  container_3d.insert(particle, true);
  Calculation<Position>::run(container_3d, 1, OPTIONS::LINKED_CELLS);

  Calculation<BoundaryConditions>::run(container_3d);
  auto &pt = container_3d[0];

  std::array<double, 3> expected_position = {8.9, 8.9, 1.9};
  auto new_position = pt.getX();
  EXPECT_NEAR(new_position[0], expected_position[0], 1e-6);
  EXPECT_NEAR(new_position[1], expected_position[1], 1e-6);
  EXPECT_NEAR(new_position[2], expected_position[2], 1e-6);

  size_t expected_cell_index = container_3d.get_cell_index(expected_position);
  size_t new_cell_index = container_3d.get_cell_index(new_position);

  EXPECT_EQ(new_cell_index, 8) << "Particle should be placed in Cell 8 after "
                                  "periodic boundary handling.";
}
