#include "simulator/particle/ParticleGenerator.h"
#include "simulator/particle/container/LinkedCellContainer.h"
#include "io/input/cli/SimParams.h"
#include "utils/logger/Logger.h"
#include "simulator/calculations/BoundaryConditions.h"
#include <array>
#include <cmath>
#include <gtest/gtest.h>


class PeriodicBoundaryTest : public ::testing::Test {
protected:
    LinkedCellContainer container;

    PeriodicBoundaryTest() : container{} {
        container.initialize({10.0, 10.0}, 2.5,
            {BoundaryCondition::Periodic, // Left
             BoundaryCondition::Periodic, // Right
             BoundaryCondition::Periodic, // Top
             BoundaryCondition::Periodic  // Bottom
            });
    }
};


TEST_F(PeriodicBoundaryTest, BasicNeighbourTest) {
    SimParams::fixed_Domain = false;
    ParticleGenerator::insertCuboid(
      std::array<double, 3>{0.0, 0.0, 0.0}, std::array<size_t, 3>{4, 4, 1}, 2.5,
      1.0, std::array<double, 3>{0.0, 0.0, 0.0}, container);

    // Test basic container setup
    ASSERT_EQ(container.size(), 16) << "Container should have 16 particles";
    ASSERT_EQ(container.halo_cell_indices.size(), 12) << "Should have 12 halo cells";

    // Verify cell contents
    for(size_t i = 0; i < container.cells.size(); i++){
        auto& cell = container.cells[i];
        ASSERT_EQ(cell.size(), 1) << "Cell " << i << " should contain exactly 1 particle";
        ASSERT_TRUE(cell.particle_ids.contains(i)) << "Cell " << i << " should contain particle " << i;
    }

    // Test corner cell properties
    auto& corner_cell = container.cells[0];
    ASSERT_TRUE(corner_cell.is_halo) << "Corner cell should be a halo cell";
    ASSERT_EQ(corner_cell.placement, Placement::BOTTOM_LEFT_CORNER) << "Corner cell should be bottom-left";

    // Get and verify neighbors for corner particle
    auto neighbours = container.get_neighbours(0);
    
    // Print debug info
    std::cout << "Number of neighbors found: " << neighbours.size() << std::endl;
    std::cout << "Neighbor particle IDs: ";
    for(const auto& n : neighbours) {
        std::cout << n->getType() << " ";
    }
    std::cout << std::endl;

    // TODO: Uncomment and adjust expected neighbor count once verified
    EXPECT_EQ(neighbours.size(), 9) << "Corner particle should have 8 neighbors with periodic conditions";

    std::vector<int> expected_ids = {0, 1, 4, 5, 12, 13, 3, 7, 15};
    std::vector<int> actual_ids;
    for (const auto& n : neighbours) {
        actual_ids.push_back(n->getType());
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
    
    std::cout << "Number of neighbors for bottom right particle: " << bottom_right_neighbours.size() << std::endl;
    std::cout << "Bottom right neighbor particle IDs: ";
    for(const auto& n : bottom_right_neighbours) {
        std::cout << n->getType() << " ";
    }
    std::cout << std::endl;

    EXPECT_EQ(bottom_right_neighbours.size(), 9) << "Bottom right particle should have 9 neighbors with periodic conditions";

    std::vector<int> bottom_right_expected_ids = {2, 3, 6, 7, 0, 4, 14, 15, 12};
    std::vector<int> bottom_right_actual_ids;
    for (const auto& n : bottom_right_neighbours) {
        bottom_right_actual_ids.push_back(n->getType());
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
    
    std::cout << "Number of neighbors for upper left particle: " << upper_left_neighbours.size() << std::endl;
    std::cout << "Upper left neighbor particle IDs: ";
    for(const auto& n : upper_left_neighbours) {
        std::cout << n->getType() << " ";
    }
    std::cout << std::endl;

    EXPECT_EQ(upper_left_neighbours.size(), 9) << "Upper left particle should have 9 neighbors with periodic conditions";

    std::vector<int> upper_left_expected_ids = {8, 9, 12, 13, 0, 1, 11, 15, 3};
    std::vector<int> upper_left_actual_ids;
    for (const auto& n : upper_left_neighbours) {
        upper_left_actual_ids.push_back(n->getType());
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
    
    std::cout << "Number of neighbors for upper right particle: " << upper_right_neighbours.size() << std::endl;
    std::cout << "Upper right neighbor particle IDs: ";
    for(const auto& n : upper_right_neighbours) {
        std::cout << n->getType() << " ";
    }
    std::cout << std::endl;

    EXPECT_EQ(upper_right_neighbours.size(), 9) << "Upper right particle should have 9 neighbors with periodic conditions";

    std::vector<int> upper_right_expected_ids = {10, 11, 14, 15, 2, 3, 8, 12, 0};
    std::vector<int> upper_right_actual_ids;
    for (const auto& n : upper_right_neighbours) {
        upper_right_actual_ids.push_back(n->getType());
    }
    std::sort(upper_right_actual_ids.begin(), upper_right_actual_ids.end());
    std::sort(upper_right_expected_ids.begin(), upper_right_expected_ids.end());
    
    EXPECT_EQ(upper_right_actual_ids, upper_right_expected_ids) << "Upper right neighbor IDs do not match expected values";

}
