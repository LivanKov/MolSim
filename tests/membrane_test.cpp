#include <gtest/gtest.h>
#include "simulator/particle/container/LinkedCellContainer.h"
#include "simulator/particle/ParticleGenerator.h"
#include "io/input/cli/SimParams.h"

class MembraneTest : public ::testing::Test {
protected:
    LinkedCellContainer container;

    MembraneTest() : container{} {
        container.initialize({10.0, 10.0}, 2.5,
            {BoundaryCondition::Outflow, // Left
             BoundaryCondition::Outflow, // Right
             BoundaryCondition::Outflow, // Top
             BoundaryCondition::Outflow  // Bottom
            });
    }
};

TEST_F(MembraneTest, SimpleMembraneTest) {
    // Test parameters
    SimParams::fixed_Domain = false;
    ParticleGenerator::insertCuboid(
      std::array<double, 3>{0.0, 0.0, 0.0}, std::array<size_t, 3>{3, 3, 1}, 3.0,
      1.0, std::array<double, 3>{0.0, 0.0, 0.0}, container, 1.0, 1.0,true);

    ASSERT_TRUE(container.size() == 9);
    
    EXPECT_EQ(container[0].membrane_neighbours.size(), 2);
    EXPECT_EQ(container[0].diagonal_membrane_neighbours.size(), 1);
    EXPECT_EQ(container[1].membrane_neighbours.size(), 3);
    EXPECT_EQ(container[1].diagonal_membrane_neighbours.size(), 2);
   
    EXPECT_EQ(container[2].membrane_neighbours.size(), 2);
    EXPECT_EQ(container[2].diagonal_membrane_neighbours.size(), 1);

    EXPECT_EQ(container[3].membrane_neighbours.size(), 3);
    EXPECT_EQ(container[3].diagonal_membrane_neighbours.size(), 2);

    EXPECT_EQ(container[4].membrane_neighbours.size(), 4);
    EXPECT_EQ(container[4].diagonal_membrane_neighbours.size(), 4);

    EXPECT_EQ(container[5].membrane_neighbours.size(), 3);
    EXPECT_EQ(container[5].diagonal_membrane_neighbours.size(), 2);

    EXPECT_EQ(container[6].membrane_neighbours.size(), 2);
    EXPECT_EQ(container[6].diagonal_membrane_neighbours.size(), 1);

    EXPECT_EQ(container[7].membrane_neighbours.size(), 3);
    EXPECT_EQ(container[7].diagonal_membrane_neighbours.size(), 2);

    EXPECT_EQ(container[8].membrane_neighbours.size(), 2);
    EXPECT_EQ(container[8].diagonal_membrane_neighbours.size(), 1);


    //check individual membrane neighbours

    for (size_t particle = 0; particle < container.size(); ++particle) {
        std::vector<int> expected_ids;
        std::vector<int> actual_ids;
        
        // Get row and column from particle index
        int row = particle / 3;
        int col = particle % 3;
        
        // Check all 8 surrounding positions
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                if (dx == 0 && dy == 0) continue;
                
                int new_row = row + dy;
                int new_col = col + dx;
                
                // Check if neighbor position is valid
                if (new_row >= 0 && new_row < 3 && new_col >= 0 && new_col < 3) {
                    int neighbor = new_row * 3 + new_col;
                    expected_ids.push_back(neighbor);
                }
            }
        }
        
        for (const auto& neighbor : container[particle].membrane_neighbours) {
            actual_ids.push_back(neighbor->getId());
        }
        for (const auto& neighbor : container[particle].diagonal_membrane_neighbours) {
            actual_ids.push_back(neighbor->getId());
        }
        
        std::sort(expected_ids.begin(), expected_ids.end());
        std::sort(actual_ids.begin(), actual_ids.end());
        
        EXPECT_EQ(actual_ids, expected_ids) << "Mismatch for particle " << particle;
    }
    

    

    









}