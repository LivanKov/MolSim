#include <gtest/gtest.h>
#include "simulator/particle/LinkedCellContainer.h"
#include "simulator/particle/ParticleGenerator.h"
#include <array>


class LinkedCellTest : public ::testing::Test {
protected:
    LinkedCellTest() : container{{9.0,9.0}, 3.0, {0.0, 0.0, 0.0}} {}

    LinkedCellContainer container;
};

TEST_F(LinkedCellTest, LocationTest) {
    ParticleGenerator::insertCuboid(std::array<double, 3>{1.5, 1.5, 0.0}, std::array<size_t, 3>{3, 3, 3}, 3.0, 1.0, std::array<double,3>{0.0, 0.0, 0.0}, 0.0, container);

    //assert that every cells contains exactly one particle
    for(size_t i = 0; i < container.cells.size(); ++i){
        EXPECT_EQ(container.cells[i].size(), 1);
    }

}