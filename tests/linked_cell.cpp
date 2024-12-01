#include <gtest/gtest.h>
#include "simulator/particle/LinkedCellContainer.h"


class LinkedCellTest : public ::testing::Test {

    
    LinkedCellTest() : container{std::vector<double>{9.0,9.0}, 3.0, std::array<double,3>{0.0, 0.0, 0.0} } {}

    LinkedCellContainer container;

};