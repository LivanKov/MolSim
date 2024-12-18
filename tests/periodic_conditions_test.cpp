class PeriodicBoundaryTest : public ::testing::Test {
protected:
    LinkedCellContainer container;

    PeriodicBoundaryTest() : container{} {
        container.initialize({10.0, 10.0}, 1.0,
            {BoundaryCondition::Periodic, // Left
             BoundaryCondition::Periodic, // Right
             BoundaryCondition::Periodic, // Top
             BoundaryCondition::Periodic  // Bottom
            });
    }
};


