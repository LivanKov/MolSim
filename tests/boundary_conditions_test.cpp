#include "../src/simulator/particle/container/LinkedCellContainer.h"
#include "../src/simulator/calculations/BoundaryConditions.h"
#include "../src/simulator/calculations/Calculation.h"
#include "gtest/gtest.h"
// #include <initializer_list>

class BoundaryConditionsTest : public ::testing::Test {
protected:
  LinkedCellContainer container;

  BoundaryConditionsTest()
      : container{}   
  {
    container.initialize({10.0, 10.0}, 1.0,
                  {BoundaryCondition::Reflecting, // Left
                   BoundaryCondition::Outflow,    // Right
                   BoundaryCondition::Outflow, // Top
                   BoundaryCondition::Reflecting // Bottom 
                   });
  }

  ~BoundaryConditionsTest() override = default;
};

// Test that particles reflect correctly off the left boundary
TEST_F(BoundaryConditionsTest, ReflectingBoundary) {
  // Particle heading towards the left boundary
  Particle p({0.5, 1.5, 0.0}, {-1.0, 0.0, 0.0}, 1.0, 0);
  container.insert(p, true);

  Calculation<BoundaryConditions>::run(container);

  auto& cell = container.cells[10];
  EXPECT_TRUE(cell.is_halo);
  EXPECT_TRUE(container.cells[10].size() == 1);
  EXPECT_TRUE(container.x == 10);
  EXPECT_EQ(container.boundary_conditions_.left, BoundaryCondition::Reflecting);

  auto p_ = container[0];
  // Position and velocity should reflect
  EXPECT_EQ(p_.getX()[0], 0.5);
  EXPECT_EQ(p_.getV()[0], 1.0); // Position unchanged
  EXPECT_EQ(p_.getV()[1], 0.0); // Velocity reversed
  EXPECT_EQ(p_.getV()[2], 0.0);
}

// Test that particles crossing the right boundary are marked for removal
TEST_F(BoundaryConditionsTest, OutflowBoundary) {
  // Particle exiting the right boundary
  Particle p({10.5, 5.0, 0.0}, {1.0, 0.0, 0.0}, 1.0, 0);
  container.insert(p, true);

  Calculation<BoundaryConditions>::run(container);

  auto& p_ = container[0];

  // Particle should be marked for removal
  EXPECT_TRUE(p_.left_domain);
  for(auto& cell : container.cells){
    EXPECT_TRUE(cell.size() == 0);
  }
}

// Test that particles reflect correctly off the bottom boundary and crossing the top boundary are marked for removal
TEST_F(BoundaryConditionsTest, BottomReflectingTopOutflow) {
  // Particle heading towards the bottom (reflecting)
  Particle p_bottom({5.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, 1.0, 0);
  container.insert(p_bottom, true);

  // Particle heading towards the top (outflow)
  Particle p_top({5.0, 10.5, 0.0}, {0.0, 1.0, 0.0}, 1.0, 1);
  container.insert(p_top, true);

  Calculation<BoundaryConditions>::run(container);

  auto& p_1 = container[0];
  auto& p_2 = container[1];

  for(size_t i = 0; i < container.cells.size(); i++){
    auto& cell = container.cells[i];
    if(i == 5){
      EXPECT_TRUE(cell.size() == 1);
    } else {
      EXPECT_TRUE(cell.size() == 0);
    }
  }

  // Check bottom boundary
  EXPECT_EQ(p_1.getX()[0], 5.0);
  EXPECT_EQ(p_1.getV()[1], 1.0);

  // Check top boundary
  EXPECT_TRUE(p_2.left_domain);
}

// Test that particle within the domain should remain unchanged
TEST_F(BoundaryConditionsTest, NoBoundaryViolation) {
  // Particle within the domain
  Particle p({5.0, 5.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0);
  container.insert(p, true);

  Calculation<BoundaryConditions>::run(container);

  for(size_t i = 0; i < container.cells.size(); i++){
    auto& cell = container.cells[i];
    if(i == 55){
      EXPECT_TRUE(cell.size() == 1);
    } else {
      EXPECT_TRUE(cell.size() == 0);
    }
  }

  auto& pt = container[0];

  // Position and velocity should remain unchanged
  EXPECT_NEAR(pt.getX()[0], 5.0, 1e-6);
  EXPECT_NEAR(pt.getX()[1], 5.0, 1e-6);
  EXPECT_NEAR(pt.getV()[0], 0.0, 1e-6);
  EXPECT_NEAR(pt.getV()[1], 0.0, 1e-6);
  EXPECT_FALSE(pt.left_domain);
}

TEST_F(BoundaryConditionsTest, CornerCrossing) {
  // Define domain size and corner location
  std::initializer_list<double> domain_size = {10.0, 10.0};
  double cutoff_radius = 1.0;

  // Reflecting conditions on all boundaries
  DomainBoundaryConditions boundary_conditions{
      BoundaryCondition::Reflecting, BoundaryCondition::Reflecting,
      BoundaryCondition::Reflecting, BoundaryCondition::Reflecting,
      BoundaryCondition::Reflecting, BoundaryCondition::Reflecting};

  LinkedCellContainer container{};
  container.initialize(domain_size, cutoff_radius,
                                boundary_conditions);

  EXPECT_EQ(container.cells[99].size(), 0);

  // Create a particle near the corner
  Particle particle({9.9, 9.9, 0}, {1.0, 1.0, 0.0}, 1.0, 0);
  container.insert(particle,true);

  // Apply boundary handling
  Calculation<BoundaryConditions>::run(container);


  auto& pt = container[0];


  // Check the particle is correctly reflected from the corner
  auto position = pt.getX();
  auto velocity = pt.getV();

  EXPECT_EQ(container.cells[99].size(), 1);
  EXPECT_TRUE(std::find(container.cells[99].particle_ids.begin(), container.cells[99].particle_ids.end(), 0) != container.cells[99].particle_ids.end());
  EXPECT_TRUE(container.cells[99].placement == Placement::TOP_RIGHT_CORNER);

  EXPECT_EQ(position[0], 9.9); // Reflected from the x boundary
  EXPECT_EQ(position[1], 9.9); // Reflected from the y boundary
  EXPECT_EQ(position[2], 0); // Reflected from the z boundary
  
  EXPECT_EQ(velocity[0], -1.0); // Velocity reversed in x
  EXPECT_EQ(velocity[1], -1.0); // Velocity reversed in y
  EXPECT_EQ(velocity[2], 0.0); // Velocity reversed in z
}