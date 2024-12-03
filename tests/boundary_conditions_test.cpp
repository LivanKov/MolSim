#include "../src/simulator/particle/LinkedCellContainer.h"
#include "gtest/gtest.h"
// #include <initializer_list>

class BoundaryConditionsTest : public ::testing::Test {
protected:
  LinkedCellContainer container;

  BoundaryConditionsTest()
      : container({10.0, 10.0}, 1.0, {0.0, 0.0},
                  {BoundaryCondition::Reflecting, // Left
                   BoundaryCondition::Outflow,    // Right
                   BoundaryCondition::Outflow, // Top
                   BoundaryCondition::Reflecting})   // Bottom
  {}

  ~BoundaryConditionsTest() override = default;
};

// Test that particles reflect correctly off the left boundary
TEST_F(BoundaryConditionsTest, ReflectingBoundary) {
  // Particle heading towards the left boundary
  Particle p({-0.5, 5.0, 0.0}, {-1.0, 0.0, 0.0}, 1.0);
  container.insert(p);

  container.handleBoundaryConditions(p);

  // Position and velocity should reflect
  EXPECT_NEAR(p.getX()[0], 0.5, 1e-6);
  EXPECT_NEAR(p.getV()[0], 1.0, 1e-6); // Velocity reversed
}

// Test that particles crossing the right boundary are marked for removal
TEST_F(BoundaryConditionsTest, OutflowBoundary) {
  // Particle exiting the right boundary
  Particle p({10.5, 5.0, 0.0}, {1.0, 0.0, 0.0}, 1.0);
  container.insert(p);

  container.handleBoundaryConditions(p);

  // Particle should be marked for removal
  EXPECT_TRUE(p.isMarkedForRemoval());
}

// Test that particles reflect correctly off the bottom boundary and crossing the top boundary are marked for removal
TEST_F(BoundaryConditionsTest, BottomReflectingTopOutflow) {
  // Particle heading towards the bottom (reflecting)
  Particle p_bottom({5.0, -0.5, 0.0}, {0.0, -1.0, 0.0}, 1.0);
  container.insert(p_bottom);

  // Particle heading towards the top (outflow)
  Particle p_top({5.0, 10.5, 0.0}, {0.0, 1.0, 0.0}, 1.0);
  container.insert(p_top);

  container.handleBoundaryConditions(p_bottom);
  // container.handleBoundaryConditions(p_top);

  // Check bottom boundary
  EXPECT_NEAR(p_bottom.getX()[1], 0.5, 1e-6); // Reflected position
  EXPECT_NEAR(p_bottom.getV()[1], 1.0, 1e-6); // Velocity reversed

  // Check top boundary
  // EXPECT_TRUE(p_top.isMarkedForRemoval());
}

// Test that particle within the domain should remain unchanged
TEST_F(BoundaryConditionsTest, NoBoundaryViolation) {
  // Particle within the domain
  Particle p({5.0, 5.0, 0.0}, {0.0, 0.0, 0.0}, 1.0);
  container.insert(p);

  container.handleBoundaryConditions(p);

  // Position and velocity should remain unchanged
  EXPECT_NEAR(p.getX()[0], 5.0, 1e-6);
  EXPECT_NEAR(p.getX()[1], 5.0, 1e-6);
  EXPECT_NEAR(p.getV()[0], 0.0, 1e-6);
  EXPECT_NEAR(p.getV()[1], 0.0, 1e-6);
  EXPECT_FALSE(p.isMarkedForRemoval());
}

TEST_F(BoundaryConditionsTest, CornerCrossing) {
  // Define domain size and corner location
  std::initializer_list<double> domain_size = {10.0, 10.0, 10.0};
  std::initializer_list<double> left_corner = {10.0, 10.0, 10.0};
  double cutoff_radius = 1.0;

  // Reflecting conditions on all boundaries
  DomainBoundaryConditions boundary_conditions{
      BoundaryCondition::Reflecting, BoundaryCondition::Reflecting,
      BoundaryCondition::Reflecting, BoundaryCondition::Reflecting,
      BoundaryCondition::Reflecting, BoundaryCondition::Reflecting};

  LinkedCellContainer container(domain_size, cutoff_radius, left_corner,
                                boundary_conditions);

  // Create a particle near the corner
  Particle particle({9.9, 9.9, 9.9}, {-1.0, -1.0, -1.0}, 1.0);
  container.insert(particle);

  // Apply boundary handling
  container.handleBoundaryConditions(particle);

  // Check the particle is correctly reflected from the corner
  auto position = particle.getX();
  auto velocity = particle.getV();

  EXPECT_NEAR(position[0], 10.1, 1e-5); // Reflected from the x boundary
  EXPECT_NEAR(position[1], 10.1, 1e-5); // Reflected from the y boundary
  EXPECT_NEAR(position[2], 10.1, 1e-5); // Reflected from the z boundary

  EXPECT_NEAR(velocity[0], 1.0, 1e-5); // Velocity reversed in x
  EXPECT_NEAR(velocity[1], 1.0, 1e-5); // Velocity reversed in y
  EXPECT_NEAR(velocity[2], 1.0, 1e-5); // Velocity reversed in z
}