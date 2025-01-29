#include "../src/simulator/particle/container/LinkedCellContainer.h"
#include "../src/simulator/calculations/BoundaryConditions.h"
#include "../src/simulator/calculations/Calculation.h"
#include "../src/simulator/calculations/Force.h"
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
  EXPECT_TRUE(container.cells[99].particle_ids.contains(0));
  EXPECT_TRUE(container.cells[99].placement == Placement::TOP_RIGHT_CORNER);

  EXPECT_EQ(position[0], 9.9); // Reflected from the x boundary
  EXPECT_EQ(position[1], 9.9); // Reflected from the y boundary
  EXPECT_EQ(position[2], 0); // Reflected from the z boundary
  
  EXPECT_EQ(velocity[0], -1.0); // Velocity reversed in x
  EXPECT_EQ(velocity[1], -1.0); // Velocity reversed in y
  EXPECT_EQ(velocity[2], 0.0); // Velocity reversed in z
}

// // verify that particles located at the left-front edge and right-back edge of a
// // 3D periodic boundary container correctly identify each other as ghost
// // neighbors using the method get_additional_neighbour_indices.
TEST_F(BoundaryConditionsTest, EdgeGhostNeighbours) {
  std::initializer_list<double> domain_size = {9.0, 9.0, 9.0};
  double cutoff_radius = 3.0;

  DomainBoundaryConditions boundary_conditions{
      BoundaryCondition::Periodic, BoundaryCondition::Periodic,
      BoundaryCondition::Periodic, BoundaryCondition::Periodic,
      BoundaryCondition::Periodic, BoundaryCondition::Periodic};

  LinkedCellContainer container_3d{};
  container_3d.initialize(domain_size, cutoff_radius, boundary_conditions);

  // Insert two particles: one at the left-front edge, one at the right-back
  // edge
  Particle left_front_particle({0.1, 5.0, 8.9}, {0.0, 0.0, 0.0}, 1.0, 0);
  Particle right_back_particle({8.9, 5.0, 0.1}, {0.0, 0.0, 0.0}, 1.0, 1);
  // Particle left_back_particle({0.1, 5.0, 0.1}, {0.0, 0.0, 0.0}, 1.0, 2);
  // Particle right_front_particle({8.9, 5.0, 8.9}, {0.0, 0.0, 0.0}, 1.0, 3);
  

  container_3d.insert(left_front_particle, true);
  container_3d.insert(right_back_particle, true);
  // container_3d.insert(left_back_particle, true);
  // container_3d.insert(right_front_particle, true);

  // Run boundary condition calculations to create ghost particles
  Calculation<BoundaryConditions>::run(container_3d);

  // Get ghost neighbors of the left-front particle
  auto ghost_neighbours_left_front =
      container_3d.get_additional_neighbour_indices(0);

  // Get ghost neighbors of the right-back particle
  auto ghost_neighbours_right_back =
      container_3d.get_additional_neighbour_indices(1);

  // Verify that the ghost neighbors are correctly identified
  bool found_right_back_as_neighbour = false;
  for (const auto &ghost : ghost_neighbours_left_front) {
    if (ghost.ptr == container_3d.cells_map[1]) {
      found_right_back_as_neighbour = true;
      break;
    }
  }

  bool found_left_front_as_neighbour = false;
  for (const auto &ghost : ghost_neighbours_right_back) {
    if (ghost.ptr == container_3d.cells_map[0]) {
      found_left_front_as_neighbour = true;
      break;
    }
  }

  // Assert that they are neighbors
  EXPECT_TRUE(found_right_back_as_neighbour)
      << "Right-back particle was not found as a ghost neighbor of left-front.";
  EXPECT_TRUE(found_left_front_as_neighbour)
      << "Left-front particle was not found as a ghost neighbor of right-back.";
}

// This test verifies that two particles located at opposite corners of a 3D
// periodic boundary container (i.e., BOTTOM_FRONT_RIGHT_CORNER and
// TOP_BACK_LEFT_CORNER) correctly identify each other as ghost neighbors
TEST_F(BoundaryConditionsTest, OppositeCornerGhostNeighbours) {
  std::initializer_list<double> domain_size = {10.0, 10.0, 10.0};
  double cutoff_radius = 2.0;

  // Periodic boundary conditions on all boundaries
  DomainBoundaryConditions boundary_conditions{
      BoundaryCondition::Periodic, BoundaryCondition::Periodic,
      BoundaryCondition::Periodic, BoundaryCondition::Periodic,
      BoundaryCondition::Periodic, BoundaryCondition::Periodic};

  // Initialize the 3D container
  LinkedCellContainer container_3d{};
  container_3d.initialize(domain_size, cutoff_radius, boundary_conditions);

  // Insert particles at opposite periodic corners
  Particle bottom_front_right({9.9, 0.1, 9.9}, {0.0, 0.0, 0.0}, 1.0, 0);
  Particle top_back_left({0.1, 9.9, 0.1}, {0.0, 0.0, 0.0}, 1.0, 1);
  container_3d.insert(bottom_front_right, true);
  container_3d.insert(top_back_left, true);

  // Apply boundary handling to generate ghost particles
  Calculation<BoundaryConditions>::run(container_3d);

  // Get ghost neighbors for both particles
  auto ghost_neighbours_bottom_front_right =
      container_3d.get_additional_neighbour_indices(0);
  auto ghost_neighbours_top_back_left =
      container_3d.get_additional_neighbour_indices(1);

  // Verify that each particle appears as a ghost neighbor in the other's list
  bool found_top_back_left_as_neighbour = false;
  for (const auto &ghost : ghost_neighbours_bottom_front_right) {
    if (ghost.ptr == container_3d.cells_map[1]) {
      found_top_back_left_as_neighbour = true;
      break;
    }
  }

  bool found_bottom_front_right_as_neighbour = false;
  for (const auto &ghost : ghost_neighbours_top_back_left) {
    if (ghost.ptr == container_3d.cells_map[0]) {
      found_bottom_front_right_as_neighbour = true;
      break;
    }
  }

  // Assertions
  EXPECT_TRUE(found_top_back_left_as_neighbour)
      << "Top-back-left particle was not found as a ghost neighbor of "
         "bottom-front-right.";
  EXPECT_TRUE(found_bottom_front_right_as_neighbour)
      << "Bottom-front-right particle was not found as a ghost neighbor of "
         "top-back-left.";
}

// Ensures that particles interact correctly across periodic boundaries and
// experience repulsive forces when close enough.
TEST_F(BoundaryConditionsTest, PeriodicRepulsiveForceZ) {
  // Define domain size (9x9x9) and cutoff radius (3.0)
  std::initializer_list<double> domain_size = {9.0, 9.0, 9.0};
  double cutoff_radius = 3.0;

  // Periodic boundary conditions on all sides
  DomainBoundaryConditions boundary_conditions{
      BoundaryCondition::Periodic, BoundaryCondition::Periodic,
      BoundaryCondition::Periodic, BoundaryCondition::Periodic,
      BoundaryCondition::Periodic, BoundaryCondition::Periodic};

  // Initialize the container
  LinkedCellContainer container_3d{};
  container_3d.initialize(domain_size, cutoff_radius, boundary_conditions);

  // Insert two particles at periodic opposite z-boundaries
  Particle particle_A({2.0, 4.0, 0.5}, {0.0, 0.0, 0.0}, 1.0, 0);
  Particle particle_B({2.0, 4.0, 8.5}, {0.0, 0.0, 0.0}, 1.0, 1);
  container_3d.insert(particle_A, true);
  container_3d.insert(particle_B, true);

  // Apply boundary conditions to create ghost particles
  Calculation<BoundaryConditions>::run(container_3d);

  // Apply force calculation (Lennard-Jones)
  Calculation<Force>::run(container_3d, ForceType::LENNARD_JONES,
                          OPTIONS::LINKED_CELLS);

  // Get forces on both particles
  auto force_A = container_3d.cells_map[0]->getF();
  auto force_B = container_3d.cells_map[1]->getF();

  // Check that the z-component of the force is repulsive
  EXPECT_GT(force_A[2], 0.0)
      << "Particle A should be pushed in the +z direction.";
  EXPECT_LT(force_B[2], 0.0)
      << "Particle B should be pushed in the -z direction.";

  // Verify that forces are equal and opposite (Newton's Third Law)
  EXPECT_NEAR(force_A[2], -force_B[2], 1e-6)
      << "Forces should be equal and opposite in z-direction.";
}
