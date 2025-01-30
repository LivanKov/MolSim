#include "../src/simulator/particle/container/LinkedCellContainer.h"
#include "../src/simulator/calculations/BoundaryConditions.h"
#include "../src/simulator/calculations/Calculation.h"
#include "../src/simulator/calculations/Force.h"
#include "../src/simulator/calculations/Position.h"
#include "../src/simulator/calculations/Velocity.h"

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
      container_3d.get_periodic_neighbours(0);

  // Get ghost neighbors of the right-back particle
  auto ghost_neighbours_right_back =
      container_3d.get_periodic_neighbours(1);

  // Verify that the ghost neighbors are correctly identified
  bool found_right_back_as_neighbour = false;
  for (const auto &ghost : ghost_neighbours_left_front) {
    if (ghost.ptr == container_3d.at(1)) {
      found_right_back_as_neighbour = true;
      break;
    }
  }

  bool found_left_front_as_neighbour = false;
  for (const auto &ghost : ghost_neighbours_right_back) {
    if (ghost.ptr == container_3d.at(0)) {
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
      container_3d.get_periodic_neighbours(0);
  auto ghost_neighbours_top_back_left =
      container_3d.get_periodic_neighbours(1);

  // Verify that each particle appears as a ghost neighbor in the other's list
  bool found_top_back_left_as_neighbour = false;
  for (const auto &ghost : ghost_neighbours_bottom_front_right) {
    if (ghost.ptr == container_3d.at(1)) {
      found_top_back_left_as_neighbour = true;
      break;
    }
  }

  bool found_bottom_front_right_as_neighbour = false;
  for (const auto &ghost : ghost_neighbours_top_back_left) {
    if (ghost.ptr == container_3d.at(0)) {
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
  std::initializer_list<double> domain_size = {10.8, 15, 10.8};
  double cutoff_radius = 3.6;

  DomainBoundaryConditions boundary_conditions{
      BoundaryCondition::Periodic, BoundaryCondition::Periodic,
      BoundaryCondition::Reflecting, BoundaryCondition::Reflecting,
      BoundaryCondition::Periodic, BoundaryCondition::Periodic};

  // Initialize the container
  LinkedCellContainer container_3d{};
  container_3d.initialize(domain_size, cutoff_radius, boundary_conditions);

  // Insert two particles at periodic opposite z-boundaries
  Particle particle_A({1.8, 3.75, 0.2}, {0.0, 0.0, 0.0}, 1.0, 0, 1.0, 1.2);
  Particle particle_B({1.8, 3.75, 10.6}, {0.0, 0.0, 0.0}, 1.0, 1, 1.0, 1.2);

  Particle particle_C({5.4, 3.75, 0.2}, {0.0, 0.0, 0.0}, 1.0, 2, 1.0, 1.2);
  Particle particle_D({5.4, 3.75, 10.6}, {0.0, 0.0, 0.0}, 1.0, 3, 1.0, 1.2);

  Particle particle_E({9, 3.75, 0.2}, {0.0, 0.0, 0.0}, 1.0, 4, 1.0, 1.2);
  Particle particle_F({9, 3.75, 10.6}, {0.0, 0.0, 0.0}, 1.0, 5, 1.0, 1.2);

  Particle particle_G({0.2, 3.75, 1.8}, {0.0, 0.0, 0.0}, 1.0, 6, 1.0, 1.2);
  Particle particle_H({10.6, 3.75, 1.8}, {0.0, 0.0, 0.0}, 1.0, 7, 1.0, 1.2);

  Particle particle_I({0.2, 3.75, 5.4}, {0.0, 0.0, 0.0}, 1.0, 8, 1.0, 1.2);
  Particle particle_J({10.6, 3.75, 5.4}, {0.0, 0.0, 0.0}, 1.0, 9, 1.0, 1.2);

  Particle particle_K({0.2, 3.75, 9}, {0.0, 0.0, 0.0}, 1.0, 10, 1.0, 1.2);
  Particle particle_L({10.6, 3.75, 9}, {0.0, 0.0, 0.0}, 1.0, 11, 1.0, 1.2);

  container_3d.insert(particle_A, true);
  container_3d.insert(particle_B, true);
  container_3d.insert(particle_C, true);
  container_3d.insert(particle_D, true);
  container_3d.insert(particle_E, true);
  container_3d.insert(particle_F, true);
  container_3d.insert(particle_G, true);
  container_3d.insert(particle_H, true);
  container_3d.insert(particle_I, true);
  container_3d.insert(particle_J, true);
  container_3d.insert(particle_K, true);
  container_3d.insert(particle_L, true);


  Calculation<BoundaryConditions>::run(container_3d);

  Calculation<Force>::run(container_3d, ForceType::LENNARD_JONES,
                          OPTIONS::LINKED_CELLS);

  // Get forces on both particles
  auto force_A = container_3d.at(0)->getF();
  auto force_B = container_3d.at(1)->getF();
  auto force_C = container_3d.at(2)->getF();
  auto force_D = container_3d.at(3)->getF();
  auto force_E = container_3d.at(4)->getF();
  auto force_F = container_3d.at(5)->getF();
  auto force_G = container_3d.at(6)->getF();
  auto force_H = container_3d.at(7)->getF();
  auto force_I = container_3d.at(8)->getF();
  auto force_J = container_3d.at(9)->getF();
  auto force_K = container_3d.at(10)->getF();
  auto force_L = container_3d.at(11)->getF();

  // Check that the z-component of the force is repulsive
  EXPECT_GT(force_A[2], 0.0)
      << "Particle A should be pushed in the +z direction.";
  EXPECT_LT(force_B[2], 0.0)
      << "Particle B should be pushed in the -z direction.";
  EXPECT_GT(force_C[2], 0.0)
      << "Particle C should be pushed in the +z direction.";
  EXPECT_LT(force_D[2], 0.0)
      << "Particle D should be pushed in the -z direction.";
  EXPECT_GT(force_E[2], 0.0)
      << "Particle E should be pushed in the +z direction.";
  EXPECT_LT(force_F[2], 0.0)
      << "Particle F should be pushed in the -z direction.";

  EXPECT_GT(force_G[0], 0.0)
      << "Particle G should be pushed in the +x direction.";
  EXPECT_LT(force_H[0], 0.0)
      << "Particle H should be pushed in the -x direction.";
  EXPECT_GT(force_I[0], 0.0)
      << "Particle I should be pushed in the +x direction.";
  EXPECT_LT(force_J[0], 0.0)
      << "Particle J should be pushed in the -x direction.";
  EXPECT_GT(force_K[0], 0.0)
      << "Particle K should be pushed in the +x direction.";
  EXPECT_LT(force_L[0], 0.0)
      << "Particle L should be pushed in the -x direction.";

  // Verify that forces are equal and opposite (Newton's Third Law)
  EXPECT_NEAR(force_A[2], -force_B[2], 1e-6)
      << "Forces should be equal and opposite in z-direction.";
  EXPECT_NEAR(force_C[2], -force_D[2], 1e-6)
      << "Forces should be equal and opposite in z-direction.";
  EXPECT_NEAR(force_E[2], -force_F[2], 1e-6)
      << "Forces should be equal and opposite in z-direction.";
  EXPECT_NEAR(force_G[0], -force_H[0], 1e-6)
      << "Forces should be equal and opposite in z-direction.";
  EXPECT_NEAR(force_I[0], -force_J[0], 1e-6)
      << "Forces should be equal and opposite in z-direction.";
  EXPECT_NEAR(force_K[0], -force_L[0], 1e-6)
      << "Forces should be equal and opposite in z-direction.";
}
