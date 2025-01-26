#include "simulator/particle/ParticleGenerator.h"
#include "simulator/particle/container/LinkedCellContainer.h"
#include "utils/logger/Logger.h"
#include <array>
#include <cmath>
#include <gtest/gtest.h>

class LinkedCellTest : public ::testing::Test {
protected:
  LinkedCellTest() : container{}, container_3d{}, uneven_container{} {
    DomainBoundaryConditions boundary_conditions{
        BoundaryCondition::Outflow, BoundaryCondition::Outflow,
        BoundaryCondition::Outflow, BoundaryCondition::Outflow,
        BoundaryCondition::Outflow, BoundaryCondition::Outflow};
    container.initialize({9.0, 9.0}, 3.0, boundary_conditions);
    container_3d.initialize({9.0, 9.0, 9.0}, 3.0, boundary_conditions);
    uneven_container.initialize({9.0, 9.0, 8.0}, 2.0, boundary_conditions);
  }

  Logger &logger = Logger::getInstance("debug");
  LinkedCellContainer container;
  LinkedCellContainer container_3d;
  LinkedCellContainer uneven_container;
};

TEST_F(LinkedCellTest, LocationTest) {
  for (size_t i = 0; i < container.cells.size(); ++i) {
    EXPECT_EQ(container.cells[i].size(), 0);
  }

  EXPECT_TRUE(container.cells.size() == 9);

  // single particles insertion test

  Particle p1(std::array<double, 3>{5.0, 4.0, 0.0},
              std::array<double, 3>{0.0, 0.0, 0.0}, 1.0, 0);

  container.insert(p1, true);
  EXPECT_TRUE(container.domain_size_.size() == 3);

  EXPECT_TRUE(container.is_within_domain(p1.getX()));

  EXPECT_EQ(container.cells[4].size(), 1);

  for (size_t i = 0; i < container.cells.size(); ++i) {
    if (i != 4) {
      EXPECT_EQ(container.cells[i].size(), 0);
    }
  }

  // update particle location, make sure it is removed from old cell and
  // inserted into new cell

  std::array<double, 3> old_position = p1.getX();

    container[0].updateX(7.0, 7.0, 0.0);

    container.update_particle_location(0, old_position);

    EXPECT_EQ(container.cells[8].size(), 1);

  for (size_t i = 0; i < container.cells.size(); ++i) {
    if (i != 8) {
      EXPECT_EQ(container.cells[i].size(), 0);
    }
  }

  // leave the domain

    std::array<double, 3>another_old_position = container[0].getX();
    container[0].updateX(10.0, 10.0, 0.0);
    container.update_particle_location(0, another_old_position);

  // Ensure that the particle is still within the container

  EXPECT_TRUE(container.size() == 1);

  // Ensure that the corresponding flag has been set

  EXPECT_TRUE(std::find(container.particles_outbound.begin(),
                        container.particles_outbound.end(), 0) !=
              container.particles_outbound.end());

  // Ensure that the particle is no longer in the cells

  for (size_t i = 0; i < container.cells.size(); ++i) {
    EXPECT_EQ(container.cells[i].size(), 0);
  }
}

TEST_F(LinkedCellTest, CuboidTest) {

  // Insert a 3x3x1 cuboid of particles with a side length of 3.0 and a mass
  // of 1.0
  ParticleGenerator::insertCuboid(
      std::array<double, 3>{1.5, 1.5, 0.0}, std::array<size_t, 3>{3, 3,
      1}, 3.0, 1.0, std::array<double, 3>{0.0, 0.0, 0.0}, container);

  EXPECT_TRUE(container.size() == 9);

  for (size_t i = 0; i < container.cells.size(); ++i) {
    EXPECT_EQ(container.cells[i].size(), 1);
  }

  container.clear();

  // Insert a 3x3x3 cuboid of particles with a side length of 3.0 and a mass
  // of 1.0

  ParticleGenerator::insertCuboid(
      std::array<double, 3>{1.5, 1.5, 0.0}, std::array<size_t, 3>{3, 3,
      3}, 3.0, 1.0, std::array<double, 3>{0.0, 0.0, 0.0}, container);

  EXPECT_TRUE(container.size() == 27);

  for (size_t i = 0; i < container.cells.size(); ++i) {
    EXPECT_EQ(container.cells[i].size(), 1);
  }
}

TEST_F(LinkedCellTest, NeighbourTest) {

  ParticleGenerator::insertCuboid(
      std::array<double, 3>{1.5, 1.5, 0.0}, std::array<size_t, 3>{3, 3,
      1}, 3.0, 1.0, std::array<double, 3>{0.0, 0.0, 0.0}, container);

  EXPECT_TRUE(container.size() == 9);

  Particle p_1 = container[0];

  EXPECT_TRUE(p_1.getX()[0] == 1.5 && p_1.getX()[1] == 1.5 &&
              p_1.getX()[2] == 0.0);

    EXPECT_TRUE(container.get_neighbours(p_1.getType()).size() == 4);

  Particle p_2 = container[4];

  EXPECT_TRUE(p_2.getX()[0] == 4.5 && p_2.getX()[1] == 4.5 &&
              p_2.getX()[2] == 0.0);

    EXPECT_TRUE(container.get_neighbours(p_2.getType()).size() == 9);

  Particle p_3 = container[7];

  EXPECT_TRUE(p_3.getX()[0] == 4.5 && p_3.getX()[1] == 7.5 &&
              p_3.getX()[2] == 0.0);

    EXPECT_TRUE(container.get_neighbours(p_3.getType()).size() == 6);

  // verify every single neigbour for posterity's sake

  auto neighbours = container.get_neighbours(p_3.getType());

  // check all surrounding coordinates
  auto it =
      std::find_if(neighbours.begin(), neighbours.end(), [](ParticlePointer
      p) {
        return p->getX()[0] == 7.5 && p->getX()[1] == 7.5 &&
               p->getX()[2] == 0.0;
      });

  EXPECT_TRUE(it != neighbours.end());
  it =
      std::find_if(neighbours.begin(), neighbours.end(), [](ParticlePointer
      p) {
        return p->getX()[0] == 7.5 && p->getX()[1] == 4.5 &&
               p->getX()[2] == 0.0;
      });
  EXPECT_TRUE(it != neighbours.end());

  it =
      std::find_if(neighbours.begin(), neighbours.end(), [](ParticlePointer
      p) {
        return p->getX()[0] == 4.5 && p->getX()[1] == 4.5 &&
               p->getX()[2] == 0.0;
      });
  EXPECT_TRUE(it != neighbours.end());

  it =
      std::find_if(neighbours.begin(), neighbours.end(), [](ParticlePointer
      p) {
        return p->getX()[0] == 1.5 && p->getX()[1] == 4.5 &&
               p->getX()[2] == 0.0;
      });
  EXPECT_TRUE(it != neighbours.end());

  it = std::find_if(neighbours.begin(), neighbours.end(),
                    [&](ParticlePointer p) {
                      return p->getX()[0] == 1.5 && p->getX()[1] == 7.5 &&
                             p->getX()[2] == 0.0;
                    });
  EXPECT_TRUE(it != neighbours.end());

  ParticleGenerator::insertCuboid(
      std::array<double, 3>{1.5, 1.5, 1.5}, std::array<size_t, 3>{3, 3,
      3}, 3.0, 1.0, std::array<double, 3>{0.0, 0.0, 0.0}, container_3d);

  EXPECT_TRUE(container_3d.size() == 27);

  Particle center_particle = container_3d[13];

  // verify that it is a middle particle
  EXPECT_TRUE(center_particle.getX()[0] == 4.5 &&
              center_particle.getX()[1] == 4.5 &&
              center_particle.getX()[2] == 4.5);

    EXPECT_TRUE(container_3d.get_neighbours(center_particle.getType()).size()
    == 27);

}

TEST_F(LinkedCellTest, UnevenDomainTest) {

  EXPECT_TRUE(uneven_container.cells.size() == 64);

  Particle p(std::array<double, 3>{8.1, 1.0, 1.0},
             std::array<double, 3>{0.0, 0.0, 0.0}, 1.0, 0);
  uneven_container.insert(p, true);

  EXPECT_TRUE(
      uneven_container.is_within_domain(std::array<double,
      3>{8.5, 1.0, 1.0}));

  EXPECT_TRUE(uneven_container.size() == 1);

  EXPECT_EQ(uneven_container.r_cutoff_x, 2.25);
  EXPECT_EQ(uneven_container.r_cutoff_y, 2.25);
  EXPECT_EQ(uneven_container.r_cutoff_z, 2.0);

  // check that the particle is placed in the correct cell. Seems wrong fix
  // later pls
  EXPECT_EQ(uneven_container.cells[3].size(), 1);

  for (size_t i = 0; i < 64; ++i) {
    if (i != 3) {
      EXPECT_EQ(uneven_container.cells[i].size(), 0);
    }
  }

  // insert another particle

  // why does this cause a segfault when set to 0?

  Particle p_2(std::array<double, 3>{1.0, 7.5, 1.0},
               std::array<double, 3>{0.0, 0.0, 0.0}, 1.0, 1);
  uneven_container.insert(p_2, true);

  EXPECT_TRUE(uneven_container.size() == 2);

  EXPECT_EQ(uneven_container.cells[12].size(), 1);

  // move the particle, ensure it is removed from the old cell and inserted
  // into the new cell

  std::array<double, 3> old_position = p_2.getX();

  uneven_container[1].updateX(2.26, 8.3, 1.0);
  uneven_container.update_particle_location(1, old_position);

  EXPECT_TRUE(uneven_container.size() == 2);

  EXPECT_EQ(uneven_container.cells[12].size(), 0);
  EXPECT_EQ(uneven_container.cells[13].size(), 1);
}

TEST_F(LinkedCellTest, RepositioningTest) {

  ParticleGenerator::insertCuboid(
      std::array<double, 3>{5.0, 5.0, 0.0}, std::array<size_t, 3>{2, 2,
      1}, 2.0, 1.0, std::array<double, 3>{0.0, 0.0, 0.0}, container);

  EXPECT_TRUE(container.size() == 4);

  EXPECT_EQ(container.cells[0].size(), 0);
  EXPECT_EQ(container.cells[1].size(), 0);
  EXPECT_EQ(container.cells[2].size(), 0);
  EXPECT_EQ(container.cells[3].size(), 0);
  EXPECT_EQ(container.cells[4].size(), 4);
  EXPECT_EQ(container.cells[5].size(), 0);
  EXPECT_EQ(container.cells[6].size(), 0);
  EXPECT_EQ(container.cells[7].size(), 0);
  EXPECT_EQ(container.cells[8].size(), 0);

  EXPECT_EQ(container.domain_size_.size(), 3);

  EXPECT_EQ(container.left_corner_coordinates[0], 1.5);
  EXPECT_EQ(container.left_corner_coordinates[1], 1.5);
  EXPECT_EQ(container.left_corner_coordinates[2], -1.5);
}

// Test the mark_halo_cells (or mark_boundary_cells) function for a 3D cube. The
// test initializes a LinkedCellContainer with specific parameters and then
// verifies the placement of each boundary cell.
TEST_F(LinkedCellTest, MarkBoundaryCells3D) {

  container_3d.mark_halo_cells();

  size_t x_cells = container_3d.x;
  size_t y_cells = container_3d.y;
  size_t z_cells = container_3d.z;

  for (size_t i = 0; i < x_cells; i++) {
    for (size_t j = 0; j < y_cells; j++) {
      for (size_t k = 0; k < z_cells; k++) {
        size_t index = i + j * x_cells + k * x_cells * y_cells;
        const auto &cell = container_3d.cells[index];

        // Only check halo cells
        if (cell.is_halo) {
          switch (cell.placement) {
          case Placement::BOTTOM_BACK_LEFT_CORNER:
            EXPECT_EQ(i, 0);
            EXPECT_EQ(j, 0);
            EXPECT_EQ(k, 0);
            break;
          case Placement::BOTTOM_FRONT_LEFT_CORNER:
            EXPECT_EQ(i, 0);
            EXPECT_EQ(j, 0);
            EXPECT_EQ(k, z_cells - 1);
            break;
          case Placement::TOP_BACK_LEFT_CORNER:
            EXPECT_EQ(i, 0);
            EXPECT_EQ(j, y_cells - 1);
            EXPECT_EQ(k, 0);
            break;
          case Placement::TOP_FRONT_LEFT_CORNER:
            EXPECT_EQ(i, 0);
            EXPECT_EQ(j, y_cells - 1);
            EXPECT_EQ(k, z_cells - 1);
            break;
          case Placement::BOTTOM_BACK_RIGHT_CORNER:
            EXPECT_EQ(i, x_cells - 1);
            EXPECT_EQ(j, 0);
            EXPECT_EQ(k, 0);
            break;
          case Placement::BOTTOM_FRONT_RIGHT_CORNER:
            EXPECT_EQ(i, x_cells - 1);
            EXPECT_EQ(j, 0);
            EXPECT_EQ(k, z_cells - 1);
            break;
          case Placement::TOP_BACK_RIGHT_CORNER:
            EXPECT_EQ(i, x_cells - 1);
            EXPECT_EQ(j, y_cells - 1);
            EXPECT_EQ(k, 0);
            break;
          case Placement::TOP_FRONT_RIGHT_CORNER:
            EXPECT_EQ(i, x_cells - 1);
            EXPECT_EQ(j, y_cells - 1);
            EXPECT_EQ(k, z_cells - 1);
            break;

          case Placement::LEFT_BACK_EDGE:
            EXPECT_EQ(i, 0);
            EXPECT_EQ(k, 0);
            EXPECT_GT(j, 0);
            EXPECT_LT(j, y_cells - 1);
            break;
          case Placement::LEFT_FRONT_EDGE:
            EXPECT_EQ(i, 0);
            EXPECT_EQ(k, z_cells - 1);
            EXPECT_GT(j, 0);
            EXPECT_LT(j, y_cells - 1);
            break;
          case Placement::RIGHT_BACK_EDGE:
            EXPECT_EQ(i, x_cells - 1);
            EXPECT_EQ(k, 0);
            EXPECT_GT(j, 0);
            EXPECT_LT(j, y_cells - 1);
            break;
          case Placement::RIGHT_FRONT_EDGE:
            EXPECT_EQ(i, x_cells - 1);
            EXPECT_EQ(k, z_cells - 1);
            EXPECT_GT(j, 0);
            EXPECT_LT(j, y_cells - 1);
            break;
          case Placement::TOP_BACK_EDGE:
            EXPECT_EQ(j, y_cells - 1);
            EXPECT_EQ(k, 0);
            EXPECT_GT(i, 0);
            EXPECT_LT(i, x_cells - 1);
            break;
          case Placement::TOP_FRONT_EDGE:
            EXPECT_EQ(j, y_cells - 1);
            EXPECT_EQ(k, z_cells - 1);
            EXPECT_GT(i, 0);
            EXPECT_LT(i, x_cells - 1);
            break;
          case Placement::BOTTOM_BACK_EDGE:
            EXPECT_EQ(j, 0);
            EXPECT_EQ(k, 0);
            EXPECT_GT(i, 0);
            EXPECT_LT(i, x_cells - 1);
            break;
          case Placement::BOTTOM_FRONT_EDGE:
            EXPECT_EQ(j, 0);
            EXPECT_EQ(k, z_cells - 1);
            EXPECT_GT(i, 0);
            EXPECT_LT(i, x_cells - 1);
            break;
          case Placement::BOTTOM_LEFT_EDGE:
            EXPECT_EQ(i, 0);
            EXPECT_EQ(j, 0);
            EXPECT_GT(k, 0);
            EXPECT_LT(k, z_cells - 1);
            break;
          case Placement::TOP_LEFT_EDGE:
            EXPECT_EQ(i, 0);
            EXPECT_EQ(j, y_cells - 1);
            EXPECT_GT(k, 0);
            EXPECT_LT(k, z_cells - 1);
            break;
          case Placement::BOTTOM_RIGHT_EDGE:
            EXPECT_EQ(i, x_cells - 1);
            EXPECT_EQ(j, 0);
            EXPECT_GT(k, 0);
            EXPECT_LT(k, z_cells - 1);
            break;
          case Placement::TOP_RIGHT_EDGE:
            EXPECT_EQ(i, x_cells - 1);
            EXPECT_EQ(j, y_cells - 1);
            EXPECT_GT(k, 0);
            EXPECT_LT(k, z_cells - 1);
            break;

          case Placement::BACK:
            EXPECT_EQ(k, 0);
            EXPECT_GT(i, 0);
            EXPECT_LT(i, x_cells - 1);
            EXPECT_GT(j, 0);
            EXPECT_LT(j, y_cells - 1);
            break;
          case Placement::FRONT:
            EXPECT_EQ(k, z_cells - 1);
            EXPECT_GT(i, 0);
            EXPECT_LT(i, x_cells - 1);
            EXPECT_GT(j, 0);
            EXPECT_LT(j, y_cells - 1);
            break;
          case Placement::BOTTOM:
            EXPECT_EQ(j, 0);
            EXPECT_GT(i, 0);
            EXPECT_LT(i, x_cells - 1);
            EXPECT_GT(k, 0);
            EXPECT_LT(k, z_cells - 1);
            break;
          case Placement::TOP:
            EXPECT_EQ(j, y_cells - 1);
            EXPECT_GT(i, 0);
            EXPECT_LT(i, x_cells - 1);
            EXPECT_GT(k, 0);
            EXPECT_LT(k, z_cells - 1);
            break;
          case Placement::LEFT:
            EXPECT_EQ(i, 0);
            EXPECT_GT(j, 0);
            EXPECT_LT(j, y_cells - 1);
            EXPECT_GT(k, 0);
            EXPECT_LT(k, z_cells - 1);
            break;
          case Placement::RIGHT:
            EXPECT_EQ(i, x_cells - 1);
            EXPECT_GT(j, 0);
            EXPECT_LT(j, y_cells - 1);
            EXPECT_GT(k, 0);
            EXPECT_LT(k, z_cells - 1);
            break;
          default:
            FAIL() << "Invalid placement for boundary cell at index " << index;
          }
        }
      }
    }
  }
}
