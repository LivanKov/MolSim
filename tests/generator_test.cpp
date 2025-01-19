#include "simulator/particle/container/LinkedCellContainer.h"
#include "simulator/particle/ParticleGenerator.h"
#include "utils/logger/Logger.h"
#include <gtest/gtest.h>
#include <iostream>
#include <numeric>

class CuboidTest : public testing::Test {
protected:
  CuboidTest() : container{} {}
  LinkedCellContainer container;
};

// simple test for the cuboid
// one simple cuboid with leftBottomCorner [0,0,0] and 3x3x3 format
TEST_F(CuboidTest, SimpleTest) {
  ParticleGenerator::insertCuboid(
      std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{3, 3, 3}, 1.0, 1.0,
      std::array<double, 3>{1.0, 1.0, 1.0}, container);

  size_t index = 0;
  int id = 0;

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      for (size_t k = 0; k < 3; k++) {
        Particle dummyParticle{std::array<double, 3>{static_cast<double>(k),
                                                     static_cast<double>(j),
                                                     static_cast<double>(i)},
                               std::array<double, 3>{1.0, 1.0, 1.0}, 1.0, id};
        id++;
        ASSERT_EQ(container[index++], dummyParticle);
      }
    }
  }
}

// test for multiple cuboids. both of them 3x3x3 format
TEST_F(CuboidTest, MultipleCuboidsTest) {
  ParticleGenerator::insertCuboid(
      std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{3, 3, 3}, 1.0, 1.0,
      std::array<double, 3>{1.0, 1.0, 1.0}, container);
  ParticleGenerator::insertCuboid(
      std::array<double, 3>{10, 10, 0}, std::array<size_t, 3>{3, 5, 4}, 1.0,
      1.0, std::array<double, 3>{1.0, 1.0, 1.0}, container);

  size_t index = 0;
  int id = 0;

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      for (size_t k = 0; k < 3; k++) {
        Particle dummyParticle{std::array<double, 3>{static_cast<double>(k),
                                                     static_cast<double>(j),
                                                     static_cast<double>(i)},
                               std::array<double, 3>{1.0, 1.0, 1.0}, 1.0, id};
        id++;
        ASSERT_EQ(container[index++], dummyParticle);
      }
    }
  }

  for (size_t i = 0; i < 4; i++) {
    for (size_t j = 10; j < 5; j++) {
      for (size_t k = 10; k < 3; k++) {
        Particle dummyParticle{std::array<double, 3>{static_cast<double>(k),
                                                     static_cast<double>(j),
                                                     static_cast<double>(i)},
                               std::array<double, 3>{1.0, 1.0, 1.0}, 1.0, id};
        id++;
        ASSERT_EQ(container[index++], dummyParticle);
      }
    }
  }
}

// ---- Disc Tests ----

class DiscTest : public testing::Test {
protected:
  DiscTest() : container{} {}
  LinkedCellContainer container;
};

// Disc test case: checks if correct amount of particles got inserted
// simple use case with radius 1
TEST_F(DiscTest, SizeSimpleDiscTest) {
  const std::array<double, 3> &center{0,0,0};
  const std::array<double, 3> &initialVelocity{0,0,0};
  size_t radius{1};
  double h{2^(1/6)};
  double mass{1.0};

  ParticleGenerator::insertDisc(center,initialVelocity,radius,h,mass, container);

  ASSERT_EQ(container.size(), 5);
}

// Disc test case: checks if correct amount of particles got inserted
// edge case: radius 0 used
TEST_F(DiscTest, SizeSimpleEdgeDiscTest) {
  const std::array<double, 3> &center{0,0,0};
  const std::array<double, 3> &initialVelocity{0,0,0};
  size_t radius{0};
  double h{2^(1/6)};
  double mass{1.0};

  ParticleGenerator::insertDisc(center,initialVelocity,radius,h,mass,container);

  ASSERT_EQ(container.size(), 1);
}

// Disc test case: checks if correct amount of particles got inserted
// use case with radius 3 -> size should be 29 using the condition realRadius = radius*h
TEST_F(DiscTest, SizeDiscTest) {
  const std::array<double, 3> &center{0,0,0};
  const std::array<double, 3> &initialVelocity{0,0,0};
  size_t radius{3};
  double h{2^(1/6)};
  double mass{1.0};

  ParticleGenerator::insertDisc(center,initialVelocity,radius,h,mass,container);

  ASSERT_EQ(container.size(), 29);
}

// checks if the right center particle got inserted
TEST_F(DiscTest, CenterSimpleDiscTest) {
  const std::array<double, 3> &center{0,0,0};
  const std::array<double, 3> &initialVelocity{0,0,0};
  size_t radius{0};
  double h{2^(1/6)};
  double mass{1.0};

  ParticleGenerator::insertDisc(center,initialVelocity,radius,h,mass,container);

  Particle dummyParticle{std::array<double, 3>{0,0,0},
                               std::array<double, 3>{0, 0, 0}, 1.0, 0};

  ASSERT_EQ(container[0], dummyParticle);
}

// checks if the right center particle got inserted using an arbitrary particle
TEST_F(DiscTest, ArbitraryCenterDiscTest) {
  const std::array<double, 3> &center{67.4,-98.2,42.0};
  const std::array<double, 3> &initialVelocity{12.4,55.9,-76.1};
  size_t radius{0};
  double h{2^(1/6)};
  double mass{5.8};

  ParticleGenerator::insertDisc(center,initialVelocity,radius,h,mass,container);

  Particle dummyParticle{std::array<double, 3>{67.4,-98.2,42.0},
                               std::array<double, 3>{12.4,55.9,-76.1}, 5.8, 0};

  ASSERT_EQ(container[0], dummyParticle);
}

// checks if other particles are corect, also checks if the correct h is used
TEST_F(DiscTest, OtherParticlesDiscTest) {
  const std::array<double, 3> &center{0,0,0};
  const std::array<double, 3> &initialVelocity{1,1,1};
  size_t radius{1};
  double h{2^(1/6)};
  double mass{1.0};

  ParticleGenerator::insertDisc(center,initialVelocity,radius,h,mass,container);

  // particle left to center
  Particle dummyParticle1{std::array<double, 3>{-h,0,0},
                               std::array<double, 3>{1,1,1}, 1.0, 0};
  ASSERT_EQ(container[0], dummyParticle1);

  // particle under center
  Particle dummyParticle2{std::array<double, 3>{0,-h,0},
                               std::array<double, 3>{1,1,1}, 1.0, 1};
  ASSERT_EQ(container[1], dummyParticle2);

  // center particle
  Particle dummyParticle3{std::array<double, 3>{0,0,0},
                               std::array<double, 3>{1,1,1}, 1.0, 2};
  ASSERT_EQ(container[2], dummyParticle3);

  // particle above center
  Particle dummyParticle4{std::array<double, 3>{0,h,0},
                               std::array<double, 3>{1,1,1}, 1.0, 3};
  ASSERT_EQ(container[3], dummyParticle4);

  // particle right to center
  Particle dummyParticle5{std::array<double, 3>{h,0,0},
                               std::array<double, 3>{1,1,1}, 1.0, 4};
  ASSERT_EQ(container[4], dummyParticle5);
}
