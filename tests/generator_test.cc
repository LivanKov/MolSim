#include "logger/Logger.h"
#include "particleSim/ParticleContainer.h"
#include "particleSim/ParticleGenerator.h"
#include <gtest/gtest.h>
#include <iostream>
#include <numeric>

class CuboidTest : public testing::Test {
protected:
  CuboidTest() : container{} {}
  ParticleContainer container;
};

TEST_F(CuboidTest, SimpleTest) {
  ParticleGenerator::insertCuboid(
      std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{3, 3, 3}, 1.0, 1.0,
      std::array<double, 3>{1.0, 1.0, 1.0}, 0, container);

  size_t index = 0;

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      for (size_t k = 0; k < 3; k++) {
        Particle dummyParticle{std::array<double, 3>{static_cast<double>(k),
                                                     static_cast<double>(j),
                                                     static_cast<double>(i)},
                               std::array<double, 3>{1.0, 1.0, 1.0}, 1.0};
        ASSERT_EQ(container[index++], dummyParticle);
      }
    }
  }
}

TEST_F(CuboidTest, MultipleCuboidsTest) {
  ParticleGenerator::insertCuboid(
      std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{3, 3, 3}, 1.0, 1.0,
      std::array<double, 3>{1.0, 1.0, 1.0}, 0, container);
  ParticleGenerator::insertCuboid(
      std::array<double, 3>{10, 10, 0}, std::array<size_t, 3>{3, 5, 4}, 1.0,
      1.0, std::array<double, 3>{1.0, 1.0, 1.0}, 0, container);

  size_t index = 0;

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      for (size_t k = 0; k < 3; k++) {
        Particle dummyParticle{std::array<double, 3>{static_cast<double>(k),
                                                     static_cast<double>(j),
                                                     static_cast<double>(i)},
                               std::array<double, 3>{1.0, 1.0, 1.0}, 1.0};
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
                               std::array<double, 3>{1.0, 1.0, 1.0}, 1.0};
        ASSERT_EQ(container[index++], dummyParticle);
      }
    }
  }
}
