#include "simulator/particle/Particle.h"
#include "simulator/particle/container/ParticleContainer.h"
#include "utils/ArrayUtils.h"
#include "utils/logger/Logger.h"
#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include <spdlog/spdlog.h>
#include <vector>
#include <algorithm>
#include <utility>

// The fixture for testing class Foo.
class BasicTest : public testing::Test {
protected:
  BasicTest() : container{}, starting_coord{1.0}, delta_t{1.0} {}

  ParticleContainer container;
  double starting_coord;
  double delta_t;
  bool calculateLJForce = true;
  Logger &logger = Logger::getInstance("debug");

  void calculateF() {
    // store the current force as the old force and reset current to 0
    for (auto &p : container) {
      p.updateOldF(p.getF());
      p.updateF(0, 0, 0);
    }

    // Iterate each pair
    for (auto it = container.pair_begin(); it != container.pair_end(); ++it) {
      ParticlePair &pair = *it;
      Particle &p1 = *(pair.first);
      Particle &p2 = *(pair.second);
      auto r12 = p2.getX() - p1.getX();
      // distance ||x_i - x_j ||
      double distance = ArrayUtils::L2Norm(r12);

      // avoid extermely small distance
      if (distance > 1e-5) {
        // switch Lennard-Jones/ Simple force
        double totalForce;
        if (calculateLJForce) {
          // Lennard-Jones parameters
          const double epsilon = 5.0;
          const double sigma = 1.0;
          // Lennard-Jones Force Formula (3)
          double term = sigma / distance;
          double term6 = pow(term, 6);
          double term12 = pow(term, 12);
          totalForce = 24 * epsilon * (term6 - 2 * term12) / distance;
        } else {
          // Simple Force Calculation Formula (14)
          totalForce = p1.getM() * p2.getM() / pow(distance, 2);
        }
        auto force = (totalForce / distance) * r12;

        p1.updateF(p1.getF() + force);
        // Newton's third law
        p2.updateF(p2.getF() - force);
      }
    }
  }

  void calculateX() {
    for (auto &p : container) {
      auto x = p.getX();
      auto v = p.getV();
      auto f = p.getF();
      double m = p.getM();

      // Velocity-Störmer-Verlet formula (8)
      x = x + delta_t * v + pow(delta_t, 2) * f / (2 * m);

      p.updateX(x);
    }
  }

  void calculateV() {
    for (auto &p : container) {
      auto v = p.getV();
      auto old_f = p.getOldF();
      auto new_f = p.getF();
      double m = p.getM();

      // Velocity-Störmer-Verlet formula (9)
      v = v + delta_t * (old_f + new_f) / (2 * m);

      p.updateV(v);
    }
  }
};

// create 3 arbitrary particles, ensure that pairs are built correctly and no
// excessive objects are created
TEST_F(BasicTest, ContainerBehaviourTest) {

  Particle p_1{std::array{1.0, 0.0, 0.0}, std::array{1.0, 0.0, 0.0}, 1.0, 0};
  Particle p_2{std::array{2.0, 0.0, 0.0}, std::array{1.0, 0.0, 0.0}, 1.0, 1};
  Particle p_3{std::array{3.0, 0.0, 0.0}, std::array{1.0, 0.0, 0.0}, 1.0, 2};

  container.insert(p_1);
  container.insert(p_2);
  container.insert(p_3);

  int id = 0;
  for (auto const &p : container) {
    Particle dummy_particle{std::array{starting_coord, 0.0, 0.0},
                            std::array{1.0, 0.0, 0.0}, 1.0, id};
    EXPECT_EQ(p, dummy_particle);
    id++;
    starting_coord++;
  }

  
  // check that the correct amount of pairs have been created

  size_t counter = 0;

  for (auto it = container.pair_begin(); it != container.pair_end(); ++it)
    counter++;

  ASSERT_EQ(counter, 3);

  // check that exactly the right particles are contained within the set

  std::vector<ParticlePair> dummy_pairs{
      ParticlePair(std::make_shared<Particle>(p_1),
                   std::make_shared<Particle>(p_2)),
      ParticlePair(std::make_shared<Particle>(p_1),
                   std::make_shared<Particle>(p_3)),
      ParticlePair(std::make_shared<Particle>(p_2),
                   std::make_shared<Particle>(p_3))};

  for (auto it = container.pair_begin(); it != container.pair_end(); ++it){
    EXPECT_TRUE(std::find(dummy_pairs.begin(), dummy_pairs.end(), *it) != dummy_pairs.end());
  }

  // check that every particle maps onto correct pair
}


//Perform a basic test to check particle movement
TEST_F(BasicTest, SimulationBehaviourTest) {
  // this in place initialization looks super ugly, perhaps we can use
  // std::initializer_list
  Particle p_1{std::array{0.0, 0.0, 0.0}, std::array{1.0, 0.0, 0.0}, 1.0, 0};

  container.insert(p_1);

  // Perform 3 Iterations

  for (size_t i{0}; i < 3; i++) {
    calculateX();
    calculateF();
    calculateV();
    ASSERT_TRUE((container[0].getX() ==
                 std::array<double, 3>{static_cast<double>(i + 1), 0.0, 0.0}));
    ASSERT_TRUE((container[0].getF() == std::array<double, 3>{0.0, 0.0, 0.0}));
    ASSERT_TRUE((container[0].getV() == std::array<double, 3>{1.0, 0.0, 0.0}));
  }
}

//Perform a basic test to check force calculation
TEST_F(BasicTest, CalculateFTest) {
  calculateLJForce = true;
  delta_t = 0.02;

  // Initialize two particles
  Particle p_1{std::array{0.0, 0.0, 0.0}, std::array{0.0, 0.0, 0.0}, 1.0, 0};
  Particle p_2{std::array{1.2, 0.0, 0.0}, std::array{0.0, 0.0, 0.0}, 1.0, 1};

  container.insert(p_1);
  container.insert(p_2);

  // iteration 1
  calculateF();
  // calculate by hand
  auto expectedForce = std::array{11.058, 0.00, 0.00};
  // Verify the forces
  ASSERT_NEAR(container[0].getF()[0], expectedForce[0], 1e-2);
  ASSERT_NEAR(container[0].getF()[1], expectedForce[1], 1e-2);
  ASSERT_NEAR(container[0].getF()[2], expectedForce[2], 1e-2);

  ASSERT_NEAR(container[1].getF()[0], -expectedForce[0], 1e-2);
  ASSERT_NEAR(container[1].getF()[1], -expectedForce[1], 1e-2);
  ASSERT_NEAR(container[1].getF()[2], -expectedForce[2], 1e-2);

  // Update positions for the next iteration
  calculateX();
  logger.debug("current p1 position: " +
               ArrayUtils::to_string(container[0].getX()));
  logger.debug("current p2 position: " +
               ArrayUtils::to_string(container[1].getX()));

  // iteration 2
  calculateF();
  // calculate by hand
  expectedForce = std::array{10.832, 0.00, 0.00};
  // Verify the forces
  ASSERT_NEAR(container[0].getF()[0], expectedForce[0], 1e-2);
  ASSERT_NEAR(container[0].getF()[1], expectedForce[1], 1e-2);
  ASSERT_NEAR(container[0].getF()[2], expectedForce[2], 1e-2);

  ASSERT_NEAR(container[1].getF()[0], -expectedForce[0], 1e-2);
  ASSERT_NEAR(container[1].getF()[1], -expectedForce[1], 1e-2);
  ASSERT_NEAR(container[1].getF()[2], -expectedForce[2], 1e-2);

  delta_t = 0.2;

  calculateX();
  logger.debug("current p1 position: " +
               ArrayUtils::to_string(container[0].getX()));
  logger.debug("current p2 position: " +
               ArrayUtils::to_string(container[1].getX()));

  // iteration 3
  calculateF();
  // calculate by hand
  expectedForce = std::array{-7376.511, 0.00, 0.00};
  // Verify the forces
  ASSERT_NEAR(container[0].getF()[0], expectedForce[0], 1);
  ASSERT_NEAR(container[0].getF()[1], expectedForce[1], 1);
  ASSERT_NEAR(container[0].getF()[2], expectedForce[2], 1);

  ASSERT_NEAR(container[1].getF()[0], -expectedForce[0], 1);
  ASSERT_NEAR(container[1].getF()[1], -expectedForce[1], 1);
  ASSERT_NEAR(container[1].getF()[2], -expectedForce[2], 1);
}
