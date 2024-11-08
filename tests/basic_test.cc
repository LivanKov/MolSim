#include "Particle.h"
#include "ParticleContainer.h"
#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <utility>

// The fixture for testing class Foo.
class BasicTest : public testing::Test {
protected:
  BasicTest() : container{}, starting_coord{1.0}, delta_t{1.0} {}

  ParticleContainer container;
  double starting_coord;
  double delta_t;

  void calculateF() {
    for (auto &p1 : container) {
      std::array<double, 3> force_copy = p1.getF();
      for (auto &p2 : container) {
        double f_x, f_y, f_z = 0;
        double distance = std::sqrt(std::pow(p1.getX()[0] - p2.getX()[0], 2) +
                                    std::pow(p1.getX()[1] - p2.getX()[1], 2) +
                                    std::pow(p1.getX()[2] - p2.getX()[2], 2));
        if (!(p1 == p2)) {
          f_x = (p2.getX()[0] - p1.getX()[0]) * (p1.getM() * p2.getM()) /
                pow(distance, 3);
          f_y = (p2.getX()[1] - p1.getX()[1]) * (p1.getM() * p2.getM()) /
                pow(distance, 3);
          f_z = (p2.getX()[2] - p1.getX()[2]) * (p1.getM() * p2.getM()) /
                pow(distance, 3);
          p1.updateF(p1.getF()[0] + f_x, p1.getF()[1] + f_y,
                     p1.getF()[2] + f_z);
        }
      }
      p1.updateOldF(force_copy[0], force_copy[1], force_copy[2]);
    }
  }

  void calculateX() {
    for (auto &p : container) {
      auto x = p.getX();
      auto v = p.getV();
      auto f = p.getF();
      double m = p.getM();

      x[0] = x[0] + delta_t * v[0] + pow(delta_t, 2) * f[0] / (2 * m);
      x[1] = x[1] + delta_t * v[1] + pow(delta_t, 2) * f[1] / (2 * m);
      x[2] = x[2] + delta_t * v[2] + pow(delta_t, 2) * f[2] / (2 * m);

      p.updateX(x[0], x[1], x[2]);
    }
  }

  void calculateV() {
    for (auto &p : container) {
      auto v = p.getV();
      auto old_f = p.getOldF();
      auto new_f = p.getF();
      double m = p.getM();

      // Velocity-Störmer-Verlet formula (9)
      v[0] = v[0] + delta_t * (old_f[0] + new_f[0]) / (2 * m);
      v[1] = v[1] + delta_t * (old_f[1] + new_f[1]) / (2 * m);
      v[2] = v[2] + delta_t * (old_f[2] + new_f[2]) / (2 * m);

      p.updateV(v[0], v[1], v[2]);
    }
  }
};

// create 3 arbitrary particles, ensure that pairs are built correctly and no
// excessive objects are created
TEST_F(BasicTest, ContainerBehaviourTest) {

  Particle p_1{std::array{1.0, 0.0, 0.0}, std::array{1.0, 0.0, 0.0}, 1.0};
  Particle p_2{std::array{2.0, 0.0, 0.0}, std::array{1.0, 0.0, 0.0}, 1.0};
  Particle p_3{std::array{3.0, 0.0, 0.0}, std::array{1.0, 0.0, 0.0}, 1.0};

  container.insert(p_1);
  container.insert(p_2);
  container.insert(p_3);

  for (auto const &p : container) {
    Particle dummy_particle{std::array{starting_coord, 0.0, 0.0},
                            std::array{1.0, 0.0, 0.0}, 1.0};
    EXPECT_EQ(p, dummy_particle);
    starting_coord++;
  }

  // check that the correct amount of pairs have been created

  size_t counter = 0;

  for (auto it = container.pair_begin(); it != container.pair_end(); ++it)
    counter++;

  ASSERT_EQ(counter, 3);

  // check that exactly the right particles are contained within the set

  std::unordered_set<ParticlePair> dummy_pairs{
      ParticlePair(std::make_shared<Particle>(p_1),
                   std::make_shared<Particle>(p_2)),
      ParticlePair(std::make_shared<Particle>(p_1),
                   std::make_shared<Particle>(p_3)),
      ParticlePair(std::make_shared<Particle>(p_2),
                   std::make_shared<Particle>(p_3))};

  for (auto it = container.pair_begin(); it != container.pair_end(); ++it)
    EXPECT_TRUE(dummy_pairs.contains(*it));

  // check that every particle maps onto correct pair

  for (auto const &p : container) {
    ASSERT_EQ(container[p].size(), 2);
    for (auto const &v : container[p]) {
      ASSERT_TRUE(*(v->first) == p && *(v->second) != p ||
                  *(v->first) != p && *(v->second) == p);
    }
  }
}

TEST_F(BasicTest, SimulationBehaviourTest) {
  // this in place initialization looks super ugly, perhaps we can use
  // std::initializer_list
  Particle p_1{std::array{0.0, 0.0, 0.0}, std::array{1.0, 0.0, 0.0}, 1.0};

  container.insert(std::forward<Particle>(p_1));

  // Perform 3 Iterations

  for (size_t i{0}; i < 3; i++) {
    calculateX();
    calculateF();
    calculateV();
    ASSERT_TRUE((container[0].getX() == std::array<double, 3>{static_cast<double>(i+1), 0.0, 0.0}));
    ASSERT_TRUE((container[0].getF() == std::array<double, 3>{0.0, 0.0, 0.0}));
    ASSERT_TRUE((container[0].getV() == std::array<double, 3>{1.0, 0.0, 0.0}));
  }
}
