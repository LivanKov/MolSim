#include <gtest/gtest.h>
#include "ParticleContainer.h"
#include "Particle.h"
#include <array>
#include <iostream>
#include <unordered_set>
#include <utility>
#include <memory>


//create 3 arbitrary particles, ensure that pairs are built correctly and no excessive objects are created
TEST(HelloTest, ParticleContainerInitializationTest) {
  ParticleContainer container{};
  //this in place initialization looks super ugly, perhaps we can use std::initializer_list

  Particle p_1{std::array{1.0, 0.0,0.0}, std::array{1.0, 0.0, 0.0}, 1.0};
  Particle p_2{std::array{2.0, 0.0,0.0}, std::array{1.0, 0.0, 0.0}, 1.0};
  Particle p_3{std::array{3.0, 0.0, 0.0}, std::array{1.0, 0.0, 0.0}, 1.0};

  container.insert(std::forward<Particle>(p_1));
  container.insert(std::forward<Particle>(p_2));
  container.insert(std::forward<Particle>(p_3));
  
  //check if the correct amount of particles have been created within the container
  ASSERT_EQ(container.size(), 3);

  // check each particle individually
  double starting_coord = 1.0;
  for(auto const& p : container){
    Particle dummy_particle{std::array{starting_coord, 0.0,0.0}, std::array{1.0, 0.0, 0.0}, 1.0};
    EXPECT_EQ(p, dummy_particle);
    starting_coord++;
  }

  //check that the correct amount of pairs have been created

  size_t counter = 0;

  for(auto it = container.pair_begin(); it != container.pair_end(); ++it)
    counter++;

  ASSERT_EQ(counter, 3);

  //check that exactly the right particles are contained within the set
  std::unordered_set<ParticlePair> dummy_pairs{
    ParticlePair(std::make_shared<Particle>(p_1),std::make_shared<Particle>(p_2)), 
    ParticlePair(std::make_shared<Particle>(p_1),std::make_shared<Particle>(p_3)), 
    ParticlePair(std::make_shared<Particle>(p_2),std::make_shared<Particle>(p_3))
  };

  for(auto it = container.pair_begin(); it != container.pair_end(); ++it)
    EXPECT_TRUE(dummy_pairs.contains(*it));


  //check that every particle maps onto correct pair

  for(auto const& p : container){
    ASSERT_EQ(container[p].size(),2);
    for(auto const& v : container[p]){
      ASSERT_TRUE(*(v->first) == p && *(v->second) != p || *(v->first) != p && *(v->second) == p);
    }
  }
}