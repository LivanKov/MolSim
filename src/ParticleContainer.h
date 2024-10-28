/*
 * ParticleContainer.h
 */

#pragma once

#include "Particle.h"
#include <vector>
#include <stdexcept>
#include <array>

struct ParticlePair{
  std::array<double, 3> f;
  std::array<double, 3> old_f;
  
  Particle first;
  Particle second;

  ParticlePair(const Particle& first, const Particle& second);

  bool operator==(const Particle& rhs) const;
};


class ParticleContainer {

public:
  ~ParticleContainer() = default;
  ParticleContainer(const ParticleContainer &lhs) = delete;
  ParticleContainer &operator=(const ParticleContainer &lhs) = delete;
  ParticleContainer(ParticleContainer &&lhs) = delete;
  ParticleContainer &operator=(ParticleContainer &&lhs) = delete;

  template <typename... Args> void emplace_back(Args... args);

  void insert(Particle &p);

  void insert(Particle &&p);

  size_t size();

private:
  std::vector<Particle> particle_container;
  
};