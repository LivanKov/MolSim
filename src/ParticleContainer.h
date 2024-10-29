/*
 * ParticleContainer.h
 */

#pragma once

#include "Particle.h"
#include <vector>
#include <stdexcept>
#include <array>
#include <unordered_set>
#include <unordered_map>

struct ParticlePair
{
  std::array<double, 3> f;
  std::array<double, 3> old_f;

  Particle first;
  Particle second;

  ParticlePair(const Particle &first, const Particle &second);

  bool operator==(ParticlePair &rhs);
};

class ParticleContainer
{

public:
  ParticleContainer() = default;
  ~ParticleContainer() = default;
  ParticleContainer(const ParticleContainer &lhs) = delete;
  ParticleContainer &operator=(const ParticleContainer &lhs) = delete;
  ParticleContainer(ParticleContainer &&lhs) = delete;
  ParticleContainer &operator=(ParticleContainer &&lhs) = delete;

  template <typename... Args>
    requires std::constructible_from<Particle, Args...>
  void emplace_back(Args... args)
  {
    _particle_container.emplace_back(std::forward<Args>(args)...);
  }

  void insert(Particle &p);

  void insert(Particle &&p);

  size_t size();

  void pairs_of(Particle &p);

  class PContainerIterator
  {
  public:
    using PValueType = Particle;
    using PPointerType = Particle *;
    using PReferenceType = Particle &;
    PContainerIterator(PPointerType p);
    PContainerIterator &operator++();
    PContainerIterator operator++(int);
    PPointerType operator->();
    bool operator==(PContainerIterator &rhs) const;
    bool operator!=(PContainerIterator &rhs) const;
    PReferenceType operator*() const;

  private:
    PPointerType _ptr;
  };

  PContainerIterator begin();
  PContainerIterator end();

private:
  std::vector<Particle> _particle_container;
  std::unordered_map<Particle, std::vector<ParticlePair>> _particle_pair_map;
};