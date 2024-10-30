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

  bool operator==(const ParticlePair &rhs) const;
};
template <>
struct std::hash<ParticlePair>
{
  std::size_t operator()(const ParticlePair &p) const;
};

class ParticleIterator
{
public:
  using PValueType = Particle;
  using PPointerType = Particle *;
  using PReferenceType = Particle &;
  ParticleIterator(PPointerType p);
  ParticleIterator &operator++();
  ParticleIterator operator++(int);
  PPointerType operator->();
  bool operator==(ParticleIterator &rhs) const;
  bool operator!=(ParticleIterator &rhs) const;
  PReferenceType operator*() const;

private:
  PPointerType _ptr;
};

class ParticlePairIterator
{
public:
  using PValueType = const Particle;
  using PPointerType = std::unordered_set<ParticlePair>::const_iterator;
  using PReferenceType = const ParticlePair &;
  ParticlePairIterator(PPointerType p);
  ParticlePairIterator &operator++();
  ParticlePairIterator operator++(int);
  PPointerType operator->();
  bool operator==(ParticlePairIterator &rhs) const;
  bool operator!=(ParticlePairIterator &rhs) const;
  PReferenceType operator*() const;

private:
  PPointerType _ptr;
};

class ParticleContainer
{

public:
  ParticleContainer();
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
    create_pairs(Particle{std::forward<Args>(args)...});
  }

  void insert(Particle &p);

  void insert(Particle &&p);

  size_t size();

  std::vector<ParticlePair> &pairs_of(const Particle &p);

  std::vector<ParticlePair> &operator[](const Particle &p);

  Particle &operator[](int index);

  ParticleIterator begin();
  ParticleIterator end();
  ParticlePairIterator pair_begin();
  ParticlePairIterator pair_end();

private:
  void create_pairs(const Particle &p);

  std::vector<Particle> _particle_container;
  std::unordered_set<ParticlePair> _particle_pair_set;
  std::unordered_map<Particle, std::vector<ParticlePair>> _particle_pair_map;
};