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
#include <memory>

using ParticlePointer = std::shared_ptr<Particle>;

struct ParticlePair
{
  std::array<double, 3> f;
  std::array<double, 3> old_f;

  ParticlePointer first;
  ParticlePointer second;

  ParticlePair(const ParticlePointer first, const ParticlePointer second);

  bool operator==(const ParticlePair &rhs) const;

  std::string toString() const;
};
template <>
struct std::hash<ParticlePair>
{
  std::size_t operator()(const ParticlePair &p) const;
};

struct ParticlePointerHash
{
  size_t operator()(const std::shared_ptr<Particle> &s) const noexcept;
};

struct ParticlePointerEqual
{
  bool operator()(const std::shared_ptr<Particle> &a,
                  const std::shared_ptr<Particle> &b) const;
};

using ParticlePairPointer = std::shared_ptr<ParticlePair>;

class ParticleIterator
{
public:
  using PValueType = Particle;
  using PPointerType = ParticlePointer *;
  using PReferenceType = Particle &;
  ParticleIterator(PPointerType p);
  ParticleIterator &operator++();
  ParticleIterator operator++(int);
  PPointerType operator->();
  bool operator==(const ParticleIterator &rhs) const;
  bool operator!=(const ParticleIterator &rhs) const;
  PReferenceType operator*() const;

private:
  PPointerType _ptr;
};

class ParticlePairIterator
{
public:
  using PValueType = ParticlePair;
  using PPointerType = std::unordered_set<ParticlePairPointer>::const_iterator;
  using PReferenceType = ParticlePair &;
  ParticlePairIterator(PPointerType p);
  ParticlePairIterator &operator++();
  ParticlePairIterator operator++(int);
  PPointerType operator->();
  bool operator==(const ParticlePairIterator &rhs) const;
  bool operator!=(const ParticlePairIterator &rhs) const;
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
    ParticlePointer p = std::make_shared<Particle>(std::forward<Args>(args)...);
    _particle_container.push_back(p);
    create_pairs(p);
  }

  void insert(Particle &p);

  void insert(Particle &&p);

  size_t size();

  std::vector<ParticlePairPointer> &pairs_of(const ParticlePointer &p);

  std::vector<ParticlePairPointer> &operator[](const Particle &p);

  Particle &operator[](int index);

  ParticleIterator begin();
  ParticleIterator end();
  ParticlePairIterator pair_begin();
  ParticlePairIterator pair_end();

private:
  void create_pairs(const ParticlePointer &new_particle);

  std::vector<ParticlePointer> _particle_container;
  std::unordered_set<ParticlePairPointer> _particle_pair_set;
  std::unordered_map<ParticlePointer, std::vector<ParticlePairPointer>,ParticlePointerHash,ParticlePointerEqual> _particle_pair_map;
};