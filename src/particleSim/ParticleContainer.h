/*
 * ParticleContainer.h
 */

#pragma once

#include "Particle.h"
#include <array>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using ParticlePointer = std::shared_ptr<Particle>;

/**
 * @struct ParticlePair.
 * @brief Manages two shared pointers to unique Particle objects. Stores the
 * force between the two particles.
 */
struct ParticlePair {
  std::array<double, 3> f;
  std::array<double, 3> old_f;

  ParticlePointer first;
  ParticlePointer second;

  /**
   * @brief Constructor.
   * @param first Pointer managing an instance of Particle class.
   * @param secondd Pointer managing an instance of Particle class.
   */
  ParticlePair(const ParticlePointer first, const ParticlePointer second);

  /**
   * @brief Override equality operator for ParticlePair class.
   * @param rhs ParticlePair object reference.
   * @return Boolean.
   */
  bool operator==(const ParticlePair &rhs) const;

  /**
   * @brief Output string representation of ParticlePair object.
   * @return Std::string representation.
   */
  std::string toString() const;
};
/**
 * @brief Custom hashing function, disregards the order of the objects stored in
 * the pair
 * @param p Reference to ParticlePair object.
 * @return std::size_t hash.
 */
template <> struct std::hash<ParticlePair> {
  std::size_t operator()(const ParticlePair &p) const;
};
/**
 * @brief Custom hashing function object, disregards the order of the objects
 * stored in the pair.
 * @param s Reference to ParticlePair object.
 * @return Std::size_t hash.
 */
struct ParticlePointerHash {
  size_t operator()(const std::shared_ptr<Particle> &s) const noexcept;
};
/**
 * @brief Struct to determine equality of two Particle objects managed by shared
 * pointers. Used for hashing in a std::unordered_map.
 * @param a Reference to a pointer managing a Particle object.
 * @param b Reference to a pointer managing a Particle object.
 * @return Boolean.
 */
struct ParticlePointerEqual {
  bool operator()(const std::shared_ptr<Particle> &a,
                  const std::shared_ptr<Particle> &b) const;
};

using ParticlePairPointer = std::shared_ptr<ParticlePair>;
/**
 * @class ParticleIterator
 * @brief Provides an iterator interface for the ParticleContainer class.
 * Based on the now deprecated std::iterator. See
 * https://en.cppreference.com/w/cpp/iterator/iterator for more details.
 */
class ParticleIterator {
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

/**
 * @class ParticlePairterator
 * @brief Provides an iterator interface that allows to iterate over unique
 * pairs stored in container. Can speed up computation by avoiding excessive
 * operations. Based on the now deprecated std::iterator. See
 * https://en.cppreference.com/w/cpp/iterator/iterator for more details.
 */

class ParticlePairIterator {
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

/**
 * @class ParticleContainer
 * @brief Main container class, contains an intuitive api for storing particles
 * and running unit tests. Stores unique particles and corresponding unique
 * pairs, managed by shared pointers.
 */

class ParticleContainer {

public:
  /**
   * @brief Constructor.
   */
  ParticleContainer();

  /**
   * @brief Create a Particle object in-place and insert it into the underlying
   * containers, allocates place for more pairs.
   * @param args Particle object constructor arguments.
   */

  template <typename... Args>
  requires std::constructible_from<Particle, Args...> void
  emplace_back(Args... args) {
    ParticlePointer p = std::make_shared<Particle>(std::forward<Args>(args)...);
    _particle_container.push_back(p);
    create_pairs(p);
  }
  /**
   * @brief Inserts a Particle object into the container. Redirects the call to
   * emplace.
   * @param p Particle object reference.
   */

  void insert(Particle &p);

  /**
   * @brief Inserts a Particle object into the container. Redirects the call to
   * emplace.
   * @param p Particle object rvalue reference.
   */

  void insert(Particle &&p);

  /**
   * @brief Returns the amount of unique particles stored in the container.
   * @return Size_t container size.
   */
  size_t size();

  /**
   * @brief Returns a Particle object stored at index.
   * @param index Index of the Particle in the underlying container.
   * @return Reference to a Particle object.
   */
  Particle &operator[](size_t index);

  /**
   * @brief Clears the container.
   */
  void clear();

  /**
   * @brief Iterator interface for the main iterator.
   */
  ParticleIterator begin();
  /**
   * @brief Iterator interface for the main iterator.
   */
  ParticleIterator end();
  /**
   * @brief Iterator interface for the pair iterator.
   */
  ParticlePairIterator pair_begin();
  /**
   * @brief Iterator interface for the pair iterator.
   */
  ParticlePairIterator pair_end();

private:
  /**
   * @brief Creates particle pairs for newly inserted particle.
   * @param new_particle Reference to a pointer managing the newly inserted
   * Particle.
   */
  void create_pairs(const ParticlePointer &new_particle);

  /**
   * @brief Main underlying container. Manages pointers to Particle objects.
   */
  std::vector<ParticlePointer> _particle_container;
  /**
   * @brief Underlying particle container. Stores unique pointers managing
   * ParticlePair objects.
   */
  std::unordered_set<ParticlePairPointer> _particle_pair_set;
};