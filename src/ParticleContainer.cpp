#include "ParticleContainer.h"
#include "Particle.h"
#include <concepts>
#include <utility>
#include "utils/ArrayUtils.h"

ParticlePair::ParticlePair(const Particle &first, const Particle &second)
    : first{first}, second{second}
{
  if (&first == &second)
  {
    throw std::invalid_argument(
        "Pair requires two different instances of Particle");
  }
}

bool ParticlePair::operator==(const ParticlePair &rhs) const
{
  return (first == rhs.first && second == rhs.second) || (first == rhs.second && second == rhs.first);
}

std::string ParticlePair::toString() const
{
  std::stringstream stream;
  stream << "ParticlePair: first: " << first.toString() << " second: " << second.toString();
  return stream.str();
}

std::size_t std::hash<ParticlePair>::operator()(const ParticlePair &p) const
{
  std::size_t hash1 = std::hash<Particle>()(p.first) ^ (std::hash<Particle>()(p.second) << 1);
  std::size_t hash2 = std::hash<Particle>()(p.second) ^ (std::hash<Particle>()(p.first) << 1);
  return hash1 ^ hash2;
}

ParticleContainer::ParticleContainer() : _particle_container{}, _particle_pair_set{}, _particle_pair_map{} {}

void ParticleContainer::insert(Particle &p)
{
  emplace_back(p);
}

void ParticleContainer::insert(Particle &&p)
{
  emplace_back(p);
}

size_t ParticleContainer::size() { return _particle_container.size(); }

std::vector<ParticlePair> &ParticleContainer::pairs_of(const Particle &p)
{
  return _particle_pair_map[p];
}

std::vector<ParticlePair> &ParticleContainer::operator[](const Particle &p)
{
  return _particle_pair_map[p];
}

Particle &ParticleContainer::operator[](int index)
{
  return _particle_container[index];
}

ParticleIterator::ParticleIterator(PPointerType p) : _ptr{p} {}

ParticleIterator &ParticleIterator::operator++()
{
  _ptr++;
  return *this;
}

ParticleIterator ParticleIterator::operator++(int)
{
  ParticleIterator _ret = *this;
  ++(*this);
  return _ret;
}

ParticleIterator::PPointerType ParticleIterator::operator->()
{
  return _ptr;
}

bool ParticleIterator::operator==(ParticleIterator &rhs) const
{
  return _ptr == rhs._ptr;
}

bool ParticleIterator::operator!=(ParticleIterator &rhs) const
{
  return !(*this == rhs);
}

ParticleIterator::PReferenceType ParticleIterator::operator*() const
{
  return *_ptr;
}

ParticleIterator ParticleContainer::begin()
{
  return ParticleIterator(_particle_container.data());
}

ParticleIterator ParticleContainer::end()
{
  return ParticleIterator(_particle_container.data() + _particle_container.size());
}

ParticlePairIterator::ParticlePairIterator(PPointerType p) : _ptr{p} {}

ParticlePairIterator &ParticlePairIterator::operator++()
{
  _ptr++;
  return *this;
}

ParticlePairIterator ParticlePairIterator::operator++(int)
{
  ParticlePairIterator _ret = *this;
  ++(*this);
  return _ret;
}

ParticlePairIterator::PPointerType
ParticlePairIterator::operator->()
{
  return _ptr;
}

bool ParticlePairIterator::operator==(ParticlePairIterator &rhs) const
{
  return _ptr == rhs._ptr;
}

bool ParticlePairIterator::operator!=(ParticlePairIterator &rhs) const
{
  return !(*this == rhs);
}

ParticlePairIterator::PReferenceType ParticlePairIterator::operator*() const
{
  return *_ptr;
}

ParticlePairIterator ParticleContainer::pair_begin()
{
  return ParticlePairIterator(_particle_pair_set.cbegin());
}

ParticlePairIterator ParticleContainer::pair_end()
{
  return ParticlePairIterator(_particle_pair_set.cend());
}

void ParticleContainer::create_pairs(const Particle &new_particle)
{
  for (const Particle &k : _particle_container)
  {
    _particle_pair_set.emplace(k, new_particle);
    _particle_pair_map[k].emplace_back(k, new_particle);
  }
}
