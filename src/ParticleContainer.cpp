#include "ParticleContainer.h"
#include "Particle.h"
#include <concepts>
#include <utility>

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

Particle& ParticleContainer::operator[](int index){
  return _particle_container[index];
}

ParticleContainer::PContainerIterator::PContainerIterator(PPointerType p) : _ptr{p} {}

ParticleContainer::PContainerIterator &ParticleContainer::PContainerIterator::operator++()
{
  _ptr++;
  return *this;
}

ParticleContainer::PContainerIterator ParticleContainer::PContainerIterator::operator++(int)
{
  PContainerIterator _ret = *this;
  ++(*this);
  return _ret;
}

ParticleContainer::PContainerIterator::PPointerType ParticleContainer::PContainerIterator::operator->()
{
  return _ptr;
}

bool ParticleContainer::PContainerIterator::operator==(ParticleContainer::PContainerIterator &rhs) const
{
  return _ptr == rhs._ptr;
}

bool ParticleContainer::PContainerIterator::operator!=(ParticleContainer::PContainerIterator &rhs) const
{
  return !(*this == rhs);
}

ParticleContainer::PContainerIterator::PReferenceType ParticleContainer::PContainerIterator::operator*() const
{
  return *_ptr;
}

ParticleContainer::PContainerIterator ParticleContainer::begin()
{
  return ParticleContainer::PContainerIterator(_particle_container.data());
}

ParticleContainer::PContainerIterator ParticleContainer::end()
{
  return ParticleContainer::PContainerIterator(_particle_container.data() + _particle_container.size());
}

void ParticleContainer::create_pairs(const Particle &p)
{

  for (const Particle &k : _particle_container)
  {
    if (!(k == p))
    {
      _particle_pair_map[p].emplace_back(k, p);
      _particle_pair_set.emplace(k,p);
    }
  }
}
