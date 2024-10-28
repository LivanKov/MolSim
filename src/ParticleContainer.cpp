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
        "Pair requires to different instances of Particle");
  }
}

bool ParticlePair::operator==(const ParticlePair &rhs) const
{
  return first == rhs.first && second == rhs.second || first == rhs.second && second == rhs.first;
}

std::size_t ParticlePairHash::operator()(const ParticlePair &p) const
{
  std::size_t hash1 = std::hash<Particle>()(p.first) ^ (std::hash<Particle>()(p.second) << 1);
  std::size_t hash2 = std::hash<Particle>()(p.second) ^ (std::hash<Particle>()(p.first) << 1);
  return hash1 ^ hash2;
}

template <typename... Args>
  requires std::constructible_from<Particle, Args...>
void ParticleContainer::emplace_back(Args... args)
{
  particle_container.emplace_back(std::forward<args>...);
}

void ParticleContainer::insert(Particle &p) { emplace_back(p); }

void ParticleContainer::insert(Particle &&p) { emplace_back(p); }

size_t ParticleContainer::size() { return particle_container.size(); }
