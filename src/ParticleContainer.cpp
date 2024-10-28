#include "ParticleContainer.h"
#include <utility>

template <typename... Args> void ParticleContainer::emplace_back(Args... args) {
  particle_container.emplace_back(std::forward<args>...);
}

void ParticleContainer::insert(Particle &p) { emplace_back(p); }

void ParticleContainer::insert(Particle &&p) { emplace_back(p); }

size_t ParticleContainer::size() { return particle_container.size(); }

