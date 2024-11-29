#include "ParticleContainer.h"
#include "Particle.h"
#include "utils/ArrayUtils.h"
#include <concepts>
#include <utility>

ParticlePair::ParticlePair(const ParticlePointer first,
                           const ParticlePointer second)
    : first{first}, second{second} {
  if (*first == *second) {
    throw std::invalid_argument(
        "Pair requires two different instances of Particle");
  }
}

bool ParticlePair::operator==(const ParticlePair &rhs) const {
  return (*first == *(rhs.first) && *second == *(rhs.second)) ||
         (*first == *(rhs.second) && *second == *(rhs.first));
}

std::string ParticlePair::toString() const {
  std::stringstream stream;
  stream << "ParticlePair: first: " << first->toString()
         << " second: " << second->toString();
  return stream.str();
}

std::size_t std::hash<ParticlePair>::operator()(const ParticlePair &p) const {
  std::size_t hash1 =
      std::hash<Particle>()(*p.first) ^ (std::hash<Particle>()(*p.second) << 1);
  std::size_t hash2 =
      std::hash<Particle>()(*p.second) ^ (std::hash<Particle>()(*p.first) << 1);
  return hash1 ^ hash2;
}

bool ParticlePointerEqual::operator()(
    const std::shared_ptr<Particle> &a,
    const std::shared_ptr<Particle> &b) const {
  return *a == *b;
}

size_t ParticlePointerHash::operator()(
    const std::shared_ptr<Particle> &s) const noexcept {
  return std::hash<Particle>()(*s);
}

ParticleContainer::ParticleContainer()
    : _particle_container{}, _particle_pair_container{} {}

void ParticleContainer::insert(Particle &p) { emplace_back(p); }

void ParticleContainer::insert(Particle &&p) { emplace_back(p); }

size_t ParticleContainer::size() { return _particle_container.size(); }

Particle &ParticleContainer::operator[](size_t index) {
  return *(_particle_container[index]);
}

ParticleIterator::ParticleIterator(PPointerType p) : _ptr(p) {}

ParticleIterator &ParticleIterator::operator++() {
  _ptr++;
  return *this;
}

ParticleIterator ParticleIterator::operator++(int) {
  ParticleIterator _ret = *this;
  ++(*this);
  return _ret;
}

ParticleIterator::PPointerType ParticleIterator::operator->() { return _ptr; }

bool ParticleIterator::operator==(const ParticleIterator &rhs) const {
  return _ptr == rhs._ptr;
}

bool ParticleIterator::operator!=(const ParticleIterator &rhs) const {
  return _ptr != rhs._ptr;
}

ParticleIterator::PReferenceType ParticleIterator::operator*() const {
  return **_ptr;
}

ParticleIterator ParticleContainer::begin() {
  return ParticleIterator(_particle_container.data());
}

ParticleIterator ParticleContainer::end() {
  return ParticleIterator(_particle_container.data() +
                          _particle_container.size());
}

ParticlePairIterator::ParticlePairIterator(PPointerType p) : _ptr{p} {}

ParticlePairIterator &ParticlePairIterator::operator++() {
  _ptr++;
  return *this;
}

ParticlePairIterator ParticlePairIterator::operator++(int) {
  ParticlePairIterator _ret = *this;
  ++(*this);
  return _ret;
}

ParticlePairIterator::PPointerType ParticlePairIterator::operator->() {
  return _ptr;
}

bool ParticlePairIterator::operator==(const ParticlePairIterator &rhs) const {
  return _ptr == rhs._ptr;
}

bool ParticlePairIterator::operator!=(const ParticlePairIterator &rhs) const {
  return _ptr != rhs._ptr;
}

ParticlePairIterator::PReferenceType ParticlePairIterator::operator*() const {
  return **_ptr;
}

void ParticleContainer::clear() {
  _particle_container.clear();
  _particle_pair_container.clear();
}

ParticlePairIterator ParticleContainer::pair_begin() {
  return ParticlePairIterator(_particle_pair_container.data());
}

ParticlePairIterator ParticleContainer::pair_end() {
  return ParticlePairIterator(_particle_pair_container.data() + _particle_pair_container.size());
}

void ParticleContainer::create_pairs(const ParticlePointer &new_particle) {
  for (auto const &p : _particle_container) {
    if (*new_particle != *p)
      _particle_pair_container.push_back(
          std::make_shared<ParticlePair>(new_particle, p));
  }
}
