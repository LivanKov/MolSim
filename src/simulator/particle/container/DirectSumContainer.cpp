#include "DirectSumContainer.h"
#include "../Particle.h"
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

DirectSumContainer::DirectSumContainer()
    : _particle_container{}, _particle_pair_container{} {}

void DirectSumContainer::insert(Particle &p) {
  ParticlePointer p_ptr = std::make_shared<Particle>(p);
  _particle_container.push_back(p_ptr);
  create_pairs(p_ptr);
}

void DirectSumContainer::insert(ParticlePointer &p) {
  _particle_container.push_back(p);
  create_pairs(p);
}

size_t DirectSumContainer::size() { return _particle_container.size(); }

Particle &DirectSumContainer::operator[](size_t index) {
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

ParticleIterator DirectSumContainer::begin() {
  return ParticleIterator(_particle_container.data());
}

ParticleIterator DirectSumContainer::end() {
  return ParticleIterator(_particle_container.data() +
                          _particle_container.size());
}

void DirectSumContainer::clear() {
  _particle_container.clear();
  _particle_pair_container.clear();
}

std::vector<ParticlePair>::iterator DirectSumContainer::pair_begin() {
  return _particle_pair_container.begin();
}

std::vector<ParticlePair>::iterator DirectSumContainer::pair_end() {
  return _particle_pair_container.end();
}

void DirectSumContainer::create_pairs(const ParticlePointer &new_particle) {
  for (auto const &p : _particle_container) {
    if (*new_particle != *p)
      _particle_pair_container.push_back(ParticlePair(new_particle, p));
  }
}