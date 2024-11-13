/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include "utils/ArrayUtils.h"
#include <iostream>

#include "logger/Logger.h"

Particle::Particle(int type_arg) {
  type = type_arg;
  Logger::getInstance().trace("Particle generated!");
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
}

Particle::Particle(const Particle &other) {
  x = other.x;
  v = other.v;
  f = other.f;
  old_f = other.old_f;
  m = other.m;
  type = other.type;
  Logger::getInstance().trace("Particle generated by copy!");
}

// Todo: maybe use initializater list instead of copy?
Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, int type_arg) {
  x = x_arg;
  v = v_arg;
  m = m_arg;
  type = type_arg;
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
  Logger::getInstance().trace("Particle generated!");
}

Particle::~Particle() {
  // std::cout << "Particle destructed!" << std::endl;
  // Logger::getInstance().info("Particle destructed!");
}

const std::array<double, 3> &Particle::getX() const { return x; }

const std::array<double, 3> &Particle::getV() const { return v; }

const std::array<double, 3> &Particle::getF() const { return f; }

const std::array<double, 3> &Particle::getOldF() const { return old_f; }

double Particle::getM() const { return m; }

int Particle::getType() const { return type; }

std::string Particle::toString() const {
  std::stringstream stream;
  stream << "Particle: X:" << x << " v: " << v << " f: " << f
         << " old_f: " << old_f << " type: " << type;
  return stream.str();
}

bool Particle::operator==(const Particle &other) const {
  return (x == other.x) and (v == other.v) and (f == other.f) and
         (type == other.type) and (m == other.m) and (old_f == other.old_f);
}

size_t std::hash<Particle>::operator()(const Particle &p) const noexcept {
  size_t h1 = std::hash<double>{}(p.getM());
  size_t h2 = std::hash<int>{}(p.getType());
  return h1 ^ (h2 << 1);
}

std::ostream &operator<<(std::ostream &stream, Particle &p) {
  stream << p.toString();
  return stream;
}

void Particle::updateX(double x_arg, double y_arg, double z_arg) {
  x = {x_arg, y_arg, z_arg};
}

void Particle::updateX(const std::array<double, 3> &position) {
  x = position;
}

void Particle::updateV(double x_arg, double y_arg, double z_arg) {
  v = {x_arg, y_arg, z_arg};
}

void Particle::updateV(const std::array<double, 3> &velocity) {
  v = velocity;
}

void Particle::updateF(double x_arg, double y_arg, double z_arg) {
  f = {x_arg, y_arg, z_arg};
}

void Particle::updateF(const std::array<double, 3> &force) {
  f = force;
}

void Particle::updateOldF(double x_arg, double y_arg, double z_arg) {
  old_f = {x_arg, y_arg, z_arg};
}

void Particle::updateOldF(const std::array<double, 3> &force) {
  old_f = force;
}
