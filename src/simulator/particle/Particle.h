/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <memory>
#include <string>
#include <vector>

/**
 * @class Particle
 * @brief Represents a particle with position, velocity, force, and mass.
 *
 * The Particle class encapsulates the properties and behaviors of a particle
 * in a simulation, including its position, velocity, force, and mass. It also
 * includes methods for updating these properties and comparing particles.
 */
class Particle {

private:
  /**
   * Position of the particle
   */
  std::array<double, 3> x;

  /*
   * Old position of the particle
   */
  std::array<double, 3> old_x;

  /**
   * Velocity of the particle
   */
  std::array<double, 3> v;

  /**
   * Force effective on this particle
   */
  std::array<double, 3> f;

  /**
   * Force which was effective on this particle
   */
  std::array<double, 3> old_f;

  /**
   * Thermal motion of this particle
   */
  std::array<double, 3> thermal_motion_;

  /**
   * Kinetic motion of this particle
   */
  std::array<double, 3> kinetic_motion_;

  /**
   * Mass of this particle
   */
  double m;

  /**
   * Type of the particle. Use it for whatever you want (e.g. to separate
   * molecules belonging to different bodies, matters, and so on)
   */
  int type;

  /**
   * @brief Lennard-Jones potential parameter epsilon
   */
  double epsilon;

  /**
   * @brief Lennard-Jones potential parameter sigma
   */
  double sigma;

  /**
   * @brief Flag to apply additional force to this particle
   */
  bool apply_fzup;
  /*
   * @brief check if particle is fixed.
   */
  bool fixed;

public:
  /**
   * @brief Constructor.
   * @param type: optional parameter to define the molecule type.
   */
  explicit Particle(int type = 0);

  /**
   * @brief Copy constructor.
   * @param other: lvalue refernce to another Particle object.
   */

  Particle(const Particle &other);

  /**
   * @brief Constructor, intializes the particle with the given position,
   * velocity, mass, and type (optional).
   * @param x_arg, v_arg, m_arg, type
   */

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
      int type, double epsilon_arg = 5.0, double sigma_arg = 1.0,
      bool fixed = false);

  /**
   * @brief Destructor
   */
  virtual ~Particle();

  /**
   * @brief access the array containing the position of the particle.
   * @return a reference to the array containing the position of the particle.
   */

  const std::array<double, 3> &getX() const;

  /**
   * @brief access the array containing the velocity of the particle.
   * @return a reference to the array containing the velocity of the particle.
   */

  const std::array<double, 3> &getV() const;

  /**
   * @brief access the array containing the force of the particle.
   * @return a reference to the array containing the force of the particle.
   */

  const std::array<double, 3> &getF() const;

  /**
   * @brief access the array containing the force of the particle in the
   * previous step.
   * @return a reference to the array containing force of the particle in the
   * previous step.
   */

  const std::array<double, 3> &getOldF() const;

  /**
   * @brief access the array containing the thermal motion of the particle.
   * @return a reference to the array containing the thermal motion of the
   * particle.
   */
  const std::array<double, 3> &getThermalMotion() const;

  /**
   * @brief access the array containing the kinetic motion of the particle.
   * @return a reference to the array containing the kinetic motion of the
   * particle.
   */
  const std::array<double, 3> &getKineticMotion() const;

  /**
   * @brief access the array containing the old position of the particle.
   * @return a reference to the array containing the old position of the
   * particle.
   */

  const std::array<double, 3> &getOldX() const;

  /**
   * @brief returns the value, that correponds to particle mmass.
   * @return double variable containing the mass of the particle.
   */
  double getM() const;

  /**
   * @brief returns the value, that correponds to particle mmass.
   * @return integer variable containing the mass of the particle.
   */

  int getId() const;

  /**
   * @brief returns the value, that correponds to particle epsilon.
   * @return double variable containing the epsilon of the particle.
   */

  double getEpsilon() const;

  /**
   * @brief returns the value, that correponds to particle sigma.
   * @return double variable containing the sigma of the particle.
   */
  double getSigma() const;

  /**
   * @brief returns the value, that correponds to particle apply_fzup.
   * @return bool variable containing the apply_fzup of the particle.
   */

  bool isApplyFZup() const;

  /**
   * @brief set the value, that correponds to particle apply_fzup.
   * @param apply_fzup_arg: bool variable containing the apply_fzup of the
   * particle.
   */

  void setAppliyFZup(bool apply_fzup_arg);

  /**
   * @brief returns the value, that correponds to particle fixed.
   * @return bool variable containing the fixed of the particle.
   */

  bool is_fixed() const;

  /**
   * @brief set the fixed variable
   * @param fixed_arg: bool variable containing new value of fixed variable.
   */

  void setFixed(bool fixed_arg);

  /**
   * @brief outbound flag, used for boundary conditions.
   */

  bool outbound;

  /**
   * @brief overload the equality (==) operator to compare particles.
   * @param other: lvalue reference to another Particle object.
   * @return true/false corresponding to whether or not the function operator
   * and the operand are equal .
   */

  bool operator==(const Particle &other) const;
  /**
   * @brief overload the inequality (!=) operator to compare particles.
   * @param other: lvalue reference to another Particle object.
   * @return true/false corresponding to whether or not the function operator
   * and the operand are not equal .
   */

  bool operator!=(const Particle &other) const;

  /**
   * @brief returns string representation of the particle.
   * @return std::string representing the particle and its attributes.
   */

  std::string toString() const;

  /**
   * @brief updates the position of the particle.
   * @param x_arg, y_arg, z_arg: new position of the particle.
   */

  void updateX(double x_arg, double y_arg, double z_arg);

  /**
   * @brief updates the position of the particle.
   * @param position: allow the method to accept array.
   */
  void updateX(const std::array<double, 3> &position);

  /**
   * @brief updates the velocity of the particle.
   * @param x_arg, y_arg, z_arg: new velocity of the particle.
   */

  void updateV(double x_arg, double y_arg, double z_arg);

  /**
   * @brief updates the velocity of the particle.
   * @param velocity: allow the method to accept array.
   */
  void updateV(const std::array<double, 3> &velocity);

  /**
   * @brief updates the force of the particle.
   * @param x_arg, y_arg, z_arg: new force of the particle.
   */

  void updateF(double x_arg, double y_arg, double z_arg);

  /**
   * @brief updates the force of the particle.
   * @param force: allowed the method to accept a std::array<double, 3>.
   */
  void updateF(const std::array<double, 3> &force);

  /**
   * @brief updates the force of the particle in the previous step.
   * @param x_arg, y_arg, z_arg: new force of the particle in the previous step.
   */

  void updateOldF(double x_arg, double y_arg, double z_arg);

  /**
   * @brief updates the force of the particle in the previous step.
   * @param force: allowed the method to accept a std::array<double, 3>.
   */
  void updateOldF(const std::array<double, 3> &force);

  /**
   * @brief updates thermal motion of particle.
   * @param x_arg, y_arg, z_arg: new thermal motion values.
   */
  void updateThermalMotion(double x_arg, double y_arg, double z_arg);

  /**
   * @brief updates thermal motion of particle.
   * @param thermal_m: allowed the method to accept a std::array<double, 3>.
   */
  void updateThermalMotion(const std::array<double, 3> &thermal_m);

  /**
   * @brief updates kinetic motion of particle.
   * @param x_arg, y_arg, z_arg: new kinetic motion values.
   */
  void updateKineticMotion(double x_arg, double y_arg, double z_arg);

  /**
   * @brief updates kinetic motion of particle.
   * @param kinetic_m: allowed the method to accept a std::array<double, 3>.
   */
  void updateKineticMotion(const std::array<double, 3> &kinetic_m);

  /**
   * @brief updates the old position of the particle.
   * @param x_arg, y_arg, z_arg: new old position of the particle.
   */

  void updateOldX(const std::array<double, 3> &position);

  /**
   * @brief updates the old position of the particle.
   * @param x_arg, y_arg, z_arg: new old position of the particle.
   */

  void updateOldX(double x_arg, double y_arg, double z_arg);

  /**
   * @brief check if particle is outside the domain.
   */

  bool left_domain;

  /**
   * @brief membrane neighbours.
   */

  std::vector<std::shared_ptr<Particle>> membrane_neighbours;

  /**
   * @brief diagonal membrane neighbours.
   */

  std::vector<std::shared_ptr<Particle>> diagonal_membrane_neighbours;

  /**
   * @brief index of the cell the particle is in.
   */

  size_t cell_index;
};

/**
 * @brief overload the output stream operator to print the particle.
 * @param stream: lvalue reference to the output stream.
 * @param p: lvalue reference to the Particle object.
 * @return lvalue reference to the output stream.
 */

std::ostream &operator<<(std::ostream &stream, Particle &p);
