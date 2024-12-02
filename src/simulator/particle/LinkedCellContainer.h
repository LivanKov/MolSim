#include "ParticleContainer.h"
#include <array>
#include <initializer_list>

#pragma once

#define DIVISION_TOLERANCE 1e-6

/**
 * @struct Cell
 * @brief Manages a vector of shared pointers to Particle objects.
 */
struct Cell {
  std::vector<ParticlePointer> particles;
  size_t size() const;
  ParticlePointer operator[](size_t index);
};

/**
 * @class LinkedCellContainer
 * @brief Class that provides a container for particles that uses linked cells
 * to speed up the computation.Inherits from ParticleContainer.
 */
class LinkedCellContainer : public ParticleContainer {
  /**
   * @brief Constructor.
   * @param domain_size Domain size of the container.
   * @param r_cutoff Cutoff radius.
   */
public:
  LinkedCellContainer(std::initializer_list<double> domain_size,
                      double r_cutoff,
                      std::initializer_list<double> left_corner_coordinates);

  void insert(Particle &p) override;

  void update_particle_location(ParticlePointer p,
                                std::array<double, 3> &old_position);

  bool is_within_domain(const std::array<double, 3> &position);

  void clear() override;

  std::vector<ParticlePointer> get_neighbours(Particle &p);

  void readjust();

  void reinitialize(ParticleContainer& container);

  void reinitialize(std::vector<Particle>& particles);

  void reinitialize(std::vector<ParticlePointer>& particles);

  const std::vector<double> domain_size_;
  std::vector<double> left_corner_coordinates;
  double r_cutoff_;
  size_t x;
  size_t y;
  size_t z;
  std::vector<Cell> cells;
  bool extend_x;
  bool extend_y;
  bool extend_z;

  Cell &get_cell(size_t index);
private:
  void readjust_coordinates(std::array<double,3> current_low_left, std::array<double,3> current_up_right);
};