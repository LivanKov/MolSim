#include "ParticleContainer.h"
#include <array>
#include <initializer_list>

#pragma once


/**
 * @struct Cell
 * @brief Manages a vector of shared pointers to Particle objects.
 */
struct Cell{
  std::vector<ParticlePointer> particles;
  size_t size() const;
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
  LinkedCellContainer(std::initializer_list<double> domain_size, double r_cutoff,
                      std::initializer_list<double> left_corner_coordinates);

  void insert(Particle& p) override; 

  void update_particle_location(Particle& p, std::array<double, 3> &old_position);

  std::vector<ParticlePointer> get_neighbours(Particle& p);

  const std::vector<double> domain_size_;
  const std::vector<double> left_corner_coordinates;
  double r_cutoff_;
  size_t x;
  size_t y;
  size_t z;
  std::vector<Cell> cells;

  Cell& get_cell(size_t index);
  };