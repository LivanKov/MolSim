#include "ParticleContainer.h"
#include <array>

#pragma once


/**
 * @struct Cell
 * @brief Manages a vector of shared pointers to Particle objects.
 */
struct Cell{
  std::vector<ParticlePointer> particles;
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
  LinkedCellContainer(const std::vector<double> domain_size, double r_cutoff,
                      std::array<double, 3> left_corner_coordinates);

  void insert(Particle& p) override; 

  void update_particle_location(Particle& p, std::array<double, 3> &old_position);

  std::vector<ParticlePointer> get_neighbours(Particle& p);
  
  private:
  const std::vector<double> domain_size_;
  std::array<double, 3> left_corner_coordinates;
  double r_cutoff_;
  size_t x;
  size_t y;
  size_t z;
  std::vector<Cell> unwrapped_cells_;

  Cell& get_cell(size_t index);
  };