#include "BoundaryCondition.h"
#include "ParticleContainer.h"
#include <array>
#include <initializer_list>

#pragma once

/**
 * @struct Cell
 * @brief Manages a vector of shared pointers to Particle objects.
 */
struct Cell {
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
public:
  LinkedCellContainer(std::initializer_list<double> domain_size,
                      double r_cutoff,
                      std::initializer_list<double> left_corner_coordinates,
                      const DomainBoundaryConditions &boundary_conditions)
      : domain_size_(domain_size),
        left_corner_coordinates(left_corner_coordinates), r_cutoff_(r_cutoff),
        boundary_conditions_(boundary_conditions) {
    if (domain_size.size() != 3 && domain_size.size() != 2) {
      throw std::invalid_argument("Domain size must have 2 or 3 elements");
    }
    x = static_cast<size_t>(domain_size_[0] / r_cutoff);
    y = static_cast<size_t>(domain_size_[1] / r_cutoff);
    z = domain_size.size() == 3
            ? static_cast<size_t>((domain_size_[2] / r_cutoff))
            : 1;
    unwrapped_cells_ = std::vector<Cell>((x + 1) * (y + 1) * (z + 1), Cell());
  }

  void insert(Particle &p) override;

  void update_particle_location(Particle &p,
                                std::array<double, 3> &old_position);

  std::vector<ParticlePointer> get_neighbours(Particle &p);

  void handleBoundaryConditions(Particle &p);

  void removeOutflowParticles();

  void updateParticles();

private:
  const std::vector<double> domain_size_;
  const std::vector<double> left_corner_coordinates;
  double r_cutoff_;
  size_t x;
  size_t y;
  size_t z;
  std::vector<Cell> unwrapped_cells_;
  DomainBoundaryConditions boundary_conditions_;

  Cell &get_cell(size_t index);
};