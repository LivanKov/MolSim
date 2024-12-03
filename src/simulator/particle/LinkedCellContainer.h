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
   * @brief Constructor for LinkedCellContainer.
   * @param domain_size The size of the simulation domain (e.g., {x, y, z}
   * dimensions).
   * @param r_cutoff The cutoff radius for interactions.
   * @param left_corner_coordinates The coordinates of the domain's lower left
   * corner.
   * @param boundary_conditions The boundary conditions for the simulation
   * domain.
   * @throws std::invalid_argument If the domain size is not 2D or 3D.
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

  /**
   * @brief Inserts a particle into the container.
   * @param p The particle to be inserted.
   */
  void insert(Particle &p) override;

  /**
   * @brief Updates the location of a particle within the container based on its
   * old position.
   * @param p The particle whose location is updated.
   * @param old_position The particle's previous position.
   */
  void update_particle_location(Particle &p,
                                std::array<double, 3> &old_position);

  /**
   * @brief Retrieves neighboring particles of a given particle within the
   * cutoff radius.
   * @param p The particle for which neighbors are retrieved.
   * @return A vector of shared pointers to neighboring particles.
   */
  std::vector<ParticlePointer> get_neighbours(Particle &p);

  /**
   * @brief Applies the specified boundary conditions to a particle.
   * @param p The particle to which boundary conditions are applied.
   */
  void handleBoundaryConditions(Particle &p);

  /**
   * @brief Removes particles that have crossed outflow boundaries.
   */
  void removeOutflowParticles();

  /**
   * @brief Updates particle positions and handles boundary conditions.
   * This function integrates particle updates and boundary management.
   */
  void updateParticles();

private:
  /**
   * @brief The size of the simulation domain.
   */
  const std::vector<double> domain_size_;

  /**
   * @brief The coordinates of the domain's lower left corner.
   */
  const std::vector<double> left_corner_coordinates;

  /**
   * @brief The cutoff radius for particle interactions.
   */
  double r_cutoff_;

  /**
   * @brief Number of cells along the x-dimension.
   */
  size_t x;

  /**
   * @brief Number of cells along the y-dimension.
   */
  size_t y;

  /**
   * @brief Number of cells along the z-dimension.
   */
  size_t z;

  /**
   * @brief A vector of all cells in the domain, stored unwrapped for easy
   * indexing.
   */
  std::vector<Cell> unwrapped_cells_;

  /**
   * @brief The boundary conditions for the simulation domain.
   */
  DomainBoundaryConditions boundary_conditions_;

  /**
   * @brief Retrieves a cell by its index in the unwrapped cell array.
   * @param index The index of the cell.
   * @return A reference to the cell at the specified index.
   */
  Cell &get_cell(size_t index);
};