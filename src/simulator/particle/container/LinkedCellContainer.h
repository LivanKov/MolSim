#include "../BoundaryCondition.h"
#include "DirectSumContainer.h"
#include "utils/logger/Logger.h"
#include <array>
#include <initializer_list>
#include <unordered_set>
#include <unordered_map>

#pragma once

#define DIVISION_TOLERANCE 1e-6

/**
 * @class LinkedCellContainer
 * @brief Class that provides a container for particles that uses linked cells
 * to speed up the computation.Inherits from DirectSumContainer.
 */
class LinkedCellContainer {

  /** @struct Cell
 *    @brief Manages a vector of shared pointers to Particle objects.
 */
struct Cell {
  std::unordered_set<int>particle_ids;
  size_t size() const;
  void insert(int id);
  void remove(int id);
};

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
  LinkedCellContainer(
      std::initializer_list<double> domain_size, double r_cutoff,
      const DomainBoundaryConditions &boundary_conditions = {
          BoundaryCondition::Outflow, BoundaryCondition::Outflow,
          BoundaryCondition::Outflow, BoundaryCondition::Outflow,
          BoundaryCondition::Outflow, BoundaryCondition::Outflow});

  /**
   * @brief Default constructor for LinkedCellContainer.
   */
  LinkedCellContainer();
  /**
   * @brief Inserts a particle into the container.
   * @param p The particle to be inserted.
   */
  void insert(Particle &p, bool placement = false);

  bool is_within_domain(const std::array<double, 3> &position);

  void clear();

  void readjust();

  void reinitialize(DirectSumContainer &container);

  void reinitialize(std::vector<Particle> &particles);

  void reinitialize(std::vector<ParticlePointer> &particles);

  /**
   * @brief Updates the location of a particle within the container based on its
   * old position.
   * @param p The particle whose location is updated.
   * @param old_position The particle's previous position.
   */
  void update_particle_location(int particle_id,
                                const std::array<double, 3> &old_position);

  /**
   * @brief Retrieves neighboring particles of a given particle within the
   * cutoff radius.
   * @param p The particle for which neighbors are retrieved.
   * @return A vector of shared pointers to neighboring particles.
   */
  std::vector<ParticlePointer> get_neighbours(int particle_id);

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

  /**
   * @brief The size of the simulation domain.
   */
  std::vector<double> domain_size_;

  /**
   * @brief The coordinates of the domain's lower left corner.
   */
  std::vector<double> left_corner_coordinates;

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

  std::vector<Cell> cells;
  DirectSumContainer particles;
  double r_cutoff_x;
  double r_cutoff_y;
  double r_cutoff_z;
  Logger &logger = Logger::getInstance();

  size_t size();

  Particle &operator[](size_t index);

   std::unordered_map<int,ParticlePointer> cells_map;

private:
  void readjust_coordinates(std::array<double, 3> current_low_left,
                            std::array<double, 3> current_up_right);

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