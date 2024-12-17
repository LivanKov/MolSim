#include "DirectSumContainer.h"
#include "utils/logger/Logger.h"
#include <array>
#include <initializer_list>
#include <unordered_map>
#include <unordered_set>

#pragma once

#define DIVISION_TOLERANCE 1e-6


/**
 *@brief Enum class for boundary conditions
 */

enum BoundaryCondition { Outflow, Reflecting, Periodic };

/**
 *@brief Struct for domain boundary conditions
 */

struct DomainBoundaryConditions {
  BoundaryCondition left, right;
  BoundaryCondition top, bottom;
  BoundaryCondition front, back;
};


enum Placement {
  TOP,
  BOTTOM,
  LEFT,
  RIGHT,
  FRONT,
  BACK,
  TOP_RIGHT_CORNER,
  TOP_LEFT_CORNER,
  BOTTOM_RIGHT_CORNER,
  BOTTOM_LEFT_CORNER
};


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
    std::unordered_set<int> particle_ids;
    size_t size() const;
    void insert(int id);
    void remove(int id);
    bool is_halo = false;
    Placement placement;
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
   * @brief Retrieves the index of a cell based on its coordinates.
   * @param position The coordinates of the cell.
   * @return The index of the cell in the unwrapped cell array.
   */
  size_t get_cell_index(const std::array<double, 3> &position) const;


  void set_boundary_conditions(DomainBoundaryConditions conditions);

  /**
   * @brief The size of the simulation domain.
   */
  std::vector<double> domain_size_;

  /**
   * @brief The coordinates of the domain's lower left corner.
   */
  std::vector<double> left_corner_coordinates;

  /**
   * @brief The boundary conditions for the simulation domain.
   */

  std::unordered_map<Placement, BoundaryCondition> placement_map;

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

  std::unordered_map<int, ParticlePointer> cells_map;

  size_t particles_left_domain;
  size_t particle_id;


  bool is_wrapper;

  size_t halo_count;

  /**
   * @brief The boundary conditions for the simulation domain.
   */
  DomainBoundaryConditions boundary_conditions_;

  bool reflective_flag;

  bool periodic_flag;

  std::vector<size_t> halo_cell_indices;

  std::unordered_set<int> particles_outbound;
    

private:
  void readjust_coordinates(std::array<double, 3> current_low_left,
                            std::array<double, 3> current_up_right);

  /**
   * @brief Retrieves a cell by its index in the unwrapped cell array.
   * @param index The index of the cell.
   * @return A reference to the cell at the specified index.
   */
  Cell &get_cell(size_t index);

  /**
   * @brief Assigns halo status to cells at the border of the array
   */
  void mark_halo_cells();
};