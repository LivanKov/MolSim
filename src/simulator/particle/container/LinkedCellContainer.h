#include "DirectSumContainer.h"
#include "utils/logger/Logger.h"
#include <array>
#include <initializer_list>
#include <unordered_map>
#include <unordered_set>
#include <optional>

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

/**
 *@brief Enum class for cell placement
 */

enum Placement {
  TOP,
  BOTTOM,
  LEFT,
  RIGHT,
  FRONT,
  BACK,
  TOP_RIGHT_CORNER_BACK,
  TOP_LEFT_CORNER_BACK,
  BOTTOM_RIGHT_CORNER_BACK,
  BOTTOM_LEFT_CORNER_BACK,
  RIGHT_FRONT,
  LEFT_FRONT,
  RIGHT_BACK,
  LEFT_BACK,
  TOP_FRONT,
  TOP_BACK,
  BOTTOM_FRONT,
  BOTTOM_BACK,
  RIGHT_TOP,
  RIGHT_BOTTOM,
  LEFT_TOP,
  LEFT_BOTTOM,
  // for 2d case
  TOP_RIGHT_CORNER,
  TOP_LEFT_CORNER,
  BOTTOM_RIGHT_CORNER,
  BOTTOM_LEFT_CORNER
};

/**
 * @struct GhostParticle
 * @brief Struct for ghost particle, stores the values necessary for the calculation
 */

struct GhostParticle {
  double sigma;
  double epsilon;
  std::array<double, 3> position;
  ParticlePointer ptr;
  int id;

  bool operator==(GhostParticle &other) const;
};

/**
 * @class LinkedCellContainer
 * @brief Class that provides a container for particles that uses linked cells
 * to speed up the computation.Inherits from DirectSumContainer.
 */
class LinkedCellContainer {

  /** @struct Cell
   *  @brief Manages a collection of particles within a spatial cell of the linked cell structure.
   *  
   *  This struct represents a single cell in the linked cell data structure, which is used
   *  to optimize particle interaction calculations by spatial partitioning. Each cell
   *  maintains a set of particle IDs that fall within its spatial boundaries.
   */
  struct Cell {
    /** @brief Set of particle IDs contained in this cell. */
    std::unordered_set<int> particle_ids;
    
    /** @brief Returns the number of particles in this cell.
     *  @return Number of particles in the cell.
     */
    size_t size() const;
    
    /** @brief Adds a particle to this cell.
     *  @param id The ID of the particle to add.
     */
    void insert(int id);
    
    /** @brief Removes a particle from this cell.
     *  @param id The ID of the particle to remove.
     */
    void remove(int id);
    
    /** @brief Flag indicating if this is a halo cell.
     *  Halo cells are used for boundary condition calculations.
     */
    bool is_halo = false;
    
    /** @brief The placement of this cell relative to the domain boundaries. */
    Placement placement;
  };

public:
  /**
   * @brief Default constructor for LinkedCellContainer.
   */
  LinkedCellContainer();

  /**
   * @brief Initialize the container with domain size, cutoff radius, and
   * boundary conditions.
   * @param domain_size The size of the simulation domain.
   * @param r_cutoff The cutoff radius for interactions.
   * @param boundary_conditions The boundary conditions for the domain.
   */
  void initialize(const std::initializer_list<double> &domain_size,
                  double r_cutoff,
                  const DomainBoundaryConditions &boundary_conditions);

  /**
   * @brief Inserts a particle into the container.
   * @param p The particle to be inserted.
   */
  void insert(Particle &p, bool placement = false);

  /**
   * @brief Checks if a particle is within the simulation domain.
   * @param position The position of the particle.
   * @return True if the particle is within the domain, false otherwise.
   */
  bool is_within_domain(const std::array<double, 3> &position);

  /*
   * @brief Clears the container of all particles.
   */
  void clear();

  /**
   * @brief Adjusts the domain placement to ensure all particles are within the domain.
   */

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

  /**
   * @brief sets the boundary conditions for the domain
   * @param conditions the boundary conditions
   */

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

  /**
   * @brief A container of all particles in the domain.
   */
  DirectSumContainer particles;

  /**
   * @brief cell size for x,y,z
   * 
   */
  /** @brief The cutoff radius for each dimension.
   *  These values determine the size of each cell in the x, y, and z directions.
   */
  double r_cutoff_x;
  double r_cutoff_y;
  double r_cutoff_z;

  /** @brief Logger instance for debugging and error reporting. */
  Logger &logger = Logger::getInstance();

  /**
   * @brief Returns the total number of particles in the container.
   * @return Size of the container.
   */
  size_t size();

  /**
   * @brief Operator overload for accessing particles by index.
   * @param index The index of the particle to access.
   * @return Reference to the particle at the given index.
   */
  Particle &operator[](size_t index);

  /**
   * @brief Maps particle IDs to their corresponding cell pointers.
   * Used for efficient particle lookup in the linked cell structure.
   */
  std::unordered_map<int, ParticlePointer> cells_map;

  /**
   * @brief Counter for particles that have left the simulation domain.
   */
  size_t particles_left_domain;
  
  /**
   * @brief Unique identifier for the next particle to be added.
   */
  size_t particle_id;

  /**
   * @brief Flag indicating if this container is a wrapper around another container.
   */
  bool is_wrapper;

  /**
   * @brief Count of halo cells in the container.
   * Halo cells are the boundary cells used for periodic boundary conditions.
   */
  size_t halo_count;

  /**
   * @brief The boundary conditions for the simulation domain.
   */
  DomainBoundaryConditions boundary_conditions_;

  /**
   * @brief Flag indicating if reflective boundary conditions are active.
   */
  bool reflective_flag;

  /**
   * @brief Flag indicating if periodic boundary conditions are active.
   */
  bool periodic_flag;

  /**
   * @brief Indices of cells that are marked as halo cells.
   */
  std::vector<size_t> halo_cell_indices;

  /**
   * @brief List of particle IDs that have moved outside the domain.
   */
  std::vector<int> particles_outbound;

  /**
   * @brief Maps cell IDs to their corresponding ghost particles.
   * Used for implementing periodic boundary conditions.
   */
  std::unordered_map<int, std::vector<GhostParticle>> cell_ghost_particles_map;

  /**
   * @brief Removes all ghost particles from the container.
   */
  void clear_ghost_particles();

  /**
   * @brief Creates ghost particles for a given particle in a specific cell.
   * @param particle_id ID of the particle to create ghosts for.
   * @param cell_index Index of the cell containing the particle.
   */
  void create_ghost_particles(int particle_id, int cell_index);

  /**
   * @brief Creates a single ghost particle with a position offset.
   * @param particle_id ID of the original particle.
   * @param position_offset The offset to apply to the ghost particle's position.
   * @return The created ghost particle.
   */
  GhostParticle
  create_ghost_particle(int particle_id,
                        const std::array<double, 3> &position_offset);

  /**
   * @brief Gets additional neighboring ghost particles for a given particle.
   * @param particle_id ID of the particle to find neighbors for.
   * @return Vector of ghost particles that are neighbors of the given particle.
   */
  std::vector<GhostParticle> get_additional_neighbour_indices(int particle_id);

private:
  /**
   * @brief Adjusts the coordinates of the domain based on new boundaries.
   * @param current_low_left The new lower left corner coordinates.
   * @param current_up_right The new upper right corner coordinates.
   */
  void readjust_coordinates(std::array<double, 3> current_low_left,
                            std::array<double, 3> current_up_right);

  /**
   * @brief Assigns halo status to cells at the border of the array
   */
  void mark_halo_cells();


  Placement determine_placement_3d(size_t index, size_t i, size_t j, size_t k);

  Placement determine_placement_2d(size_t index, size_t i, size_t j);


  void assign_placements();
};