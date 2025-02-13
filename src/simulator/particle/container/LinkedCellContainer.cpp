#include "LinkedCellContainer.h"
#include "../../Simulation.h"
#include "io/input/cli/SimParams.h"
#include <cmath>
#include <iostream>

size_t LinkedCellContainer::Cell::size() const { return particle_ids.size(); }

void LinkedCellContainer::Cell::insert(int id) { particle_ids.push_back(id); }

void LinkedCellContainer::Cell::remove(int id) {
  particle_ids.erase(std::remove(particle_ids.begin(), particle_ids.end(), id),
                     particle_ids.end());
}

void LinkedCellContainer::initialize(
    const std::vector<double> &domain_size, double r_cutoff,
    const DomainBoundaryConditions &boundary_conditions) {
  if (SimParams::fixed_Domain) {
    left_corner_coordinates = {SimParams::lower_left_corner[0],
                               SimParams::lower_left_corner[1],
                               SimParams::lower_left_corner[2]};
  }

  logger.info("Initializing LinkedCellContainer");
  domain_size_.assign(domain_size.begin(), domain_size.end());

  r_cutoff_ = r_cutoff;
  boundary_conditions_ = boundary_conditions;
  r_cutoff_x = r_cutoff;
  r_cutoff_y = r_cutoff;
  r_cutoff_z = r_cutoff;

  if (domain_size.size() != 3 && domain_size.size() != 2) {
    throw std::invalid_argument("Domain size must have 2 or 3 elements");
  }

  double remainder_x = std::fmod(domain_size_[0], r_cutoff);
  double remainder_y = std::fmod(domain_size_[1], r_cutoff);
  double remainder_z =
      domain_size.size() == 3 ? std::fmod(domain_size_[2], r_cutoff) : 0.0;
  if (std::abs(remainder_x) > DIVISION_TOLERANCE)
    r_cutoff_x += remainder_x / (std::floor(domain_size_[0] / r_cutoff));

  if (std::abs(remainder_y) > DIVISION_TOLERANCE)
    r_cutoff_y += remainder_y / (std::floor(domain_size_[1] / r_cutoff));

  if (domain_size.size() == 3) {
    if (std::abs(remainder_z) > DIVISION_TOLERANCE)
      r_cutoff_z += remainder_z / (std::floor(domain_size_[2] / r_cutoff));
  }

  x = static_cast<size_t>(domain_size_[0] / r_cutoff);
  y = static_cast<size_t>(domain_size_[1] / r_cutoff);
  z = domain_size.size() == 3
          ? (static_cast<size_t>((domain_size_[2] / r_cutoff)))
          : 1;
  cells = std::vector<Cell>(x * y * z, Cell());

  placement_map[Placement::TOP] = boundary_conditions.top;
  placement_map[Placement::BOTTOM] = boundary_conditions.bottom;
  placement_map[Placement::LEFT] = boundary_conditions.left;
  placement_map[Placement::RIGHT] = boundary_conditions.right;
  placement_map[Placement::FRONT] = boundary_conditions.front;
  placement_map[Placement::BACK] = boundary_conditions.back;

  if (domain_size.size() == 2) {
    std::cout << "Domain size: " << domain_size.size() << std::endl;
    std::cout << "Domain size_1: " << domain_size_.size() << std::endl;

    if (placement_map[Placement::TOP] == placement_map[Placement::RIGHT])
      placement_map[Placement::TOP_RIGHT_CORNER] =
          placement_map[Placement::RIGHT];
    else
      placement_map[Placement::TOP_RIGHT_CORNER] =
          placement_map[Placement::TOP];

    if (placement_map[Placement::TOP] == placement_map[Placement::LEFT])
      placement_map[Placement::TOP_LEFT_CORNER] =
          placement_map[Placement::LEFT];
    else
      placement_map[Placement::TOP_LEFT_CORNER] = placement_map[Placement::TOP];

    if (placement_map[Placement::BOTTOM] == placement_map[Placement::RIGHT])
      placement_map[Placement::BOTTOM_RIGHT_CORNER] =
          placement_map[Placement::RIGHT];
    else
      placement_map[Placement::BOTTOM_RIGHT_CORNER] =
          placement_map[Placement::BOTTOM];

    if (placement_map[Placement::BOTTOM] == placement_map[Placement::LEFT])
      placement_map[Placement::BOTTOM_LEFT_CORNER] =
          placement_map[Placement::LEFT];
    else
      placement_map[Placement::BOTTOM_LEFT_CORNER] =
          placement_map[Placement::BOTTOM];
  } else if (domain_size.size() == 3) {
    // initialize 3D edges
    placement_map[Placement::TOP_FRONT_EDGE] =
        (boundary_conditions.top == boundary_conditions.front)
            ? boundary_conditions.front
            : boundary_conditions.top;

    placement_map[Placement::TOP_BACK_EDGE] =
        (boundary_conditions.top == boundary_conditions.back)
            ? boundary_conditions.back
            : boundary_conditions.top;

    placement_map[Placement::BOTTOM_FRONT_EDGE] =
        (boundary_conditions.bottom == boundary_conditions.front)
            ? boundary_conditions.front
            : boundary_conditions.bottom;

    placement_map[Placement::BOTTOM_BACK_EDGE] =
        (boundary_conditions.bottom == boundary_conditions.back)
            ? boundary_conditions.back
            : boundary_conditions.bottom;

    placement_map[Placement::LEFT_FRONT_EDGE] =
        (boundary_conditions.left == boundary_conditions.front)
            ? boundary_conditions.front
            : boundary_conditions.left;

    placement_map[Placement::LEFT_BACK_EDGE] =
        (boundary_conditions.left == boundary_conditions.back)
            ? boundary_conditions.back
            : boundary_conditions.left;

    placement_map[Placement::RIGHT_FRONT_EDGE] =
        (boundary_conditions.right == boundary_conditions.front)
            ? boundary_conditions.front
            : boundary_conditions.right;

    placement_map[Placement::RIGHT_BACK_EDGE] =
        (boundary_conditions.right == boundary_conditions.back)
            ? boundary_conditions.back
            : boundary_conditions.right;

    placement_map[Placement::TOP_LEFT_EDGE] =
        (boundary_conditions.top == boundary_conditions.left)
            ? boundary_conditions.left
            : boundary_conditions.top;

    placement_map[Placement::TOP_RIGHT_EDGE] =
        (boundary_conditions.top == boundary_conditions.right)
            ? boundary_conditions.right
            : boundary_conditions.top;

    placement_map[Placement::BOTTOM_LEFT_EDGE] =
        (boundary_conditions.bottom == boundary_conditions.left)
            ? boundary_conditions.left
            : boundary_conditions.bottom;

    placement_map[Placement::BOTTOM_RIGHT_EDGE] =
        (boundary_conditions.bottom == boundary_conditions.right)
            ? boundary_conditions.right
            : boundary_conditions.bottom;

    // Initialize 3D corners (inferred from adjacent edges or faces)
    placement_map[Placement::TOP_FRONT_RIGHT_CORNER] =
        (boundary_conditions.top == boundary_conditions.front &&
         boundary_conditions.top == boundary_conditions.right)
            ? boundary_conditions.right
            : BoundaryCondition::Reflecting;

    placement_map[Placement::TOP_FRONT_LEFT_CORNER] =
        (boundary_conditions.top == boundary_conditions.front &&
         boundary_conditions.top == boundary_conditions.left)
            ? boundary_conditions.left
            : BoundaryCondition::Reflecting;

    placement_map[Placement::TOP_BACK_RIGHT_CORNER] =
        (boundary_conditions.top == boundary_conditions.back &&
         boundary_conditions.top == boundary_conditions.right)
            ? boundary_conditions.right
            : BoundaryCondition::Reflecting;

    placement_map[Placement::TOP_BACK_LEFT_CORNER] =
        (boundary_conditions.top == boundary_conditions.back &&
         boundary_conditions.top == boundary_conditions.left)
            ? boundary_conditions.left
            : BoundaryCondition::Reflecting;

    placement_map[Placement::BOTTOM_FRONT_RIGHT_CORNER] =
        (boundary_conditions.bottom == boundary_conditions.front &&
         boundary_conditions.bottom == boundary_conditions.right)
            ? boundary_conditions.right
            : BoundaryCondition::Reflecting;

    placement_map[Placement::BOTTOM_FRONT_LEFT_CORNER] =
        (boundary_conditions.bottom == boundary_conditions.front &&
         boundary_conditions.bottom == boundary_conditions.left)
            ? boundary_conditions.left
            : BoundaryCondition::Reflecting;

    placement_map[Placement::BOTTOM_BACK_RIGHT_CORNER] =
        (boundary_conditions.bottom == boundary_conditions.back &&
         boundary_conditions.bottom == boundary_conditions.right)
            ? boundary_conditions.right
            : BoundaryCondition::Reflecting;

    placement_map[Placement::BOTTOM_BACK_LEFT_CORNER] =
        (boundary_conditions.bottom == boundary_conditions.back &&
         boundary_conditions.bottom == boundary_conditions.left)
            ? boundary_conditions.left
            : BoundaryCondition::Reflecting;
  }

  mark_halo_cells();
}

LinkedCellContainer::LinkedCellContainer()
    : particle_id{0}, domain_size_{0, 0, 0}, left_corner_coordinates{0.0, 0.0,
                                                                     0.0},
      placement_map{}, r_cutoff_{0}, x{0}, y{0}, z{0}, particles_left_domain{0},
      is_wrapper{false}, halo_count{0}, boundary_conditions_{},
      reflective_flag{false}, periodic_flag{false}, halo_cell_indices{},
      particles_outbound{}, cell_ghost_particles_map{} {}

void LinkedCellContainer::insert(Particle &p, bool placement) {
  ParticlePointer p_ptr = std::make_shared<Particle>(p);
  if (placement && is_within_domain(p_ptr->getX())) {
    size_t index = get_cell_index(p_ptr->getX());
    cells[index].insert(p_ptr->getId());
  } else if (!is_within_domain(p_ptr->getX())) {
    p_ptr->left_domain = true;
    particles_left_domain++;
  }
  particles.insert(p_ptr);
}

bool LinkedCellContainer::is_within_domain(
    const std::array<double, 3> &position) {
  return position[0] >= left_corner_coordinates[0] &&
         position[0] < left_corner_coordinates[0] + domain_size_[0] &&
         position[1] >= left_corner_coordinates[1] &&
         position[1] < left_corner_coordinates[1] + domain_size_[1] &&
         position[2] >= left_corner_coordinates[2] &&
         (domain_size_.size() == 3
              ? position[2] < left_corner_coordinates[2] + domain_size_[2]
              : position[2] < left_corner_coordinates[2] + r_cutoff_z);
}

void LinkedCellContainer::update_particle_location(
    int particle_id, const std::array<double, 3> &old_position) {

  size_t old_index = get_cell_index(old_position);
  size_t current_index = get_cell_index(particles[particle_id].getX());

  bool current_within_domain = is_within_domain(particles[particle_id].getX());

  bool old_within_domain = is_within_domain(old_position);

  if (current_index != old_index ||
      current_within_domain != old_within_domain) {
    if (old_within_domain) {
      cells[old_index].remove(particle_id);
    }
    if (current_within_domain) {
      cells[current_index].insert(particle_id);

    } else {
      particles_outbound.push_back(particle_id);
    }
  }
}

std::vector<ParticlePointer>
LinkedCellContainer::get_neighbours(int particle_id) {
  std::vector<ParticlePointer> neighbours{};
  if (particles[particle_id].left_domain || particles[particle_id].outbound) {
    return neighbours;
  }
  std::array<double, 3> position = particles[particle_id].getX();
  int cell_index = get_cell_index(position);

  for (auto &i : cells[cell_index].particle_ids) {

    if (i != particle_id) {
      neighbours.push_back(particles.at(i));
    }
  }

  size_t i = static_cast<size_t>((position[0] - left_corner_coordinates[0]) /
                                 r_cutoff_x);
  size_t j = static_cast<size_t>((position[1] - left_corner_coordinates[1]) /
                                 r_cutoff_y);
  size_t k = domain_size_.size() == 3
                 ? static_cast<size_t>(
                       (position[2] - left_corner_coordinates[2]) / r_cutoff_z)
                 : 0;

  for (int di = -1; di <= 1; ++di) {
    for (int dj = -1; dj <= 1; ++dj) {
      for (int dk = -1; dk <= 1; ++dk) {
        if (di == 0 && dj == 0 && dk == 0) {
          continue;
        }
        int ni = i + di;
        int nj = j + dj;
        int nk = k + dk;

        if (ni >= 0 && static_cast<size_t>(ni) < x && nj >= 0 &&
            static_cast<size_t>(nj) < y && nk >= 0 &&
            static_cast<size_t>(nk) < z) {
          int neighborIndex = ni + (nj * x) + nk * x * y;
          for (auto &s : cells[neighborIndex].particle_ids) {
            neighbours.push_back(particles.at(s));
          }
        }
      }
    }
  }
  return neighbours;
}

std::vector<GhostParticle>
LinkedCellContainer::get_periodic_neighbours(int particle_id) {

  auto cell_index = get_cell_index(particles[particle_id].getX());

  std::vector<GhostParticle> ghost_neighbours =
      cell_ghost_particles_map[cell_index];

  return ghost_neighbours;
}

void LinkedCellContainer::clear() {
  for (size_t i = 0; i < cells.size(); ++i) {
    cells[i].particle_ids.clear();
  }
  particles.clear();
}

size_t LinkedCellContainer::get_cell_index(
    const std::array<double, 3> &position) const {
  size_t i = static_cast<size_t>((position[0] - left_corner_coordinates[0]) /
                                 r_cutoff_x);
  size_t j = static_cast<size_t>((position[1] - left_corner_coordinates[1]) /
                                 r_cutoff_y);
  size_t k = domain_size_.size() == 3
                 ? static_cast<size_t>(
                       (position[2] - left_corner_coordinates[2]) / r_cutoff_z)
                 : 0;
  return i + j * x + k * x * y;
}

size_t LinkedCellContainer::size() { return particles.size(); }

Particle &LinkedCellContainer::operator[](size_t index) {
  return particles[index];
}

ParticlePointer &LinkedCellContainer::at(size_t index) {
  return particles.at(index);
}

void LinkedCellContainer::mark_halo_cells() {
  for (size_t i = 0; i < x; i++) {
    for (size_t j = 0; j < y; j++) {
      for (size_t k = 0; k < z; k++) {
        if (i == 0 || i == x - 1 || j == 0 || j == y - 1 ||
            (z > 1 && (k == 0 || k == z - 1))) {

          size_t index = i + j * x + k * x * y;
          cells[index].is_halo = true;
          halo_cell_indices.push_back(index);
          halo_count++;

          if (domain_size_.size() == 2) {
            if (index == 0) {
              cells[index].placement = Placement::BOTTOM_LEFT_CORNER;
              cells[index].boundary_condition =
                  placement_map[BOTTOM_LEFT_CORNER];
            } else if (index == x - 1) {
              cells[index].placement = Placement::BOTTOM_RIGHT_CORNER;
              cells[index].boundary_condition =
                  placement_map[BOTTOM_RIGHT_CORNER];
            } else if (index == x * y - x) {
              cells[index].placement = Placement::TOP_LEFT_CORNER;
              cells[index].boundary_condition = placement_map[TOP_LEFT_CORNER];
            } else if (index == x * y - 1) {
              cells[index].placement = Placement::TOP_RIGHT_CORNER;
              cells[index].boundary_condition = placement_map[TOP_RIGHT_CORNER];
            } else if (index < x) {
              cells[index].placement = Placement::BOTTOM;
              cells[index].boundary_condition = placement_map[BOTTOM];
            } else if (index >= x * (y - 1)) {
              cells[index].placement = Placement::TOP;
              cells[index].boundary_condition = placement_map[TOP];
            } else if (index % x == 0) {
              cells[index].placement = Placement::LEFT;
              cells[index].boundary_condition = placement_map[LEFT];
            } else if ((index + 1) % x == 0) {
              cells[index].placement = Placement::RIGHT;
              cells[index].boundary_condition = placement_map[RIGHT];
            }
          } else if (domain_size_.size() == 3) {
            // Corner : Most specific
            if (i == 0 && j == 0 && k == 0) {
              cells[index].placement = Placement::BOTTOM_BACK_LEFT_CORNER;
              continue;
            }
            if (i == 0 && j == 0 && k == z - 1) {
              cells[index].placement = Placement::BOTTOM_FRONT_LEFT_CORNER;
              continue;
            }
            if (i == 0 && j == y - 1 && k == 0) {
              cells[index].placement = Placement::TOP_BACK_LEFT_CORNER;
              continue;
            }
            if (i == 0 && j == y - 1 && k == z - 1) {
              cells[index].placement = Placement::TOP_FRONT_LEFT_CORNER;
              continue;
            }
            if (i == x - 1 && j == 0 && k == 0) {
              cells[index].placement = Placement::BOTTOM_BACK_RIGHT_CORNER;
              continue;
            }
            if (i == x - 1 && j == 0 && k == z - 1) {
              cells[index].placement = Placement::BOTTOM_FRONT_RIGHT_CORNER;
              continue;
            }
            if (i == x - 1 && j == y - 1 && k == 0) {
              cells[index].placement = Placement::TOP_BACK_RIGHT_CORNER;
              continue;
            }
            if (i == x - 1 && j == y - 1 && k == z - 1) {
              cells[index].placement = Placement::TOP_FRONT_RIGHT_CORNER;
              continue;
            }

            // Edge: Less specific than corners
            if (i == 0 && j == 0) {
              cells[index].placement = Placement::BOTTOM_LEFT_EDGE;
              continue;
            }
            if (i == 0 && j == y - 1) {
              cells[index].placement = Placement::TOP_LEFT_EDGE;
              continue;
            }
            if (i == x - 1 && j == 0) {
              cells[index].placement = Placement::BOTTOM_RIGHT_EDGE;
              continue;
            }
            if (i == x - 1 && j == y - 1) {
              cells[index].placement = Placement::TOP_RIGHT_EDGE;
              continue;
            }
            if (j == 0 && k == 0) {
              cells[index].placement = Placement::BOTTOM_BACK_EDGE;
              continue;
            }
            if (j == 0 && k == z - 1) {
              cells[index].placement = Placement::BOTTOM_FRONT_EDGE;
              continue;
            }
            if (j == y - 1 && k == 0) {
              cells[index].placement = Placement::TOP_BACK_EDGE;
              continue;
            }
            if (j == y - 1 && k == z - 1) {
              cells[index].placement = Placement::TOP_FRONT_EDGE;
              continue;
            }
            if (i == 0 && k == 0) {
              cells[index].placement = Placement::LEFT_BACK_EDGE;
              continue;
            }
            if (i == 0 && k == z - 1) {
              cells[index].placement = Placement::LEFT_FRONT_EDGE;
              continue;
            }
            if (i == x - 1 && k == 0) {
              cells[index].placement = Placement::RIGHT_BACK_EDGE;
              continue;
            }
            if (i == x - 1 && k == z - 1) {
              cells[index].placement = Placement::RIGHT_FRONT_EDGE;
              continue;
            }

            // Face: Least specific
            if (k == 0) {
              cells[index].placement = Placement::BACK;
              continue;
            }
            if (k == z - 1) {
              cells[index].placement = Placement::FRONT;
              continue;
            }
            if (j == 0) {
              cells[index].placement = Placement::BOTTOM;
              continue;
            }
            if (j == y - 1) {
              cells[index].placement = Placement::TOP;
              continue;
            }
            if (i == 0) {
              cells[index].placement = Placement::LEFT;
              continue;
            }
            if (i == x - 1) {
              cells[index].placement = Placement::RIGHT;
              continue;
            }
          }
        }
      }
    }
  }
}

void LinkedCellContainer::readjust_coordinates(
    std::array<double, 3> current_low_left,
    std::array<double, 3> current_up_right) {
  std::array<double, 3> midpoint{};
  for (size_t i = 0; i < 3; ++i) {
    midpoint[i] = (current_low_left[i] + current_up_right[i]) / 2.0;
  }
  for (size_t i = 0; i < 3; ++i) {
    left_corner_coordinates[i] = midpoint[i] - domain_size_[i] / 2.0;
  }
}

void LinkedCellContainer::readjust() {
  std::array<double, 3> current_low_left = particles[0].getX();
  std::array<double, 3> current_up_right = particles[0].getX();
  for (auto &p : particles) {
    if (p.getX()[0] < current_low_left[0] ||
        p.getX()[1] < current_low_left[1] ||
        p.getX()[2] < current_low_left[2]) {
      current_low_left = p.getX();
    }
    if (p.getX()[0] > current_up_right[0] ||
        p.getX()[1] > current_up_right[1] ||
        p.getX()[2] > current_up_right[2]) {
      current_up_right = p.getX();
    }
  }
  readjust_coordinates(current_low_left, current_up_right);
  auto particles_ = particles;
  clear();
  for (auto &p : particles_) {
    insert(p, true);
  }
}

void LinkedCellContainer::set_boundary_conditions(

    DomainBoundaryConditions conditions) {
  this->boundary_conditions_ = conditions;
  placement_map[Placement::TOP] = conditions.top;
  placement_map[Placement::BOTTOM] = conditions.bottom;
  placement_map[Placement::LEFT] = conditions.left;
  placement_map[Placement::RIGHT] = conditions.right;
  placement_map[Placement::FRONT] = conditions.front;
  placement_map[Placement::BACK] = conditions.back;
}

void LinkedCellContainer::clear_ghost_particles() {
  cell_ghost_particles_map.clear();
}

GhostParticle LinkedCellContainer::create_ghost_particle(
    int particle_id, const std::array<double, 3> &position_offset) {
  GhostParticle ghost;
  ghost.sigma = particles[particle_id].getSigma();
  ghost.epsilon = particles[particle_id].getEpsilon();
  ghost.position = {particles[particle_id].getX()[0] + position_offset[0],
                    particles[particle_id].getX()[1] + position_offset[1],
                    particles[particle_id].getX()[2] + position_offset[2]};
  ghost.id = particle_id;
  ghost.ptr = particles.at(particle_id);
  return ghost;
}

ParticleIterator LinkedCellContainer::begin() { return particles.begin(); }

ParticleIterator LinkedCellContainer::end() { return particles.end(); }

void LinkedCellContainer::create_ghost_particles(int particle_id,
                                                 int cell_index) {
  const auto &placement = cells[cell_index].placement;
  if (domain_size_.size() == 2) {
    // Helper arrays for offsets
    const std::array<double, 3> right_offset = {domain_size_[0], 0, 0};
    const std::array<double, 3> left_offset = {-domain_size_[0], 0, 0};
    const std::array<double, 3> top_offset = {0, -domain_size_[1], 0};
    const std::array<double, 3> bottom_offset = {0, domain_size_[1], 0};

    switch (placement) {
    case Placement::LEFT: {
      auto ghost = create_ghost_particle(particle_id, right_offset);
      cell_ghost_particles_map[cell_index + x - 1].push_back(ghost);
      cell_ghost_particles_map[cell_index + x + x - 1].push_back(ghost);
      cell_ghost_particles_map[cell_index - 1].push_back(ghost);
      break;
    }

    case Placement::RIGHT: {
      auto ghost = create_ghost_particle(particle_id, left_offset);
      cell_ghost_particles_map[cell_index - x + 1].push_back(ghost);
      cell_ghost_particles_map[cell_index - (x + x - 1)].push_back(ghost);
      cell_ghost_particles_map[cell_index + 1].push_back(ghost);
      break;
    }
    case Placement::TOP: {
      auto ghost = create_ghost_particle(particle_id, top_offset);
      cell_ghost_particles_map[cell_index - (y - 1) * x].push_back(ghost);
      cell_ghost_particles_map[cell_index - ((y - 1) * x) - 1].push_back(ghost);
      cell_ghost_particles_map[cell_index - ((y - 1) * x) + 1].push_back(ghost);
      break;
    }
    case Placement::BOTTOM: {
      auto ghost = create_ghost_particle(particle_id, bottom_offset);
      cell_ghost_particles_map[cell_index + (y - 1) * x].push_back(ghost);
      cell_ghost_particles_map[cell_index + ((y - 1) * x) - 1].push_back(ghost);
      cell_ghost_particles_map[cell_index + ((y - 1) * x) + 1].push_back(ghost);
      break;
    }
    case Placement::BOTTOM_LEFT_CORNER: {
      // Corner ghost
      auto ghost_corner = create_ghost_particle(
          particle_id, {right_offset[0], bottom_offset[1], 0});
      cell_ghost_particles_map[cell_index + y * x - 1].push_back(ghost_corner);

      // Right side ghost
      auto ghost_right = create_ghost_particle(particle_id, right_offset);
      cell_ghost_particles_map[cell_index + x - 1].push_back(ghost_right);
      cell_ghost_particles_map[cell_index + x + x - 1].push_back(ghost_right);

      // Bottom side ghost
      auto ghost_bottom = create_ghost_particle(particle_id, bottom_offset);
      cell_ghost_particles_map[cell_index + ((y - 1) * x)].push_back(
          ghost_bottom);
      cell_ghost_particles_map[cell_index + ((y - 1) * x) + 1].push_back(
          ghost_bottom);
      break;
    }
    case Placement::TOP_RIGHT_CORNER: {
      // Corner ghost
      auto ghost_corner = create_ghost_particle(
          particle_id, {left_offset[0], top_offset[1], 0});
      cell_ghost_particles_map[cell_index - (y * x - 1)].push_back(
          ghost_corner);

      // Left side ghost
      auto ghost_left = create_ghost_particle(particle_id, left_offset);
      cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_left);
      cell_ghost_particles_map[cell_index - (x + x - 1)].push_back(ghost_left);

      // Top side ghost
      auto ghost_top = create_ghost_particle(particle_id, top_offset);
      cell_ghost_particles_map[cell_index - ((y - 1) * x)].push_back(ghost_top);
      cell_ghost_particles_map[cell_index - ((y - 1) * x) - 1].push_back(
          ghost_top);
      break;
    }
    case Placement::TOP_LEFT_CORNER: {
      // Corner ghost
      auto ghost_corner = create_ghost_particle(
          particle_id, {right_offset[0], top_offset[1], 0});
      cell_ghost_particles_map[cell_index - (y - 2) * x - 1].push_back(
          ghost_corner);

      // Right side ghost
      auto ghost_right = create_ghost_particle(particle_id, right_offset);
      cell_ghost_particles_map[cell_index + x - 1].push_back(ghost_right);
      cell_ghost_particles_map[cell_index - 1].push_back(ghost_right);

      // Top side ghost
      auto ghost_top = create_ghost_particle(particle_id, top_offset);
      cell_ghost_particles_map[cell_index - ((y - 1) * x)].push_back(ghost_top);
      cell_ghost_particles_map[cell_index - ((y - 1) * x) + 1].push_back(
          ghost_top);
      break;
    }
    case Placement::BOTTOM_RIGHT_CORNER: {
      // Corner ghost
      auto ghost_corner = create_ghost_particle(
          particle_id, {left_offset[0], bottom_offset[1], 0});
      cell_ghost_particles_map[cell_index + (y - 2) * x + 1].push_back(
          ghost_corner);

      // Left side ghost
      auto ghost_left = create_ghost_particle(particle_id, left_offset);
      cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_left);
      cell_ghost_particles_map[cell_index + 1].push_back(ghost_left);

      // Bottom side ghost
      auto ghost_bottom = create_ghost_particle(particle_id, bottom_offset);
      cell_ghost_particles_map[cell_index + ((y - 1) * x)].push_back(
          ghost_bottom);
      cell_ghost_particles_map[cell_index + ((y - 1) * x) - 1].push_back(
          ghost_bottom);
      break;
    }
    default:
      break;
    }
  } else if (domain_size_.size() == 3) {
    // Helper arrays for offsets
    const std::array<double, 3> left_offset = {domain_size_[0], 0, 0};
    const std::array<double, 3> right_offset = {-domain_size_[0], 0, 0};
    const std::array<double, 3> top_offset = {0, -domain_size_[1], 0};
    const std::array<double, 3> bottom_offset = {0, domain_size_[1], 0};
    const std::array<double, 3> front_offset = {0, 0, -domain_size_[2]};
    const std::array<double, 3> back_offset = {0, 0, domain_size_[2]};
    switch (placement) {
    case Placement::LEFT: {
      auto ghost = create_ghost_particle(particle_id, left_offset);
      // Neighbors of the halo ghost particle mirrored to the RIGHT of the
      // domain
      cell_ghost_particles_map[cell_index + x - 1].push_back(
          ghost); // Direct Left Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x].push_back(
          ghost); // Left-Top Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x].push_back(
          ghost); // Left-bottom Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x * y].push_back(
          ghost); // Left-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x * y].push_back(
          ghost); // Left-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x + x * y].push_back(
          ghost); // Left-Top-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x - x * y].push_back(
          ghost); // Left-Top-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x - x * y].push_back(
          ghost); // Left-Bottom-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x + x * y].push_back(
          ghost); // Left-Bottom-Front Neighbour
      break;
    }

    case Placement::RIGHT: {
      auto ghost = create_ghost_particle(particle_id, right_offset);
      // Neighbors of the halo ghost particle mirrored to the LEFT of the
      // domain
      cell_ghost_particles_map[cell_index - x + 1].push_back(ghost);
      // Right
      cell_ghost_particles_map[cell_index - x + 1 + x].push_back(
          ghost); // Right-Top
      cell_ghost_particles_map[cell_index - x + 1 - x].push_back(
          ghost); // Right-Bottom
      cell_ghost_particles_map[cell_index - x + 1 - x * y].push_back(
          ghost); // Right-Back
      cell_ghost_particles_map[cell_index - x + 1 + x * y].push_back(
          ghost); // Right-Front
      cell_ghost_particles_map[cell_index - x + 1 + x + x * y].push_back(
          ghost); // Right-Top-Front
      cell_ghost_particles_map[cell_index - x + 1 + x - x * y].push_back(
          ghost); // Right-Top-Back
      cell_ghost_particles_map[cell_index - x + 1 - x - x * y].push_back(
          ghost); // Right-Bottom-Back
      cell_ghost_particles_map[cell_index - x + 1 - x + x * y].push_back(
          ghost); // Right-Bottom-Front
      break;
    }

    case Placement::TOP: {
      auto ghost = create_ghost_particle(particle_id, top_offset);
      // Neighbors of the halo ghost particle mirrored to the TOP of the
      // domain
      cell_ghost_particles_map[cell_index - (y - 1) * x].push_back(
          ghost); // Top
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1].push_back(
          ghost); // Top-Left
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1].push_back(
          ghost); // Top-Right
      cell_ghost_particles_map[cell_index - (y - 1) * x - x * y].push_back(
          ghost); // Top-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x + x * y].push_back(
          ghost); // Top-Front
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1 - x * y].push_back(
          ghost); // Top-Right-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1 + x * y].push_back(
          ghost); // Top-Right-Front
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1 - x * y].push_back(
          ghost); // Top-Left-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1 + x * y].push_back(
          ghost); // Top-Left-Front
      break;
    }

    case Placement::BOTTOM: {
      auto ghost = create_ghost_particle(particle_id, bottom_offset);
      // Neighbors of the halo ghost particle mirrored to the BOTTOM of the
      // domain
      cell_ghost_particles_map[cell_index + (y - 1) * x].push_back(
          ghost); // Bottom
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1].push_back(
          ghost); // Bottom-Left
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1].push_back(
          ghost); // Bottom-Right
      cell_ghost_particles_map[cell_index + (y - 1) * x - x * y].push_back(
          ghost); // Bottom-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x + x * y].push_back(
          ghost); // Bottom-Front
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1 - x * y].push_back(
          ghost); // Bottom-Right-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1 + x * y].push_back(
          ghost); // Bottom-Right-Front
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1 - x * y].push_back(
          ghost); // Bottom-Left-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1 + x * y].push_back(
          ghost); // Bottom-Left-Front
      break;
    }

    case Placement::FRONT: {
      auto ghost = create_ghost_particle(particle_id, front_offset);
      // Neighbors of the halo ghost particle mirrored to the BACK of the
      // domain
      cell_ghost_particles_map[cell_index - (z - 1) * x * y].push_back(
          ghost); // Front
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1].push_back(
          ghost); // Front-Left
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1].push_back(
          ghost); // Front-Right
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + x].push_back(
          ghost); // Front-Top
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - x].push_back(
          ghost); // Front-Bottom
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1 + x].push_back(
          ghost); // Front-Right-Top
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1 - x].push_back(
          ghost); // Front-Right-Bottom
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1 + x].push_back(
          ghost); // Front-Left-Top
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1 - x].push_back(
          ghost); // Front-Left-Bottom
      break;
    }

    case Placement::BACK: {
      auto ghost = create_ghost_particle(particle_id, back_offset);
      // Neighbors of the halo ghost particle mirrored to the FRONT of the
      // domain
      cell_ghost_particles_map[cell_index + (z - 1) * x * y].push_back(
          ghost); // Back
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1].push_back(
          ghost); // Back-Left
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1].push_back(
          ghost); // Back-Right
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + x].push_back(
          ghost); // Back-Top
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - x].push_back(
          ghost); // Back-Bottom
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1 + x].push_back(
          ghost); // Back-Right-Top
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1 - x].push_back(
          ghost); // Back-Right-Bottom
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1 + x].push_back(
          ghost); // Back-Left-Top
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1 - x].push_back(
          ghost); // Back-Left-Bottom
      break;
    }

    case Placement::LEFT_FRONT_EDGE: {
      auto ghost_1 = create_ghost_particle(particle_id, left_offset);
      // map 6 neighbours
      cell_ghost_particles_map[cell_index + x - 1].push_back(ghost_1); // Left
      cell_ghost_particles_map[cell_index + x + x - 1].push_back(
          ghost_1); // Left-Top
      cell_ghost_particles_map[cell_index - 1].push_back(
          ghost_1); // Left-Bottom
      cell_ghost_particles_map[cell_index + x - 1 - x * y].push_back(
          ghost_1); // Left-Back
      cell_ghost_particles_map[cell_index + x + x - 1 - x * y].push_back(
          ghost_1); // Left-Top-Back
      cell_ghost_particles_map[cell_index - 1 - x * y].push_back(
          ghost_1); // Left-Bottom-Back

      auto ghost_2 = create_ghost_particle(particle_id, front_offset);
      // map 6 neighbours
      cell_ghost_particles_map[cell_index - (z - 1) * x * y].push_back(
          ghost_2); // Front
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1].push_back(
          ghost_2); // Front-Right
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + x].push_back(
          ghost_2); // Front-Top
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - x].push_back(
          ghost_2); // Front-Bottom
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1 + x].push_back(
          ghost_2); // Front-Right-Top
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1 - x].push_back(
          ghost_2); // Front-Right-Bottom

      auto ghost_3 = create_ghost_particle(
          particle_id, {left_offset[0], 0, front_offset[2]});
      // map 3 neighbours
      cell_ghost_particles_map[cell_index + x - 1 - (z - 1) * x * y].push_back(
          ghost_3); // Left-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x - (z - 1) * x * y]
          .push_back(ghost_3); // Left-Top-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x - (z - 1) * x * y]
          .push_back(ghost_3); // Left-Bottom-Front Neighbour
      break;
    }

    case Placement::LEFT_BACK_EDGE: {
      auto ghost_1 = create_ghost_particle(particle_id, left_offset);
      cell_ghost_particles_map[cell_index + x - 1].push_back(
          ghost_1); // Direct Left Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x].push_back(
          ghost_1); // Left-Top Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x].push_back(
          ghost_1); // Left-bottom Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x * y].push_back(
          ghost_1); // Left-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x + x * y].push_back(
          ghost_1); // Left-Top-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x + x * y].push_back(
          ghost_1); // Left-Bottom-Front Neighbour

      auto ghost_2 = create_ghost_particle(particle_id, back_offset);
      cell_ghost_particles_map[cell_index + (z - 1) * x * y].push_back(
          ghost_2); // Back
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1].push_back(
          ghost_2); // Back-Right
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + x].push_back(
          ghost_2); // Back-Top
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - x].push_back(
          ghost_2); // Back-Bottom
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1 + x].push_back(
          ghost_2); // Back-Right-Top
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1 - x].push_back(
          ghost_2); // Back-Right-Bottom

      auto ghost_3 = create_ghost_particle(particle_id,
                                           {left_offset[0], 0, back_offset[2]});
      cell_ghost_particles_map[cell_index + x - 1 + (z - 1) * x * y].push_back(
          ghost_3); // Left-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x + (z - 1) * x * y]
          .push_back(ghost_3); // Left-Top-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x + (z - 1) * x * y]
          .push_back(ghost_3); // Left-Bottom-Back Neighbour
      break;
    }

    case Placement::RIGHT_FRONT_EDGE: {
      auto ghost_1 = create_ghost_particle(particle_id, right_offset);
      cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_1); // Right
      cell_ghost_particles_map[cell_index - x + 1 + x].push_back(
          ghost_1); // Right-Top
      cell_ghost_particles_map[cell_index - x + 1 - x].push_back(
          ghost_1); // Right-Bottom
      cell_ghost_particles_map[cell_index - x + 1 - x * y].push_back(
          ghost_1); // Right-Back
      cell_ghost_particles_map[cell_index - x + 1 + x - x * y].push_back(
          ghost_1); // Right-Top-Back
      cell_ghost_particles_map[cell_index - x + 1 - x - x * y].push_back(
          ghost_1); // Right-Bottom-Back

      auto ghost_2 = create_ghost_particle(particle_id, front_offset);
      cell_ghost_particles_map[cell_index - (z - 1) * x * y].push_back(
          ghost_2); // Front
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1].push_back(
          ghost_2); // Front-Left
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + x].push_back(
          ghost_2); // Front-Top
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - x].push_back(
          ghost_2); // Front-Bottom
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1 + x].push_back(
          ghost_2); // Front-Left-Top
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1 - x].push_back(
          ghost_2); // Front-Left-Bottom

      auto ghost_3 = create_ghost_particle(
          particle_id, {right_offset[0], 0, front_offset[2]});
      cell_ghost_particles_map[cell_index - x + 1 - (z - 1) * x * y].push_back(
          ghost_3); // Right-Front
      cell_ghost_particles_map[cell_index - x + 1 + x - (z - 1) * x * y]
          .push_back(ghost_3); // Right-Top-Front
      cell_ghost_particles_map[cell_index - x + 1 - x - (z - 1) * x * y]
          .push_back(ghost_3); // Right-Bottom-Front
      break;
    }

    case Placement::RIGHT_BACK_EDGE: {
      auto ghost_1 = create_ghost_particle(particle_id, right_offset);
      cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_1); // Right
      cell_ghost_particles_map[cell_index - x + 1 + x].push_back(
          ghost_1); // Right-Top
      cell_ghost_particles_map[cell_index - x + 1 - x].push_back(
          ghost_1); // Right-Bottom
      cell_ghost_particles_map[cell_index - x + 1 + x * y].push_back(
          ghost_1); // Right-Front
      cell_ghost_particles_map[cell_index - x + 1 + x + x * y].push_back(
          ghost_1); // Right-Top-Front
      cell_ghost_particles_map[cell_index - x + 1 - x + x * y].push_back(
          ghost_1); // Right-Bottom-Front

      auto ghost_2 = create_ghost_particle(particle_id, back_offset);
      cell_ghost_particles_map[cell_index + (z - 1) * x * y].push_back(
          ghost_2); // Back
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1].push_back(
          ghost_2); // Back-Left
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + x].push_back(
          ghost_2); // Back-Top
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - x].push_back(
          ghost_2); // Back-Bottom
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1 + x].push_back(
          ghost_2); // Back-Left-Top
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1 - x].push_back(
          ghost_2); // Back-Left-Bottom

      auto ghost_3 = create_ghost_particle(
          particle_id, {right_offset[0], 0, back_offset[2]});
      cell_ghost_particles_map[cell_index - x + 1 + (z - 1) * x * y].push_back(
          ghost_3); // Right-Back
      cell_ghost_particles_map[cell_index - x + 1 + x + (z - 1) * x * y]
          .push_back(ghost_3); // Right-Top-Back
      cell_ghost_particles_map[cell_index - x + 1 - x + (z - 1) * x * y]
          .push_back(ghost_3); // Right-Bottom-Back
      break;
    }

    case Placement::TOP_LEFT_EDGE: {
      auto ghost_1 = create_ghost_particle(particle_id, left_offset);
      cell_ghost_particles_map[cell_index + x - 1].push_back(
          ghost_1); // Direct Left Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x].push_back(
          ghost_1); // Left-bottom Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x * y].push_back(
          ghost_1); // Left-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x * y].push_back(
          ghost_1); // Left-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x - x * y].push_back(
          ghost_1); // Left-Bottom-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x + x * y].push_back(
          ghost_1); // Left-Bottom-Front Neighbour

      auto ghost_2 = create_ghost_particle(particle_id, top_offset);
      cell_ghost_particles_map[cell_index - (y - 1) * x].push_back(
          ghost_2); // Top
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1].push_back(
          ghost_2); // Top-Right
      cell_ghost_particles_map[cell_index - (y - 1) * x - x * y].push_back(
          ghost_2); // Top-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x + x * y].push_back(
          ghost_2); // Top-Front
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1 - x * y].push_back(
          ghost_2); // Top-Right-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1 + x * y].push_back(
          ghost_2); // Top-Right-Front

      auto ghost_3 = create_ghost_particle(particle_id,
                                           {left_offset[0], top_offset[1], 0});
      cell_ghost_particles_map[cell_index + x - 1 - (y - 1) * x].push_back(
          ghost_3); // Left-Top Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - (y - 1) * x + x * y]
          .push_back(ghost_3); // Left-Top-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - (y - 1) * x - x * y]
          .push_back(ghost_3); // Left-Top-Back Neighbour
      break;
    }

    case Placement::TOP_RIGHT_EDGE: {
      auto ghost_1 = create_ghost_particle(particle_id, right_offset);
      cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_1); // Right
      cell_ghost_particles_map[cell_index - x + 1 - x].push_back(
          ghost_1); // Right-Bottom
      cell_ghost_particles_map[cell_index - x + 1 - x * y].push_back(
          ghost_1); // Right-Back
      cell_ghost_particles_map[cell_index - x + 1 + x * y].push_back(
          ghost_1); // Right-Front
      cell_ghost_particles_map[cell_index - x + 1 - x - x * y].push_back(
          ghost_1); // Right-Bottom-Back
      cell_ghost_particles_map[cell_index - x + 1 - x + x * y].push_back(
          ghost_1); // Right-Bottom-Front

      auto ghost_2 = create_ghost_particle(particle_id, top_offset);
      cell_ghost_particles_map[cell_index - (y - 1) * x].push_back(
          ghost_2); // Top
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1].push_back(
          ghost_2); // Top-Left
      cell_ghost_particles_map[cell_index - (y - 1) * x - x * y].push_back(
          ghost_2); // Top-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x + x * y].push_back(
          ghost_2); // Top-Front
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1 - x * y].push_back(
          ghost_2); // Top-Left-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1 + x * y].push_back(
          ghost_2); // Top-Left-Front
      auto ghost_3 = create_ghost_particle(particle_id,
                                           {right_offset[0], top_offset[1], 0});
      cell_ghost_particles_map[cell_index - x + 1 - (y - 1) * x].push_back(
          ghost_3); // Right-Top
      cell_ghost_particles_map[cell_index - x + 1 - (y - 1) * x + x * y]
          .push_back(ghost_3); // Right-Top-Front
      cell_ghost_particles_map[cell_index - x + 1 - (y - 1) * x - x * y]
          .push_back(ghost_3); // Right-Top-Back
      break;
    }

    case Placement::BOTTOM_RIGHT_EDGE: {
      auto ghost_1 = create_ghost_particle(particle_id, right_offset);
      cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_1); // Right
      cell_ghost_particles_map[cell_index - x + 1 + x].push_back(
          ghost_1); // Right-Top
      cell_ghost_particles_map[cell_index - x + 1 - x * y].push_back(
          ghost_1); // Right-Back
      cell_ghost_particles_map[cell_index - x + 1 + x * y].push_back(
          ghost_1); // Right-Front
      cell_ghost_particles_map[cell_index - x + 1 + x + x * y].push_back(
          ghost_1); // Right-Top-Front
      cell_ghost_particles_map[cell_index - x + 1 + x - x * y].push_back(
          ghost_1); // Right-Top-Back

      auto ghost_2 = create_ghost_particle(particle_id, bottom_offset);
      cell_ghost_particles_map[cell_index + (y - 1) * x].push_back(
          ghost_2); // Bottom
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1].push_back(
          ghost_2); // Bottom-Left
      cell_ghost_particles_map[cell_index + (y - 1) * x - x * y].push_back(
          ghost_2); // Bottom-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x + x * y].push_back(
          ghost_2); // Bottom-Front
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1 - x * y].push_back(
          ghost_2); // Bottom-Left-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1 + x * y].push_back(
          ghost_2); // Bottom-Left-Front
      auto ghost_3 = create_ghost_particle(
          particle_id, {right_offset[0], bottom_offset[1], 0});
      cell_ghost_particles_map[cell_index - x + 1 + (y - 1) * x].push_back(
          ghost_3); // Right-Bottom
      cell_ghost_particles_map[cell_index - x + 1 + (y - 1) * x - x * y]
          .push_back(ghost_3); // Right-Bottom-Back
      cell_ghost_particles_map[cell_index - x + 1 + (y - 1) * x + x * y]
          .push_back(ghost_3); // Right-Bottom-Front
      break;
    }

    case Placement::BOTTOM_LEFT_EDGE: {
      auto ghost_1 = create_ghost_particle(particle_id, left_offset);
      cell_ghost_particles_map[cell_index + x - 1].push_back(
          ghost_1); // Direct Left Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x].push_back(
          ghost_1); // Left-Top Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x * y].push_back(
          ghost_1); // Left-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x * y].push_back(
          ghost_1); // Left-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x + x * y].push_back(
          ghost_1); // Left-Top-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x - x * y].push_back(
          ghost_1); // Left-Top-Back Neighbour
      auto ghost_2 = create_ghost_particle(particle_id, bottom_offset);
      cell_ghost_particles_map[cell_index + (y - 1) * x].push_back(
          ghost_2); // Bottom
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1].push_back(
          ghost_2); // Bottom-Right
      cell_ghost_particles_map[cell_index + (y - 1) * x - x * y].push_back(
          ghost_2); // Bottom-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x + x * y].push_back(
          ghost_2); // Bottom-Front
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1 - x * y].push_back(
          ghost_2); // Bottom-Right-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1 + x * y].push_back(
          ghost_2); // Bottom-Right-Front

      auto ghost_3 = create_ghost_particle(
          particle_id, {left_offset[0], bottom_offset[1], 0});
      cell_ghost_particles_map[cell_index + x - 1 + (y - 1) * x].push_back(
          ghost_3); // Left-bottom Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + (y - 1) * x - x * y]
          .push_back(ghost_3); // Left-Bottom-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + (y - 1) * x + x * y]
          .push_back(ghost_3); // Left-Bottom-Front Neighbour
      break;
    }

    case Placement::TOP_FRONT_EDGE: {
      auto ghost_1 = create_ghost_particle(particle_id, top_offset);
      cell_ghost_particles_map[cell_index - (y - 1) * x].push_back(
          ghost_1); // Top
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1].push_back(
          ghost_1); // Top-Left
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1].push_back(
          ghost_1); // Top-Right
      cell_ghost_particles_map[cell_index - (y - 1) * x - x * y].push_back(
          ghost_1); // Top-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1 - x * y].push_back(
          ghost_1); // Top-Right-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1 - x * y].push_back(
          ghost_1); // Top-Left-Back

      auto ghost_2 = create_ghost_particle(particle_id, front_offset);
      cell_ghost_particles_map[cell_index - (z - 1) * x * y].push_back(
          ghost_2); // Front
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1].push_back(
          ghost_2); // Front-Left
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1].push_back(
          ghost_2); // Front-Right
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - x].push_back(
          ghost_2); // Front-Bottom
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1 - x].push_back(
          ghost_2); // Front-Right-Bottom
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1 - x].push_back(
          ghost_2); // Front-Left-Bottom

      auto ghost_3 = create_ghost_particle(particle_id,
                                           {0, top_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index - (y - 1) * x - (z - 1) * x * y]
          .push_back(ghost_3); // Top-Front
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1 - (z - 1) * x * y]
          .push_back(ghost_3); // Top-Right-Front
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1 - (z - 1) * x * y]
          .push_back(ghost_3); // Top-Left-Front
      break;
    }

    case Placement::BOTTOM_FRONT_EDGE: {
      auto ghost_1 = create_ghost_particle(particle_id, bottom_offset);
      cell_ghost_particles_map[cell_index + (y - 1) * x].push_back(
          ghost_1); // Bottom
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1].push_back(
          ghost_1); // Bottom-Left
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1].push_back(
          ghost_1); // Bottom-Right
      cell_ghost_particles_map[cell_index + (y - 1) * x - x * y].push_back(
          ghost_1); // Bottom-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1 - x * y].push_back(
          ghost_1); // Bottom-Right-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1 - x * y].push_back(
          ghost_1); // Bottom-Left-Back
      auto ghost_2 = create_ghost_particle(particle_id, front_offset);
      cell_ghost_particles_map[cell_index - (z - 1) * x * y].push_back(
          ghost_2); // Front
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1].push_back(
          ghost_2); // Front-Left
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1].push_back(
          ghost_2); // Front-Right
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + x].push_back(
          ghost_2); // Front-Top
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1 + x].push_back(
          ghost_2); // Front-Right-Top
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1 + x].push_back(
          ghost_2); // Front-Left-Top
      auto ghost_3 = create_ghost_particle(
          particle_id, {0, bottom_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + (y - 1) * x - (z - 1) * x * y]
          .push_back(ghost_3); // Bottom-Front
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1 - (z - 1) * x * y]
          .push_back(ghost_3); // Bottom-Right-Front
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1 - (z - 1) * x * y]
          .push_back(ghost_3); // Bottom-Left-Front
      break;
    }

    case Placement::BOTTOM_BACK_EDGE: {
      auto ghost_1 = create_ghost_particle(particle_id, bottom_offset);
      cell_ghost_particles_map[cell_index + (y - 1) * x].push_back(
          ghost_1); // Bottom
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1].push_back(
          ghost_1); // Bottom-Left
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1].push_back(
          ghost_1); // Bottom-Right
      cell_ghost_particles_map[cell_index + (y - 1) * x + x * y].push_back(
          ghost_1); // Bottom-Front
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1 + x * y].push_back(
          ghost_1); // Bottom-Right-Front
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1 + x * y].push_back(
          ghost_1); // Bottom-Left-Front

      auto ghost_2 = create_ghost_particle(particle_id, back_offset);
      cell_ghost_particles_map[cell_index + (z - 1) * x * y].push_back(
          ghost_2); // Back
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1].push_back(
          ghost_2); // Back-Left
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1].push_back(
          ghost_2); // Back-Right
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + x].push_back(
          ghost_2); // Back-Top
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1 + x].push_back(
          ghost_2); // Back-Right-Top
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1 + x].push_back(
          ghost_2); // Back-Left-Top

      auto ghost_3 = create_ghost_particle(
          particle_id, {0, bottom_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index + (y - 1) * x + (z - 1) * x * y]
          .push_back(ghost_3); // Bottom-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1 + (z - 1) * x * y]
          .push_back(ghost_3); // Bottom-Right-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1 + (z - 1) * x * y]
          .push_back(ghost_3); // Bottom-Left-Back
      break;
    }

    case Placement::TOP_BACK_EDGE: {
      auto ghost_1 = create_ghost_particle(particle_id, top_offset);
      cell_ghost_particles_map[cell_index - (y - 1) * x].push_back(
          ghost_1); // Top
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1].push_back(
          ghost_1); // Top-Left
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1].push_back(
          ghost_1); // Top-Right
      cell_ghost_particles_map[cell_index - (y - 1) * x + x * y].push_back(
          ghost_1); // Top-Front
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1 + x * y].push_back(
          ghost_1); // Top-Right-Front
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1 + x * y].push_back(
          ghost_1); // Top-Left-Front

      auto ghost_2 = create_ghost_particle(particle_id, back_offset);
      cell_ghost_particles_map[cell_index + (z - 1) * x * y].push_back(
          ghost_2); // Back
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1].push_back(
          ghost_2); // Back-Left
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1].push_back(
          ghost_2); // Back-Right
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - x].push_back(
          ghost_2); // Back-Bottom
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1 - x].push_back(
          ghost_2); // Back-Right-Bottom
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1 - x].push_back(
          ghost_2); // Back-Left-Bottom
      auto ghost_3 = create_ghost_particle(particle_id,
                                           {0, top_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index - (y - 1) * x + (z - 1) * x * y]
          .push_back(ghost_3); // Top-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1 + (z - 1) * x * y]
          .push_back(ghost_3); // Top-Right-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1 + (z - 1) * x * y]
          .push_back(ghost_3); // Top-Left-Back
      break;
    }

    case Placement::TOP_FRONT_RIGHT_CORNER: {
      auto ghost_1 = create_ghost_particle(particle_id, top_offset);
      cell_ghost_particles_map[cell_index - (y - 1) * x].push_back(
          ghost_1); // Top
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1].push_back(
          ghost_1); // Top-Left
      cell_ghost_particles_map[cell_index - (y - 1) * x - x * y].push_back(
          ghost_1); // Top-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1 - x * y].push_back(
          ghost_1); // Top-Left-Back
      auto ghost_2 = create_ghost_particle(particle_id, front_offset);
      cell_ghost_particles_map[cell_index - (z - 1) * x * y].push_back(
          ghost_2); // Front
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1].push_back(
          ghost_2); // Front-Left
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - x].push_back(
          ghost_2); // Front-Bottom
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1 - x].push_back(
          ghost_2); // Front-Left-Bottom
      auto ghost_3 = create_ghost_particle(particle_id, right_offset);
      cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_3); // Right
      cell_ghost_particles_map[cell_index - x + 1 - x].push_back(
          ghost_3); // Right-Bottom
      cell_ghost_particles_map[cell_index - x + 1 - x * y].push_back(
          ghost_3); // Right-Back
      cell_ghost_particles_map[cell_index - x + 1 - x - x * y].push_back(
          ghost_3); // Right-Bottom-Back
      auto ghost_4 = create_ghost_particle(particle_id,
                                           {right_offset[0], top_offset[1], 0});
      cell_ghost_particles_map[cell_index - x + 1 - (y - 1) * x].push_back(
          ghost_4); // Right-Top
      cell_ghost_particles_map[cell_index - x + 1 - (y - 1) * x - x * y]
          .push_back(ghost_4); // Right-Top-Back
      auto ghost_5 = create_ghost_particle(
          particle_id, {right_offset[0], 0, front_offset[2]});
      cell_ghost_particles_map[cell_index - x + 1 - (z - 1) * x * y].push_back(
          ghost_5); // Right-Front
      cell_ghost_particles_map[cell_index - x + 1 - x - (z - 1) * x * y]
          .push_back(ghost_5); // Right-Bottom-Front
      auto ghost_6 = create_ghost_particle(particle_id,
                                           {0, top_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index - (y - 1) * x - (z - 1) * x * y]
          .push_back(ghost_6); // Top-Front
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1 - (z - 1) * x * y]
          .push_back(ghost_6); // Top-Left-Front
      auto ghost_7 = create_ghost_particle(
          particle_id, {right_offset[0], top_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index - x + 1 - (y - 1) * x + x * y]
          .push_back(ghost_7); // Right-Top-Front
      break;
    }

    case Placement::TOP_FRONT_LEFT_CORNER: {
      auto ghost_1 = create_ghost_particle(particle_id, top_offset);
      cell_ghost_particles_map[cell_index - (y - 1) * x].push_back(
          ghost_1); // Top
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1].push_back(
          ghost_1); // Top-Right
      cell_ghost_particles_map[cell_index - (y - 1) * x - x * y].push_back(
          ghost_1); // Top-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1 - x * y].push_back(
          ghost_1); // Top-Right-Back
      auto ghost_2 = create_ghost_particle(particle_id, front_offset);
      cell_ghost_particles_map[cell_index - (z - 1) * x * y].push_back(
          ghost_2); // Front
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1].push_back(
          ghost_2); // Front-Right
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - x].push_back(
          ghost_2); // Front-Bottom
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1 - x].push_back(
          ghost_2); // Front-Right-Bottom
      auto ghost_3 = create_ghost_particle(particle_id, left_offset);
      cell_ghost_particles_map[cell_index + x - 1].push_back(
          ghost_3); // Direct Left Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x].push_back(
          ghost_3); // Left-bottom Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x * y].push_back(
          ghost_3); // Left-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x - x * y].push_back(
          ghost_3); // Left-Bottom-Back Neighbour
      auto ghost_4 = create_ghost_particle(particle_id,
                                           {left_offset[0], top_offset[1], 0});
      cell_ghost_particles_map[cell_index + x - 1 - (y - 1) * x].push_back(
          ghost_4); // Left-Top Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - (y - 1) * x - x * y]
          .push_back(ghost_4); // Left-Top-Back Neighbour
      auto ghost_5 = create_ghost_particle(
          particle_id, {left_offset[0], 0, front_offset[2]});
      cell_ghost_particles_map[cell_index + x - 1 - (z - 1) * x * y].push_back(
          ghost_5); // Left-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x - (z - 1) * x * y]
          .push_back(ghost_5); // Left-Bottom-Front Neighbour
      auto ghost_6 = create_ghost_particle(particle_id,
                                           {0, top_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index - (y - 1) * x - (z - 1) * x * y]
          .push_back(ghost_6); // Top-Front
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1 - (z - 1) * x * y]
          .push_back(ghost_6); // Top-Right-Front
      auto ghost_7 = create_ghost_particle(
          particle_id, {left_offset[0], top_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index - (y - 1) * x + x - 1 -
                               (z - 1) * x * y]
          .push_back(ghost_7); // Top-Left-Front
      break;
    }

    case Placement::TOP_BACK_LEFT_CORNER: {
      auto ghost_1 = create_ghost_particle(particle_id, top_offset);
      cell_ghost_particles_map[cell_index - (y - 1) * x].push_back(
          ghost_1); // Top
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1].push_back(
          ghost_1); // Top-Right
      cell_ghost_particles_map[cell_index - (y - 1) * x + x * y].push_back(
          ghost_1); // Top-Front
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1 + x * y].push_back(
          ghost_1); // Top-Right-Front
      auto ghost_2 = create_ghost_particle(particle_id, back_offset);
      cell_ghost_particles_map[cell_index + (z - 1) * x * y].push_back(
          ghost_2); // Back
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1].push_back(
          ghost_2); // Back-Right
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - x].push_back(
          ghost_2); // Back-Bottom
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1 - x].push_back(
          ghost_2); // Back-Right-Bottom
      auto ghost_3 = create_ghost_particle(particle_id, left_offset);
      cell_ghost_particles_map[cell_index + x - 1].push_back(
          ghost_3); // Direct Left Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x].push_back(
          ghost_3); // Left-bottom Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x * y].push_back(
          ghost_3); // Left-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x + x * y].push_back(
          ghost_3); // Left-Bottom-Front Neighbour
      auto ghost_4 = create_ghost_particle(particle_id,
                                           {left_offset[0], top_offset[1], 0});
      cell_ghost_particles_map[cell_index + x - 1 - (y - 1) * x].push_back(
          ghost_4); // Left-Top Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - (y - 1) * x + x * y]
          .push_back(ghost_4); // Left-Top-Front Neighbour
      auto ghost_5 = create_ghost_particle(particle_id,
                                           {left_offset[0], 0, back_offset[2]});
      cell_ghost_particles_map[cell_index + x - 1 + (z - 1) * x * y].push_back(
          ghost_5); // Left-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x + (z - 1) * x * y]
          .push_back(ghost_5); // Left-Bottom-Back Neighbour

      auto ghost_6 = create_ghost_particle(particle_id,
                                           {0, top_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index - (y - 1) * x + (z - 1) * x * y]
          .push_back(ghost_6); // Top-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x + 1 + (z - 1) * x * y]
          .push_back(ghost_6); // Top-Right-Back
      auto ghost_7 = create_ghost_particle(
          particle_id, {left_offset[0], top_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index - (y - 1) * x + x - 1 +
                               (z - 1) * x * y]
          .push_back(ghost_7); // Top-Left-Back
      break;
    }

    case Placement::TOP_BACK_RIGHT_CORNER: {
      auto ghost_1 = create_ghost_particle(particle_id, top_offset);
      cell_ghost_particles_map[cell_index - (y - 1) * x].push_back(
          ghost_1); // Top
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1].push_back(
          ghost_1); // Top-Left
      cell_ghost_particles_map[cell_index - (y - 1) * x + x * y].push_back(
          ghost_1); // Top-Front
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1 + x * y].push_back(
          ghost_1); // Top-Left-Front

      auto ghost_2 = create_ghost_particle(particle_id, back_offset);
      cell_ghost_particles_map[cell_index + (z - 1) * x * y].push_back(
          ghost_2); // Back
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1].push_back(
          ghost_2); // Back-Left
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - x].push_back(
          ghost_2); // Back-Bottom
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1 - x].push_back(
          ghost_2); // Back-Left-Bottom

      auto ghost_3 = create_ghost_particle(particle_id, right_offset);
      cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_3); // Right
      cell_ghost_particles_map[cell_index - x + 1 - x].push_back(
          ghost_3); // Right-Bottom
      cell_ghost_particles_map[cell_index - x + 1 + x * y].push_back(
          ghost_3); // Right-Front
      cell_ghost_particles_map[cell_index - x + 1 - x + x * y].push_back(
          ghost_3); // Right-Bottom-Front

      auto ghost_4 = create_ghost_particle(particle_id,
                                           {right_offset[0], top_offset[1], 0});
      cell_ghost_particles_map[cell_index - x + 1 - (y - 1) * x].push_back(
          ghost_4); // Right-Top
      cell_ghost_particles_map[cell_index - x + 1 - (y - 1) * x + x * y]
          .push_back(ghost_4); // Right-Top-Front

      auto ghost_5 = create_ghost_particle(
          particle_id, {right_offset[0], 0, back_offset[2]});
      cell_ghost_particles_map[cell_index - x + 1 + (z - 1) * x * y].push_back(
          ghost_5); // Right-Back
      cell_ghost_particles_map[cell_index - x + 1 - x + (z - 1) * x * y]
          .push_back(ghost_5); // Right-Bottom-Back

      auto ghost_6 = create_ghost_particle(particle_id,
                                           {0, top_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index - (y - 1) * x + (z - 1) * x * y]
          .push_back(ghost_6); // Top-Back
      cell_ghost_particles_map[cell_index - (y - 1) * x - 1 + (z - 1) * x * y]
          .push_back(ghost_6); // Top-Left-Back

      auto ghost_7 = create_ghost_particle(
          particle_id, {right_offset[0], top_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index - (y - 1) * x - x + 1 +
                               (z - 1) * x * y]
          .push_back(ghost_7); // Top-Right-Back
      break;
    }

    case Placement::BOTTOM_FRONT_RIGHT_CORNER: {
      auto ghost_1 = create_ghost_particle(particle_id, bottom_offset);
      cell_ghost_particles_map[cell_index + (y - 1) * x].push_back(
          ghost_1); // Bottom
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1].push_back(
          ghost_1); // Bottom-Left
      cell_ghost_particles_map[cell_index + (y - 1) * x - x * y].push_back(
          ghost_1); // Bottom-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1 - x * y].push_back(
          ghost_1); // Bottom-Left-Back

      auto ghost_2 = create_ghost_particle(particle_id, front_offset);
      cell_ghost_particles_map[cell_index - (z - 1) * x * y].push_back(
          ghost_2); // Front
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1].push_back(
          ghost_2); // Front-Left
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + x].push_back(
          ghost_2); // Front-Top
      cell_ghost_particles_map[cell_index - (z - 1) * x * y - 1 + x].push_back(
          ghost_2); // Front-Left-Top

      auto ghost_3 = create_ghost_particle(particle_id, right_offset);
      cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_3); // Right
      cell_ghost_particles_map[cell_index - x + 1 + x].push_back(
          ghost_3); // Right-Top
      cell_ghost_particles_map[cell_index - x + 1 - x * y].push_back(
          ghost_3); // Right-Back
      cell_ghost_particles_map[cell_index - x + 1 + x - x * y].push_back(
          ghost_3); // Right-Top-Back

      auto ghost_4 = create_ghost_particle(
          particle_id, {right_offset[0], bottom_offset[1], 0});
      cell_ghost_particles_map[cell_index - x + 1 + (y - 1) * x].push_back(
          ghost_4); // Right-Bottom
      cell_ghost_particles_map[cell_index - x + 1 + (y - 1) * x - x * y]
          .push_back(ghost_4); // Right-Bottom-Back

      auto ghost_5 = create_ghost_particle(
          particle_id, {right_offset[0], 0, front_offset[2]});
      cell_ghost_particles_map[cell_index - x + 1 - (z - 1) * x * y].push_back(
          ghost_5); // Right-Front
      cell_ghost_particles_map[cell_index - x + 1 + x - (z - 1) * x * y]
          .push_back(ghost_5); // Right-Top-Front

      auto ghost_6 = create_ghost_particle(
          particle_id, {0, bottom_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + (y - 1) * x - (z - 1) * x * y]
          .push_back(ghost_6); // Bottom-Front
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1 - (z - 1) * x * y]
          .push_back(ghost_6); // Bottom-Left-Front

      auto ghost_7 = create_ghost_particle(
          particle_id, {right_offset[0], bottom_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + (y - 1) * x - x + 1 -
                               (z - 1) * x * y]
          .push_back(ghost_7); // Bottom-Right-Front
      break;
    }

    case Placement::BOTTOM_FRONT_LEFT_CORNER: {
      auto ghost_1 = create_ghost_particle(particle_id, bottom_offset);
      cell_ghost_particles_map[cell_index + (y - 1) * x].push_back(
          ghost_1); // Bottom
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1].push_back(
          ghost_1); // Bottom-Right
      cell_ghost_particles_map[cell_index + (y - 1) * x - x * y].push_back(
          ghost_1); // Bottom-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1 - x * y].push_back(
          ghost_1); // Bottom-Right-Back

      auto ghost_2 = create_ghost_particle(particle_id, front_offset);
      cell_ghost_particles_map[cell_index - (z - 1) * x * y].push_back(
          ghost_2); // Front
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1].push_back(
          ghost_2); // Front-Right
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + x].push_back(
          ghost_2); // Front-Top
      cell_ghost_particles_map[cell_index - (z - 1) * x * y + 1 + x].push_back(
          ghost_2); // Front-Right-Top

      auto ghost_3 = create_ghost_particle(particle_id, left_offset);
      cell_ghost_particles_map[cell_index + x - 1].push_back(
          ghost_3); // Direct Left Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x].push_back(
          ghost_3); // Left-Top Neighbour
      cell_ghost_particles_map[cell_index + x - 1 - x * y].push_back(
          ghost_3); // Left-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x - x * y].push_back(
          ghost_3); // Left-Top-Back Neighbour

      auto ghost_4 = create_ghost_particle(
          particle_id, {left_offset[0], bottom_offset[1], 0});
      cell_ghost_particles_map[cell_index + x - 1 + (y - 1) * x].push_back(
          ghost_4); // Left-bottom Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + (y - 1) * x - x * y]
          .push_back(ghost_4); // Left-Bottom-Back Neighbour

      auto ghost_5 = create_ghost_particle(
          particle_id, {left_offset[0], 0, front_offset[2]});
      cell_ghost_particles_map[cell_index + x - 1 - (z - 1) * x * y].push_back(
          ghost_5); // Left-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x - (z - 1) * x * y]
          .push_back(ghost_5); // Left-Top-Front Neighbour

      auto ghost_6 = create_ghost_particle(
          particle_id, {0, bottom_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + (y - 1) * x - (z - 1) * x * y]
          .push_back(ghost_6); // Bottom-Front
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1 - (z - 1) * x * y]
          .push_back(ghost_6); // Bottom-Right-Front

      auto ghost_7 = create_ghost_particle(
          particle_id, {left_offset[0], bottom_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + (y - 1) * x + x - 1 -
                               (z - 1) * x * y]
          .push_back(ghost_7); // Bottom-Left-Front
      break;
    }

    case Placement::BOTTOM_BACK_LEFT_CORNER: {
      auto ghost_1 = create_ghost_particle(particle_id, bottom_offset);
      cell_ghost_particles_map[cell_index + (y - 1) * x].push_back(
          ghost_1); // Bottom
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1].push_back(
          ghost_1); // Bottom-Right
      cell_ghost_particles_map[cell_index + (y - 1) * x + x * y].push_back(
          ghost_1); // Bottom-Front
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1 + x * y].push_back(
          ghost_1); // Bottom-Right-Front

      auto ghost_2 = create_ghost_particle(particle_id, back_offset);
      cell_ghost_particles_map[cell_index + (z - 1) * x * y].push_back(
          ghost_2); // Back
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1].push_back(
          ghost_2); // Back-Right
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + x].push_back(
          ghost_2); // Back-Top
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + 1 + x].push_back(
          ghost_2); // Back-Right-Top

      auto ghost_3 = create_ghost_particle(particle_id, left_offset);
      cell_ghost_particles_map[cell_index + x - 1].push_back(
          ghost_3); // Direct Left Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x].push_back(
          ghost_3); // Left-Top Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x * y].push_back(
          ghost_3); // Left-Front Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x + x * y].push_back(
          ghost_3); // Left-Top-Front Neighbour

      auto ghost_4 = create_ghost_particle(
          particle_id, {left_offset[0], bottom_offset[1], 0});
      cell_ghost_particles_map[cell_index + x - 1 + (y - 1) * x].push_back(
          ghost_4); // Left-bottom Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + (y - 1) * x + x * y]
          .push_back(ghost_4); // Left-Bottom-Front Neighbour

      auto ghost_5 = create_ghost_particle(particle_id,
                                           {left_offset[0], 0, back_offset[2]});
      cell_ghost_particles_map[cell_index + x - 1 + (z - 1) * x * y].push_back(
          ghost_5); // Left-Back Neighbour
      cell_ghost_particles_map[cell_index + x - 1 + x + (z - 1) * x * y]
          .push_back(ghost_5); // Left-Top-Back Neighbour

      auto ghost_6 = create_ghost_particle(
          particle_id, {0, bottom_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index + (y - 1) * x + (z - 1) * x * y]
          .push_back(ghost_6); // Bottom-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x + 1 + (z - 1) * x * y]
          .push_back(ghost_6); // Bottom-Right-Back

      auto ghost_7 = create_ghost_particle(
          particle_id, {left_offset[0], bottom_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index + (y - 1) * x + x - 1 +
                               (z - 1) * x * y]
          .push_back(ghost_7); // Bottom-Left-Back
      break;
    }

    case Placement::BOTTOM_BACK_RIGHT_CORNER: {
      auto ghost_1 = create_ghost_particle(particle_id, bottom_offset);
      cell_ghost_particles_map[cell_index + (y - 1) * x].push_back(
          ghost_1); // Bottom
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1].push_back(
          ghost_1); // Bottom-Left
      cell_ghost_particles_map[cell_index + (y - 1) * x + x * y].push_back(
          ghost_1); // Bottom-Front
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1 + x * y].push_back(
          ghost_1); // Bottom-Left-Front

      auto ghost_2 = create_ghost_particle(particle_id, back_offset);
      cell_ghost_particles_map[cell_index + (z - 1) * x * y].push_back(
          ghost_2); // Back
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1].push_back(
          ghost_2); // Back-Left
      cell_ghost_particles_map[cell_index + (z - 1) * x * y + x].push_back(
          ghost_2); // Back-Top
      cell_ghost_particles_map[cell_index + (z - 1) * x * y - 1 + x].push_back(
          ghost_2); // Back-Left-Top

      auto ghost_3 = create_ghost_particle(particle_id, right_offset);
      cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_3); // Right
      cell_ghost_particles_map[cell_index - x + 1 + x].push_back(
          ghost_3); // Right-Top
      cell_ghost_particles_map[cell_index - x + 1 + x * y].push_back(
          ghost_3); // Right-Front
      cell_ghost_particles_map[cell_index - x + 1 + x + x * y].push_back(
          ghost_3); // Right-Top-Front

      auto ghost_4 = create_ghost_particle(
          particle_id, {right_offset[0], bottom_offset[1], 0});
      cell_ghost_particles_map[cell_index - x + 1 + (y - 1) * x].push_back(
          ghost_4); // Right-Bottom
      cell_ghost_particles_map[cell_index - x + 1 + (y - 1) * x + x * y]
          .push_back(ghost_4); // Right-Bottom-Front

      auto ghost_5 = create_ghost_particle(
          particle_id, {right_offset[0], 0, back_offset[2]});
      cell_ghost_particles_map[cell_index - x + 1 + (z - 1) * x * y].push_back(
          ghost_5); // Right-Back
      cell_ghost_particles_map[cell_index - x + 1 + x + (z - 1) * x * y]
          .push_back(ghost_5); // Right-Top-Back

      auto ghost_6 = create_ghost_particle(
          particle_id, {0, bottom_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index + (y - 1) * x + (z - 1) * x * y]
          .push_back(ghost_6); // Bottom-Back
      cell_ghost_particles_map[cell_index + (y - 1) * x - 1 + (z - 1) * x * y]
          .push_back(ghost_6); // Bottom-Left-Back

      auto ghost_7 = create_ghost_particle(
          particle_id, {right_offset[0], bottom_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index + (y - 1) * x - x + 1 +
                               (z - 1) * x * y]
          .push_back(ghost_7); // Bottom-Right-Back
      break;
    }
    default:
      break;
    }
  }
}
