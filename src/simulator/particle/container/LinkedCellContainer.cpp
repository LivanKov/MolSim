#include "LinkedCellContainer.h"
#include "../../Simulation.h"
#include "io/input/cli/SimParams.h"
#include <cmath>
#include <iostream>


bool GhostParticle::operator==(GhostParticle &other) const {
  return id == other.id && position == other.position;
}

size_t LinkedCellContainer::Cell::size() const { return particle_ids.size(); }

void LinkedCellContainer::Cell::insert(int id) { particle_ids.insert(id); }

void LinkedCellContainer::Cell::remove(int id) { particle_ids.erase(id); }

void LinkedCellContainer::initialize(
    const std::initializer_list<double> &domain_size, double r_cutoff,
    const DomainBoundaryConditions &boundary_conditions) {
  if (SimParams::fixed_Domain) {
    left_corner_coordinates = {SimParams::lower_left_corner[0],
                               SimParams::lower_left_corner[1],
                               SimParams::lower_left_corner[2]};
  }

  logger.info("Initializing LinkedCellContainer");
  domain_size_ = std::vector<double>(domain_size);
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

  if (std::abs(remainder_z) > DIVISION_TOLERANCE)
    r_cutoff_z += remainder_z / (std::floor(domain_size_[1] / r_cutoff));

  if (domain_size.size() == 2) {
    domain_size_.push_back(r_cutoff_z);
  }

  // TODO rework this
  x = static_cast<size_t>(domain_size_[0] / r_cutoff);
  y = static_cast<size_t>(domain_size_[1] / r_cutoff);
  z = domain_size.size() == 3
          ? (static_cast<size_t>((domain_size_[2] / r_cutoff)))
          : 1;
  cells = std::vector<Cell>(x * y * z, Cell());

  for (auto &p : particles.get_all_particles()) {
    cells_map[p->getType()] = p;
  }

  set_boundary_conditions(boundary_conditions);
  assign_placements();
  mark_halo_cells();
}

LinkedCellContainer::LinkedCellContainer()
    : domain_size_{0, 0, 0}, r_cutoff_{0},
      left_corner_coordinates{0.0, 0.0, 0.0}, x{0}, y{0}, z{0},
      boundary_conditions_{}, cells_map{}, particle_id{0},
      particles_left_domain{0}, is_wrapper{false}, halo_count{0},
      reflective_flag{false}, periodic_flag{false}, halo_cell_indices{},
      particles_outbound{}, placement_map{}, cell_ghost_particles_map{} {}

void LinkedCellContainer::insert(Particle &p, bool placement) {
  ParticlePointer p_ptr = std::make_shared<Particle>(p);
  if (placement && is_within_domain(p_ptr->getX())) {
    size_t index = get_cell_index(p_ptr->getX());
    cells[index].insert(p_ptr->getType());
  } else if (!is_within_domain(p_ptr->getX())) {
    p_ptr->left_domain = true;
    particles_left_domain++;
  }
  particles.insert(p_ptr);
  cells_map[p_ptr->getType()] = p_ptr;
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
  size_t current_index = get_cell_index(cells_map[particle_id]->getX());

  bool current_within_domain = is_within_domain(cells_map[particle_id]->getX());

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
  if (cells_map[particle_id]->left_domain || cells_map[particle_id]->outbound) {
    return neighbours;
  }
  std::array<double, 3> position = cells_map[particle_id]->getX();
  int cell_index = get_cell_index(position);

  for (auto &i : cells[cell_index].particle_ids) {
    neighbours.push_back(cells_map[i]);
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

        if (ni >= 0 && ni < x && nj >= 0 && nj < y && nk >= 0 && nk < z) {
          int neighborIndex = ni + (nj * x) + nk * x * y;
          for (auto &s : cells[neighborIndex].particle_ids) {
            neighbours.push_back(cells_map[s]);
          }
        }
      }
    }
  }
  return neighbours;
}

std::vector<GhostParticle>
LinkedCellContainer::get_additional_neighbour_indices(int particle_id) {

  auto cell_index = get_cell_index(cells_map[particle_id]->getX());

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

void LinkedCellContainer::mark_halo_cells() {
  for (size_t k = 0; k < z; k++) {
    for (size_t j = 0; j < y; j++) {
      for (size_t i = 0; i < x; i++) {
        if (i == 0 || i == x - 1 || j == 0 || j == y - 1 ||
            (z > 1 && (k == 0 || k == z - 1))) {

          size_t index = i + j * x + k * x * y;
          cells[index].is_halo = true;
          halo_cell_indices.push_back(index);
          halo_count++;
          cells[index].placement = z == 1
                                       ? determine_placement_2d(index, i, j)
                                       : determine_placement_3d(index, i, j, k);
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
    DomainBoundaryConditions boundary_conditions) {
  placement_map[Placement::TOP] = boundary_conditions.top;
  placement_map[Placement::BOTTOM] = boundary_conditions.bottom;
  placement_map[Placement::LEFT] = boundary_conditions.left;
  placement_map[Placement::RIGHT] = boundary_conditions.right;
  if (z != 1) {
    placement_map[Placement::FRONT] = boundary_conditions.front;
    placement_map[Placement::BACK] = boundary_conditions.back;
  }
}

void LinkedCellContainer::clear_ghost_particles() {
  cell_ghost_particles_map.clear();
}

GhostParticle LinkedCellContainer::create_ghost_particle(
    int particle_id, const std::array<double, 3> &position_offset) {
  GhostParticle ghost;
  ghost.sigma = cells_map[particle_id]->getSigma();
  ghost.epsilon = cells_map[particle_id]->getEpsilon();
  ghost.position = {cells_map[particle_id]->getX()[0] + position_offset[0],
                    cells_map[particle_id]->getX()[1] + position_offset[1],
                    cells_map[particle_id]->getX()[2] + position_offset[2]};
  ghost.id = particle_id;
  ghost.ptr = cells_map[particle_id];
  return ghost;
}

void LinkedCellContainer::create_ghost_particles(int particle_id,
                                                 int cell_index) {
  const auto &placement = cells[cell_index].placement;

  // Helper arrays for offsets
  const std::array<double, 3> left_offset = {domain_size_[0], 0, 0};
  const std::array<double, 3> right_offset = {-domain_size_[0], 0, 0};
  const std::array<double, 3> top_offset = {0, -domain_size_[1], 0};
  const std::array<double, 3> bottom_offset = {0, domain_size_[1], 0};
  const std::array<double, 3> front_offset = {0, 0, -domain_size_[2]};
  const std::array<double, 3> back_offset = {0, 0, domain_size_[2]};

  int layer_size = x * y;

  switch (placement) {
  case Placement::LEFT: {
    auto ghost = create_ghost_particle(particle_id, left_offset);
    cell_ghost_particles_map[cell_index + x - 1].push_back(ghost);
    cell_ghost_particles_map[cell_index + x + x - 1].push_back(ghost);
    cell_ghost_particles_map[cell_index - 1].push_back(ghost);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + x - 1 + layer_size].push_back(
          ghost);
      cell_ghost_particles_map[cell_index + x + x - 1 + layer_size].push_back(
          ghost);
      cell_ghost_particles_map[cell_index - 1 + layer_size].push_back(ghost);

      cell_ghost_particles_map[cell_index + x - 1 - layer_size].push_back(
          ghost);
      cell_ghost_particles_map[cell_index + x + x - 1 - layer_size].push_back(
          ghost);
      cell_ghost_particles_map[cell_index - 1 - layer_size].push_back(ghost);
    }
    break;
  }
  case Placement::RIGHT: {
    auto ghost = create_ghost_particle(particle_id, right_offset);
    cell_ghost_particles_map[cell_index - x + 1].push_back(ghost);
    cell_ghost_particles_map[cell_index - (x + x - 1)].push_back(ghost);
    cell_ghost_particles_map[cell_index + 1].push_back(ghost);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - x + 1 + layer_size].push_back(
          ghost);
      cell_ghost_particles_map[cell_index - (x + x - 1) + layer_size].push_back(
          ghost);
      cell_ghost_particles_map[cell_index + 1 + layer_size].push_back(ghost);

      cell_ghost_particles_map[cell_index - x + 1 - layer_size].push_back(
          ghost);
      cell_ghost_particles_map[cell_index - (x + x - 1) - layer_size].push_back(
          ghost);
      cell_ghost_particles_map[cell_index + 1 - layer_size].push_back(ghost);
    }
    break;
  }
  case Placement::TOP: {
    auto ghost = create_ghost_particle(particle_id, top_offset);
    cell_ghost_particles_map[cell_index - (y - 1) * x].push_back(ghost);
    cell_ghost_particles_map[cell_index - ((y - 1) * x) - 1].push_back(ghost);
    cell_ghost_particles_map[cell_index - ((y - 1) * x) + 1].push_back(ghost);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - (y - 1) * x + layer_size].push_back(
          ghost);
      cell_ghost_particles_map[cell_index - ((y - 1) * x) - 1 + layer_size]
          .push_back(ghost);
      cell_ghost_particles_map[cell_index - ((y - 1) * x) + 1 + layer_size]
          .push_back(ghost);

      cell_ghost_particles_map[cell_index - (y - 1) * x - layer_size].push_back(
          ghost);
      cell_ghost_particles_map[cell_index - ((y - 1) * x) - 1 - layer_size]
          .push_back(ghost);
      cell_ghost_particles_map[cell_index - ((y - 1) * x) + 1 - layer_size]
          .push_back(ghost);
    }
    break;
  }
  case Placement::BOTTOM: {
    auto ghost = create_ghost_particle(particle_id, bottom_offset);
    cell_ghost_particles_map[cell_index + (y - 1) * x].push_back(ghost);
    cell_ghost_particles_map[cell_index + ((y - 1) * x) - 1].push_back(ghost);
    cell_ghost_particles_map[cell_index + ((y - 1) * x) + 1].push_back(ghost);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + (y - 1) * x + layer_size].push_back(
          ghost);
      cell_ghost_particles_map[cell_index + ((y - 1) * x) - 1 + layer_size]
          .push_back(ghost);
      cell_ghost_particles_map[cell_index + ((y - 1) * x) + 1 + layer_size]
          .push_back(ghost);

      cell_ghost_particles_map[cell_index + (y - 1) * x - layer_size].push_back(
          ghost);
      cell_ghost_particles_map[cell_index + ((y - 1) * x) - 1 - layer_size]
          .push_back(ghost);
      cell_ghost_particles_map[cell_index + ((y - 1) * x) + 1 - layer_size]
          .push_back(ghost);
    }
    break;
  }
  case Placement::FRONT: {
    if (z != 1) {
      auto ghost_front = create_ghost_particle(particle_id, front_offset);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1)].push_back(
          ghost_front);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1].push_back(
          ghost_front);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 1].push_back(
          ghost_front);

      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x].push_back(
          ghost_front);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1 + x]
          .push_back(ghost_front);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 1 + x]
          .push_back(ghost_front);

      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x].push_back(
          ghost_front);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1 - x]
          .push_back(ghost_front);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 1 - x]
          .push_back(ghost_front);
    }
    break;
  }
  case Placement::BACK: {
    if (z != 1) {
      auto ghost_back = create_ghost_particle(particle_id, back_offset);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1)].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + 1].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - 1].push_back(
          ghost_back);

      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + 1 + x]
          .push_back(ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - 1 + x]
          .push_back(ghost_back);

      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + 1 - x]
          .push_back(ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - 1 - x]
          .push_back(ghost_back);
    }
    break;
  }
  case Placement::BOTTOM_LEFT_CORNER: {
    // Corner ghost
    auto ghost_corner = create_ghost_particle(
        particle_id, {left_offset[0], bottom_offset[1], 0});
    cell_ghost_particles_map[cell_index + y * x - 1].push_back(ghost_corner);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + y * x - 1 + layer_size].push_back(
          ghost_corner);
    }

    // Right side ghost
    auto ghost_right = create_ghost_particle(particle_id, left_offset);
    cell_ghost_particles_map[cell_index + x - 1].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x + x - 1].push_back(ghost_right);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + x + layer_size - 1].push_back(
          ghost_right);
      cell_ghost_particles_map[cell_index + x + x + layer_size - 1].push_back(
          ghost_right);
    }

    // Bottom side ghost
    auto ghost_bottom = create_ghost_particle(particle_id, bottom_offset);
    cell_ghost_particles_map[cell_index + ((y - 1) * x)].push_back(
        ghost_bottom);
    cell_ghost_particles_map[cell_index + ((y - 1) * x) + 1].push_back(
        ghost_bottom);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + ((y - 1) * x) + layer_size]
          .push_back(ghost_bottom);
      cell_ghost_particles_map[cell_index + ((y - 1) * x) + 1 + layer_size]
          .push_back(ghost_bottom);
    }

    if (z != 1) {
      auto ghost_back =
          create_ghost_particle(particle_id, {0, 0, front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1)].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1 + x]
          .push_back(ghost_back);

      auto ghost_diag_right = create_ghost_particle(
          particle_id, {left_offset[0], 0, front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x - 1]
          .push_back(ghost_diag_right);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x - 1 + x]
          .push_back(ghost_diag_right);

      auto ghost_diag_up = create_ghost_particle(
          particle_id, {0, bottom_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x * (y - 1)]
          .push_back(ghost_diag_up);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x * (y - 1) +
                               1]
          .push_back(ghost_diag_up);

      auto ghost_diag = create_ghost_particle(
          particle_id, {left_offset[0], bottom_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * z - 1].push_back(
          ghost_diag);
    }
    break;
  }
  case Placement::TOP_RIGHT_CORNER: {
    // Corner ghost
    auto ghost_corner =
        create_ghost_particle(particle_id, {right_offset[0], top_offset[1], 0});
    cell_ghost_particles_map[cell_index - (y * x - 1)].push_back(ghost_corner);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - (y * x - 1) + layer_size].push_back(
          ghost_corner);
    }

    // Left side ghost
    auto ghost_left = create_ghost_particle(particle_id, right_offset);
    cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - (x + x - 1)].push_back(ghost_left);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - x + 1 + layer_size].push_back(
          ghost_left);
      cell_ghost_particles_map[cell_index - (x + x - 1) + layer_size].push_back(
          ghost_left);
    }

    // Top side ghost
    auto ghost_top = create_ghost_particle(particle_id, top_offset);
    cell_ghost_particles_map[cell_index - ((y - 1) * x)].push_back(ghost_top);
    cell_ghost_particles_map[cell_index - ((y - 1) * x) - 1].push_back(
        ghost_top);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - ((y - 1) * x) + layer_size]
          .push_back(ghost_top);
      cell_ghost_particles_map[cell_index - ((y - 1) * x) - 1 + layer_size]
          .push_back(ghost_top);
    }

    if (z != 1) {
      auto ghost_back =
          create_ghost_particle(particle_id, {0, 0, front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1)].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 1].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 1 - x]
          .push_back(ghost_back);

      auto ghost_diag_left = create_ghost_particle(
          particle_id, {right_offset[0], 0, front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x + 1]
          .push_back(ghost_diag_left);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x + 1 - x]
          .push_back(ghost_diag_left);

      auto ghost_diag_down = create_ghost_particle(
          particle_id, {0, top_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x * (y - 1)]
          .push_back(ghost_diag_down);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x * (y - 1) -
                               1]
          .push_back(ghost_diag_down);

      auto ghost_diag = create_ghost_particle(
          particle_id, {right_offset[0], top_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - y * x + 1]
          .push_back(ghost_diag);
    }
    break;
  }
  case Placement::TOP_LEFT_CORNER: {
    // Corner ghost
    auto ghost_corner =
        create_ghost_particle(particle_id, {left_offset[0], top_offset[1], 0});
    cell_ghost_particles_map[cell_index - (y - 2) * x - 1].push_back(
        ghost_corner);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - (y - 2) * x - 1 + layer_size]
          .push_back(ghost_corner);
    }

    // Right side ghost
    auto ghost_right = create_ghost_particle(particle_id, left_offset);
    cell_ghost_particles_map[cell_index + x - 1].push_back(ghost_right);
    cell_ghost_particles_map[cell_index - 1].push_back(ghost_right);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + x - 1 + layer_size].push_back(
          ghost_right);
      cell_ghost_particles_map[cell_index - 1 + layer_size].push_back(
          ghost_right);
    }

    // Top side ghost
    auto ghost_top = create_ghost_particle(particle_id, top_offset);
    cell_ghost_particles_map[cell_index - ((y - 1) * x)].push_back(ghost_top);
    cell_ghost_particles_map[cell_index - ((y - 1) * x) + 1].push_back(
        ghost_top);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - ((y - 1) * x) + layer_size]
          .push_back(ghost_top);
      cell_ghost_particles_map[cell_index - ((y - 1) * x) + 1 + layer_size]
          .push_back(ghost_top);
    }

    if (z != 1) {
      auto ghost_back =
          create_ghost_particle(particle_id, {0, 0, front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1)].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x + 1]
          .push_back(ghost_back);

      auto ghost_diag_right = create_ghost_particle(
          particle_id, {left_offset[0], 0, front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x - 1]
          .push_back(ghost_diag_right);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x - 1 - x]
          .push_back(ghost_diag_right);

      auto ghost_diag_down = create_ghost_particle(
          particle_id, {0, top_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x * (y - 1)]
          .push_back(ghost_diag_down);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x * (y - 1) +
                               1]
          .push_back(ghost_diag_down);

      auto ghost_diag = create_ghost_particle(
          particle_id, {left_offset[0], top_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - (y - 1) * x +
                               x - 1]
          .push_back(ghost_diag);
    }

    break;
  }
  case Placement::BOTTOM_RIGHT_CORNER: {
    // Corner ghost
    auto ghost_corner = create_ghost_particle(
        particle_id, {right_offset[0], bottom_offset[1], 0});
    cell_ghost_particles_map[cell_index + (y - 2) * x + 1].push_back(
        ghost_corner);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + (y - 2) * x + 1 + layer_size]
          .push_back(ghost_corner);
    }

    // Left side ghost
    auto ghost_left = create_ghost_particle(particle_id, right_offset);
    cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_left);
    cell_ghost_particles_map[cell_index + 1].push_back(ghost_left);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - x + 1 + layer_size].push_back(
          ghost_left);
      cell_ghost_particles_map[cell_index + 1 + layer_size].push_back(
          ghost_left);
    }

    // Bottom side ghost
    auto ghost_bottom = create_ghost_particle(particle_id, bottom_offset);
    cell_ghost_particles_map[cell_index + ((y - 1) * x)].push_back(
        ghost_bottom);
    cell_ghost_particles_map[cell_index + ((y - 1) * x) - 1].push_back(
        ghost_bottom);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + ((y - 1) * x) + layer_size]
          .push_back(ghost_bottom);
      cell_ghost_particles_map[cell_index + ((y - 1) * x) - 1 + layer_size]
          .push_back(ghost_bottom);
    }

    if (z != 1) {
      auto ghost_back =
          create_ghost_particle(particle_id, {0, 0, front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1)].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 1].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x - 1]
          .push_back(ghost_back);

      auto ghost_diag_left = create_ghost_particle(
          particle_id, {right_offset[0], 0, front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x + 1]
          .push_back(ghost_diag_left);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1].push_back(
          ghost_diag_left);

      auto ghost_diag_up = create_ghost_particle(
          particle_id, {0, bottom_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x * (y - 1)]
          .push_back(ghost_diag_up);
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x * (y - 1) -
                               1]
          .push_back(ghost_diag_up);

      auto ghost_diag = create_ghost_particle(
          particle_id, {right_offset[0], bottom_offset[1], front_offset[2]});
      cell_ghost_particles_map[cell_index + layer_size * (z - 1) + (y - 1) * x -
                               x + 1]
          .push_back(ghost_diag);
    }
    break;
  }
  case Placement::TOP_LEFT_CORNER_BACK: {
    auto ghost_corner =
        create_ghost_particle(particle_id, {left_offset[0], top_offset[1], 0});
    cell_ghost_particles_map[cell_index - (y - 2) * x - 1].push_back(
        ghost_corner);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - (y - 2) * x - 1 - layer_size]
          .push_back(ghost_corner);
    }

    // Right side ghost
    auto ghost_right = create_ghost_particle(particle_id, left_offset);
    cell_ghost_particles_map[cell_index + x - 1].push_back(ghost_right);
    cell_ghost_particles_map[cell_index - 1].push_back(ghost_right);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + x - 1 - layer_size].push_back(
          ghost_right);
      cell_ghost_particles_map[cell_index - 1 - layer_size].push_back(
          ghost_right);
    }

    // bottom side ghost
    auto ghost_top = create_ghost_particle(particle_id, top_offset);
    cell_ghost_particles_map[cell_index - ((y - 1) * x)].push_back(ghost_top);
    cell_ghost_particles_map[cell_index - ((y - 1) * x) + 1].push_back(
        ghost_top);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - ((y - 1) * x) - layer_size]
          .push_back(ghost_top);
      cell_ghost_particles_map[cell_index - ((y - 1) * x) + 1 - layer_size]
          .push_back(ghost_top);
    }

    if (z != 1) {
      auto ghost_back =
          create_ghost_particle(particle_id, {0, 0, back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1)].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + 1].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x + 1]
          .push_back(ghost_back);

      auto ghost_diag_right = create_ghost_particle(
          particle_id, {left_offset[0], 0, back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x - 1]
          .push_back(ghost_diag_right);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - 1].push_back(
          ghost_diag_right);

      auto ghost_diag_down = create_ghost_particle(
          particle_id, {0, top_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x * (y - 1)]
          .push_back(ghost_diag_down);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x * (y - 1) +
                               1]
          .push_back(ghost_diag_down);

      auto ghost_diag = create_ghost_particle(
          particle_id, {left_offset[0], top_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - (y - 1) * x +
                               x - 1]
          .push_back(ghost_diag);
    }

    break;
  }
  case Placement::BOTTOM_LEFT_CORNER_BACK: {
    auto ghost_corner = create_ghost_particle(
        particle_id, {left_offset[0], bottom_offset[1], 0});
    cell_ghost_particles_map[cell_index + y * x - 1].push_back(ghost_corner);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + y * x - 1 - layer_size].push_back(
          ghost_corner);
    }

    // Right side ghost
    auto ghost_right = create_ghost_particle(particle_id, left_offset);
    cell_ghost_particles_map[cell_index + x - 1].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x + x - 1].push_back(ghost_right);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + x - 1 - layer_size].push_back(
          ghost_right);
      cell_ghost_particles_map[cell_index + x + x - 1 - layer_size].push_back(
          ghost_right);
    }

    // Upper side ghost
    auto ghost_bottom = create_ghost_particle(particle_id, bottom_offset);
    cell_ghost_particles_map[cell_index + ((y - 1) * x)].push_back(
        ghost_bottom);
    cell_ghost_particles_map[cell_index + ((y - 1) * x) + 1].push_back(
        ghost_bottom);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + ((y - 1) * x) - layer_size]
          .push_back(ghost_bottom);
      cell_ghost_particles_map[cell_index + ((y - 1) * x) + 1 - layer_size]
          .push_back(ghost_bottom);
    }

    if (z != 1) {
      auto ghost_front =
          create_ghost_particle(particle_id, {0, 0, back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1)].push_back(
          ghost_front);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + 1].push_back(
          ghost_front);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x].push_back(
          ghost_front);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + 1 + x]
          .push_back(ghost_front);

      auto ghost_diag_right = create_ghost_particle(
          particle_id, {left_offset[0], 0, back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x - 1]
          .push_back(ghost_diag_right);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x - 1 + x]
          .push_back(ghost_diag_right);

      auto ghost_diag_up = create_ghost_particle(
          particle_id, {0, bottom_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x * (y - 1)]
          .push_back(ghost_diag_up);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x * (y - 1) +
                               1]
          .push_back(ghost_diag_up);

      auto ghost_diag = create_ghost_particle(
          particle_id, {left_offset[0], bottom_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x * y - 1]
          .push_back(ghost_diag);
    }
    break;
  }
  case Placement::BOTTOM_RIGHT_CORNER_BACK: {
    auto ghost_corner = create_ghost_particle(
        particle_id, {right_offset[0], bottom_offset[1], 0});
    cell_ghost_particles_map[cell_index + (y - 2) * x + 1].push_back(
        ghost_corner);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + (y - 2) * x + 1 - layer_size]
          .push_back(ghost_corner);
    }

    // Left side ghost
    auto ghost_left = create_ghost_particle(particle_id, right_offset);
    cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_left);
    cell_ghost_particles_map[cell_index + 1].push_back(ghost_left);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - x + 1 - layer_size].push_back(
          ghost_left);
      cell_ghost_particles_map[cell_index + 1 - layer_size].push_back(
          ghost_left);
    }

    // Bottom side ghost
    auto ghost_bottom = create_ghost_particle(particle_id, bottom_offset);
    cell_ghost_particles_map[cell_index + ((y - 1) * x)].push_back(
        ghost_bottom);
    cell_ghost_particles_map[cell_index + ((y - 1) * x) - 1].push_back(
        ghost_bottom);
    if (z != 1) {
      cell_ghost_particles_map[cell_index + ((y - 1) * x) - layer_size]
          .push_back(ghost_bottom);
      cell_ghost_particles_map[cell_index + ((y - 1) * x) - 1 - layer_size]
          .push_back(ghost_bottom);
    }

    if (z != 1) {
      auto ghost_back =
          create_ghost_particle(particle_id, {0, 0, back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1)].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - 1].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x - 1]
          .push_back(ghost_back);

      auto ghost_diag_left = create_ghost_particle(
          particle_id, {right_offset[0], 0, back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x + 1]
          .push_back(ghost_diag_left);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + 1].push_back(
          ghost_diag_left);

      auto ghost_diag_up = create_ghost_particle(
          particle_id, {0, bottom_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x * (y - 1)]
          .push_back(ghost_diag_up);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x * (y - 1) -
                               1]
          .push_back(ghost_diag_up);

      auto ghost_diag = create_ghost_particle(
          particle_id, {right_offset[0], bottom_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) + (y - 1) * x -
                               x + 1]
          .push_back(ghost_diag);
    }
    break;
  }

  case Placement::TOP_RIGHT_CORNER_BACK: {

    auto ghost_corner =
        create_ghost_particle(particle_id, {right_offset[0], top_offset[1], 0});
    cell_ghost_particles_map[cell_index - (y * x - 1)].push_back(ghost_corner);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - (y * x - 1) - layer_size].push_back(
          ghost_corner);
    }

    // Left side ghost
    auto ghost_left = create_ghost_particle(particle_id, right_offset);
    cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - (x + x - 1)].push_back(ghost_left);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - x + 1 - layer_size].push_back(
          ghost_left);
      cell_ghost_particles_map[cell_index - (x + x - 1) - layer_size].push_back(
          ghost_left);
    }

    // Top side ghost
    auto ghost_top = create_ghost_particle(particle_id, top_offset);
    cell_ghost_particles_map[cell_index - ((y - 1) * x)].push_back(ghost_top);
    cell_ghost_particles_map[cell_index - ((y - 1) * x) - 1].push_back(
        ghost_top);
    if (z != 1) {
      cell_ghost_particles_map[cell_index - ((y - 1) * x) - layer_size]
          .push_back(ghost_top);
      cell_ghost_particles_map[cell_index - ((y - 1) * x) - 1 - layer_size]
          .push_back(ghost_top);
    }

    if (z != 1) {
      auto ghost_back =
          create_ghost_particle(particle_id, {0, 0, back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1)].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - 1].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x].push_back(
          ghost_back);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - 1 - x]
          .push_back(ghost_back);

      auto ghost_diag_left = create_ghost_particle(
          particle_id, {right_offset[0], 0, back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x + 1]
          .push_back(ghost_diag_left);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x + 1 - x]
          .push_back(ghost_diag_left);

      auto ghost_diag_down = create_ghost_particle(
          particle_id, {0, top_offset[1], back_offset[2]});
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x * (y - 1)]
          .push_back(ghost_diag_down);
      cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x * (y - 1) -
                               1]
          .push_back(ghost_diag_down);

      auto ghost_diag = create_ghost_particle(
          particle_id, {right_offset[0], top_offset[1], back_offset[2]});
      cell_ghost_particles_map[0].push_back(ghost_diag);
    }
    break;
  }

  case Placement::LEFT_TOP: {
    auto ghost_right = create_ghost_particle(particle_id, left_offset);
    cell_ghost_particles_map[cell_index + x - 1].push_back(ghost_right);
    cell_ghost_particles_map[cell_index - 1].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 + layer_size].push_back(
        ghost_right);
    cell_ghost_particles_map[cell_index - 1 + layer_size].push_back(
        ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 - layer_size].push_back(
        ghost_right);
    cell_ghost_particles_map[cell_index - 1 - layer_size].push_back(
        ghost_right);

    auto ghost_bottom = create_ghost_particle(particle_id, top_offset);
    cell_ghost_particles_map[cell_index - x * (y - 1)].push_back(ghost_bottom);
    cell_ghost_particles_map[cell_index - x * (y - 1) + 1].push_back(
        ghost_bottom);
    cell_ghost_particles_map[cell_index - x * (y - 1) + layer_size].push_back(
        ghost_bottom);
    cell_ghost_particles_map[cell_index - x * (y - 1) + 1 + layer_size]
        .push_back(ghost_bottom);
    cell_ghost_particles_map[cell_index - x * (y - 1) - layer_size].push_back(
        ghost_bottom);
    cell_ghost_particles_map[cell_index - x * (y - 1) + 1 - layer_size]
        .push_back(ghost_bottom);

    auto ghost_corner =
        create_ghost_particle(particle_id, {left_offset[0], top_offset[1], 0});
    cell_ghost_particles_map[cell_index - x * (y - 1) + x - 1].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index - x * (y - 1) + x - 1 + layer_size]
        .push_back(ghost_corner);
    cell_ghost_particles_map[cell_index - x * (y - 1) + x - 1 - layer_size]
        .push_back(ghost_corner);
    break;
  }
  case Placement::LEFT_BOTTOM: {

    auto ghost_right = create_ghost_particle(particle_id, left_offset);
    cell_ghost_particles_map[cell_index + x - 1].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + 2 * x - 1].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 + layer_size].push_back(
        ghost_right);
    cell_ghost_particles_map[cell_index + 2 * x - 1 + layer_size].push_back(
        ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 - layer_size].push_back(
        ghost_right);
    cell_ghost_particles_map[cell_index + 2 * x - 1 - layer_size].push_back(
        ghost_right);

    auto ghost_top = create_ghost_particle(particle_id, bottom_offset);
    cell_ghost_particles_map[cell_index + x * (y - 1)].push_back(ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) + 1].push_back(ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) + layer_size].push_back(
        ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) + 1 + layer_size]
        .push_back(ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) - layer_size].push_back(
        ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) + 1 - layer_size]
        .push_back(ghost_top);

    auto ghost_corner = create_ghost_particle(
        particle_id, {left_offset[0], bottom_offset[1], 0});
    cell_ghost_particles_map[cell_index + x * y - 1].push_back(ghost_corner);
    cell_ghost_particles_map[cell_index + x * y - 1 + layer_size].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index + x * y - 1 - layer_size].push_back(
        ghost_corner);
    break;
  }
  case Placement::RIGHT_BOTTOM: {
    auto ghost_left = create_ghost_particle(particle_id, right_offset);
    cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_left);
    cell_ghost_particles_map[cell_index + 1].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 + layer_size].push_back(
        ghost_left);
    cell_ghost_particles_map[cell_index + 1 + layer_size].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 - layer_size].push_back(
        ghost_left);
    cell_ghost_particles_map[cell_index + 1 - layer_size].push_back(ghost_left);

    auto ghost_top = create_ghost_particle(particle_id, bottom_offset);
    cell_ghost_particles_map[cell_index + x * (y - 1)].push_back(ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) - 1].push_back(ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) + layer_size].push_back(
        ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) - 1 + layer_size]
        .push_back(ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) - layer_size].push_back(
        ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) - 1 - layer_size]
        .push_back(ghost_top);

    auto ghost_corner = create_ghost_particle(
        particle_id, {right_offset[0], bottom_offset[1], 0});
    cell_ghost_particles_map[cell_index + 1 + x * (y - 2)].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index + 1 + x * (y - 2) + layer_size]
        .push_back(ghost_corner);
    cell_ghost_particles_map[cell_index + 1 + x * (y - 2) - layer_size]
        .push_back(ghost_corner);
    break;
  }
  case Placement::RIGHT_TOP: {
    auto ghost_left = create_ghost_particle(particle_id, right_offset);
    cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - 2 * x + 1].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 + layer_size].push_back(
        ghost_left);
    cell_ghost_particles_map[cell_index - 2 * x + 1 + layer_size].push_back(
        ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 - layer_size].push_back(
        ghost_left);
    cell_ghost_particles_map[cell_index - 2 * x + 1 - layer_size].push_back(
        ghost_left);

    auto ghost_top = create_ghost_particle(particle_id, top_offset);
    cell_ghost_particles_map[cell_index - x * (y - 1)].push_back(ghost_top);
    cell_ghost_particles_map[cell_index - x * (y - 1) - 1].push_back(ghost_top);
    cell_ghost_particles_map[cell_index - x * (y - 1) + layer_size].push_back(
        ghost_top);
    cell_ghost_particles_map[cell_index - x * (y - 1) - 1 + layer_size]
        .push_back(ghost_top);
    cell_ghost_particles_map[cell_index - x * (y - 1) - layer_size].push_back(
        ghost_top);
    cell_ghost_particles_map[cell_index - x * (y - 1) - 1 - layer_size]
        .push_back(ghost_top);

    auto ghost_corner = create_ghost_particle(
        particle_id, {right_offset[0], top_offset[1], 0});
    cell_ghost_particles_map[cell_index - x * (y - 1) - x + 1].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index - x * (y - 1) - x + 1 + layer_size]
        .push_back(ghost_corner);
    cell_ghost_particles_map[cell_index - x * (y - 1) - x + 1 - layer_size]
        .push_back(ghost_corner);
    break;
  }
  case Placement::TOP_FRONT: {
    auto ghost_back = create_ghost_particle(particle_id, front_offset);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1)].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x + 1].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 1].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x - 1].push_back(ghost_back);

    auto ghost_bottom = create_ghost_particle(particle_id, top_offset);
    cell_ghost_particles_map[cell_index - x * (y - 1)].push_back(ghost_bottom);
    cell_ghost_particles_map[cell_index - x * (y - 1) + layer_size].push_back(ghost_bottom);
    cell_ghost_particles_map[cell_index - x * (y - 1) + 1].push_back(ghost_bottom);
    cell_ghost_particles_map[cell_index - x * (y - 1) + layer_size + 1].push_back(ghost_bottom);
    cell_ghost_particles_map[cell_index - x * (y - 1) - 1].push_back(ghost_bottom);
    cell_ghost_particles_map[cell_index - x * (y - 1) + layer_size - 1].push_back(ghost_bottom);

    auto ghost_corner = create_ghost_particle(
        particle_id, {0, top_offset[1], front_offset[2]});
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x * (y - 1)].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x * (y - 1) + 1]
        .push_back(ghost_corner);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x * (y - 1) - 1]
        .push_back(ghost_corner);
    break;
  }
  case Placement::BOTTOM_FRONT: {
    auto ghost_back = create_ghost_particle(particle_id, front_offset);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1)].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x + 1].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 1].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x - 1].push_back(ghost_back);

    auto ghost_top = create_ghost_particle(particle_id, bottom_offset);
    cell_ghost_particles_map[cell_index + x * (y - 1)].push_back(ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) + layer_size].push_back(ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) + 1].push_back(ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) + layer_size + 1].push_back(ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) - 1].push_back(ghost_top);
    cell_ghost_particles_map[cell_index + x * (y - 1) + layer_size - 1].push_back(ghost_top);

    auto ghost_corner = create_ghost_particle(
        particle_id, {0, bottom_offset[1], front_offset[2]});
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x * (y - 1)].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x * (y - 1) + 1]
        .push_back(ghost_corner);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x * (y - 1) - 1]
        .push_back(ghost_corner);
    break;
  }
  case Placement::BOTTOM_BACK: {
    auto ghost_front = create_ghost_particle(particle_id, back_offset);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1)].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + 1].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x + 1].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - 1].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x - 1].push_back(ghost_front);

    auto ghost_up = create_ghost_particle(particle_id, bottom_offset);
    cell_ghost_particles_map[cell_index + x * (y - 1)].push_back(ghost_up);
    cell_ghost_particles_map[cell_index + x * (y - 1) - layer_size].push_back(ghost_up);
    cell_ghost_particles_map[cell_index + x * (y - 1) + 1].push_back(ghost_up);
    cell_ghost_particles_map[cell_index + x * (y - 1) - layer_size + 1].push_back(ghost_up);
    cell_ghost_particles_map[cell_index + x * (y - 1) - 1].push_back(ghost_up);
    cell_ghost_particles_map[cell_index + x * (y - 1) - layer_size - 1].push_back(ghost_up);

    auto ghost_corner = create_ghost_particle(
        particle_id, {0, bottom_offset[1], back_offset[2]});
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x * (y - 1)].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x * (y - 1) + 1]
        .push_back(ghost_corner);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x * (y - 1) - 1]
        .push_back(ghost_corner);
    break;
  }
  case Placement::TOP_BACK: {
    auto ghost_front = create_ghost_particle(particle_id, back_offset);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1)].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + 1].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x + 1].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - 1].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x - 1].push_back(ghost_front);

    auto ghost_down = create_ghost_particle(particle_id, top_offset);
    cell_ghost_particles_map[cell_index - x * (y - 1)].push_back(ghost_down);
    cell_ghost_particles_map[cell_index - x * (y - 1) - layer_size].push_back(ghost_down);
    cell_ghost_particles_map[cell_index - x * (y - 1) + 1].push_back(ghost_down);
    cell_ghost_particles_map[cell_index - x * (y - 1) - layer_size + 1].push_back(ghost_down);
    cell_ghost_particles_map[cell_index - x * (y - 1) - 1].push_back(ghost_down);
    cell_ghost_particles_map[cell_index - x * (y - 1) - layer_size - 1].push_back(ghost_down);

    auto ghost_corner = create_ghost_particle(
        particle_id, {0, bottom_offset[1], back_offset[2]});
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x * (y - 1)].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x * (y - 1) + 1]
        .push_back(ghost_corner);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x * (y - 1) - 1]
        .push_back(ghost_corner);
    break;
  }

  case Placement::RIGHT_FRONT: {
    auto ghost_back = create_ghost_particle(particle_id, front_offset);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1)].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 1].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 1 + x].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 1 - x].push_back(ghost_back);

    auto ghost_left = create_ghost_particle(particle_id, right_offset);
    cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 + layer_size].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 + x].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 + layer_size + x].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 - x].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 + layer_size - x].push_back(ghost_left);

    auto ghost_corner = create_ghost_particle(
        particle_id, {0, right_offset[1], front_offset[2]});
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x + 1].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1]
        .push_back(ghost_corner);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 2 * x + 1]
        .push_back(ghost_corner);
    break;
  }
  case Placement::LEFT_FRONT: {
    auto ghost_back = create_ghost_particle(particle_id, front_offset);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1)].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1 + x].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - x].push_back(ghost_back);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + 1 - x].push_back(ghost_back);

    auto ghost_right = create_ghost_particle(particle_id, left_offset);
    cell_ghost_particles_map[cell_index + x - 1].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 + layer_size].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 + x].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 + layer_size + x].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 - x].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 + layer_size - x].push_back(ghost_right);

    auto ghost_corner = create_ghost_particle(
        particle_id, {0, left_offset[1], front_offset[2]});
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x - 1].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) + x - 1 + x].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index + layer_size * (z - 1) - 1].push_back(
        ghost_corner);
    break;
  }
  case Placement::LEFT_BACK: {
    auto ghost_front = create_ghost_particle(particle_id, back_offset);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1)].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + 1].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + 1 + x].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + 1 - x].push_back(ghost_front);

    auto ghost_right = create_ghost_particle(particle_id, left_offset);
    cell_ghost_particles_map[cell_index + x - 1].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 - layer_size].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 + x].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 - layer_size + x].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 - x].push_back(ghost_right);
    cell_ghost_particles_map[cell_index + x - 1 - layer_size - x].push_back(ghost_right);

    auto ghost_corner = create_ghost_particle(
        particle_id, {0, left_offset[1], back_offset[2]});
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x - 1].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x - 1 + x].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x - 1 - x].push_back(
        ghost_corner);
    break;
  }
  case Placement::RIGHT_BACK: {
    auto ghost_front = create_ghost_particle(particle_id, back_offset);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1)].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - 1].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) + x].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - 1 + x].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x].push_back(ghost_front);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - 1 - x].push_back(ghost_front);

    auto ghost_left = create_ghost_particle(particle_id, right_offset);
    cell_ghost_particles_map[cell_index - x + 1].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 - layer_size].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 + x].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 - layer_size + x].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 - x].push_back(ghost_left);
    cell_ghost_particles_map[cell_index - x + 1 - layer_size - x].push_back(ghost_left);

    auto ghost_corner = create_ghost_particle(
        particle_id, {0, right_offset[1], back_offset[2]});
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x + 1].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x + 1 + x].push_back(
        ghost_corner);
    cell_ghost_particles_map[cell_index - layer_size * (z - 1) - x + 1 - x].push_back(
        ghost_corner);
    break;
  }
  default:
    break;
  }
}

Placement LinkedCellContainer::determine_placement_3d(size_t index, size_t i,
                                                      size_t j, size_t k) {
  bool x_min = (i == 0);
  bool x_max = (i == x - 1);
  bool y_min = (j == 0);
  bool y_max = (j == y - 1);
  bool z_min = (k == 0);
  bool z_max = (k == z - 1);

  int count_extremes = (x_min || x_max) + (y_min || y_max) + (z_min || z_max);

  if (count_extremes == 3) {
    // Corner case
    if (z_min) {
      if (x_min) {
        return y_min ? BOTTOM_LEFT_CORNER : TOP_LEFT_CORNER;
      } else {
        return y_min ? BOTTOM_RIGHT_CORNER : TOP_RIGHT_CORNER;
      }
    } else {
      if (x_min) {
        return y_min ? BOTTOM_LEFT_CORNER_BACK : TOP_LEFT_CORNER_BACK;
      } else {
        return y_min ? BOTTOM_RIGHT_CORNER_BACK : TOP_RIGHT_CORNER_BACK;
      }
    }
  } else if (count_extremes == 2) {
    // Edge case
    bool x_extreme = x_min || x_max;
    bool y_extreme = y_min || y_max;
    bool z_extreme = z_min || z_max;

    if (x_extreme && z_extreme) {
      if (x_max) {
        return z_min ? RIGHT_FRONT : RIGHT_BACK;
      } else {
        return z_min ? LEFT_FRONT : LEFT_BACK;
      }
    } else if (x_extreme && y_extreme) {
      if (x_max) {
        return y_max ? RIGHT_TOP : RIGHT_BOTTOM;
      } else {
        return y_max ? LEFT_TOP : LEFT_BOTTOM;
      }
    } else {
      if (y_max) {
        return z_min ? TOP_FRONT : TOP_BACK;
      } else {
        return z_min ? BOTTOM_FRONT : BOTTOM_BACK;
      }
    }
  } else {
    if (x_min)
      return LEFT;
    if (x_max)
      return RIGHT;
    if (y_min)
      return BOTTOM;
    if (y_max)
      return TOP;
    if (z_min)
      return FRONT;
    return BACK;
  }
}

Placement LinkedCellContainer::determine_placement_2d(size_t index, size_t i,
                                                      size_t j) {
  bool x_min = (i == 0);
  bool x_max = (i == x - 1);
  bool y_min = (j == 0);
  bool y_max = (j == y - 1);

  int count_extremes = (x_min || x_max) + (y_min || y_max);

  if (count_extremes == 2) {

    if (x_min) {
      if (y_max)
        return TOP_LEFT_CORNER;
      else
        return BOTTOM_LEFT_CORNER;
    } else if (x_max) {
      if (y_max)
        return TOP_RIGHT_CORNER;
      else
        return BOTTOM_RIGHT_CORNER;
    }
  } else {
    if (x_min)
      return LEFT;
    if (x_max)
      return RIGHT;
    if (y_min)
      return BOTTOM;
    if (y_max)
      return TOP;
  }
}

void LinkedCellContainer::assign_placements() {

  if (z == 1) {

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

  } else {

    if (placement_map[Placement::TOP] == placement_map[Placement::RIGHT] ==
        placement_map[Placement::FRONT])
      placement_map[Placement::TOP_RIGHT_CORNER] =
          placement_map[Placement::RIGHT];
    else
      placement_map[Placement::TOP_RIGHT_CORNER] =
          placement_map[Placement::TOP];

    if (placement_map[Placement::TOP] == placement_map[Placement::LEFT] ==
        placement_map[Placement::FRONT])
      placement_map[Placement::TOP_LEFT_CORNER] =
          placement_map[Placement::LEFT];
    else
      placement_map[Placement::TOP_LEFT_CORNER] = placement_map[Placement::TOP];

    if (placement_map[Placement::BOTTOM] == placement_map[Placement::RIGHT] ==
        placement_map[Placement::FRONT])
      placement_map[Placement::BOTTOM_RIGHT_CORNER] =
          placement_map[Placement::RIGHT];
    else
      placement_map[Placement::BOTTOM_RIGHT_CORNER] =
          placement_map[Placement::BOTTOM];

    if (placement_map[Placement::BOTTOM] == placement_map[Placement::LEFT] ==
        placement_map[Placement::FRONT])
      placement_map[Placement::BOTTOM_LEFT_CORNER] =
          placement_map[Placement::LEFT];
    else
      placement_map[Placement::BOTTOM_LEFT_CORNER] =
          placement_map[Placement::BOTTOM];

    if (placement_map[Placement::TOP] == placement_map[Placement::RIGHT] ==
        placement_map[Placement::BACK])
      placement_map[Placement::TOP_RIGHT_CORNER_BACK] =
          placement_map[Placement::RIGHT];
    else
      placement_map[Placement::TOP_RIGHT_CORNER_BACK] =
          placement_map[Placement::TOP];

    if (placement_map[Placement::TOP] == placement_map[Placement::LEFT] ==
        placement_map[Placement::BACK])
      placement_map[Placement::TOP_LEFT_CORNER_BACK] =
          placement_map[Placement::LEFT];
    else
      placement_map[Placement::TOP_LEFT_CORNER_BACK] =
          placement_map[Placement::TOP];

    if (placement_map[Placement::BOTTOM] == placement_map[Placement::RIGHT] ==
        placement_map[Placement::BACK])
      placement_map[Placement::BOTTOM_RIGHT_CORNER_BACK] =
          placement_map[Placement::RIGHT];
    else
      placement_map[Placement::BOTTOM_RIGHT_CORNER_BACK] =
          placement_map[Placement::BOTTOM];

    if (placement_map[Placement::BOTTOM] == placement_map[Placement::LEFT] ==
        placement_map[Placement::BACK])
      placement_map[Placement::BOTTOM_LEFT_CORNER_BACK] =
          placement_map[Placement::LEFT];
    else
      placement_map[Placement::BOTTOM_LEFT_CORNER_BACK] =
          placement_map[Placement::BOTTOM];

    if (placement_map[Placement::RIGHT] == placement_map[Placement::TOP])
      placement_map[Placement::RIGHT_TOP] = placement_map[Placement::RIGHT];
    else
      placement_map[Placement::RIGHT_TOP] = placement_map[Placement::TOP];

    if (placement_map[Placement::RIGHT] == placement_map[Placement::BOTTOM])
      placement_map[Placement::RIGHT_BOTTOM] = placement_map[Placement::RIGHT];
    else
      placement_map[Placement::RIGHT_BOTTOM] = placement_map[Placement::BOTTOM];

    if (placement_map[Placement::LEFT] == placement_map[Placement::TOP])
      placement_map[Placement::LEFT_TOP] = placement_map[Placement::LEFT];
    else
      placement_map[Placement::LEFT_TOP] = placement_map[Placement::TOP];

    if (placement_map[Placement::LEFT] == placement_map[Placement::BOTTOM])
      placement_map[Placement::LEFT_BOTTOM] = placement_map[Placement::LEFT];
    else
      placement_map[Placement::LEFT_BOTTOM] = placement_map[Placement::BOTTOM];
    
    if (placement_map[Placement::TOP] == placement_map[Placement::FRONT])
      placement_map[Placement::TOP_FRONT] = placement_map[Placement::TOP];
    else
      placement_map[Placement::TOP_FRONT] = placement_map[Placement::TOP];
    
    if(placement_map[Placement::BOTTOM] == placement_map[Placement::FRONT])
      placement_map[Placement::BOTTOM_FRONT] = placement_map[Placement::BOTTOM];
    else
      placement_map[Placement::BOTTOM_FRONT] = placement_map[Placement::BOTTOM];
    
    if(placement_map[Placement::TOP] == placement_map[Placement::BACK])
      placement_map[Placement::TOP_BACK] = placement_map[Placement::TOP];
    else
      placement_map[Placement::TOP_BACK] = placement_map[Placement::TOP];
    
    if(placement_map[Placement::BOTTOM] == placement_map[Placement::BACK])
      placement_map[Placement::BOTTOM_BACK] = placement_map[Placement::BOTTOM];
    else
      placement_map[Placement::BOTTOM_BACK] = placement_map[Placement::BOTTOM];

    if(placement_map[Placement::RIGHT] == placement_map[Placement::FRONT])
      placement_map[Placement::RIGHT_FRONT] = placement_map[Placement::RIGHT];
    else
      placement_map[Placement::RIGHT_FRONT] = placement_map[Placement::FRONT];

    if(placement_map[Placement::LEFT] == placement_map[Placement::FRONT])
      placement_map[Placement::LEFT_FRONT] = placement_map[Placement::LEFT];
    else
      placement_map[Placement::LEFT_FRONT] = placement_map[Placement::FRONT];

    if(placement_map[Placement::RIGHT] == placement_map[Placement::BACK])
      placement_map[Placement::RIGHT_BACK] = placement_map[Placement::RIGHT];
    else
      placement_map[Placement::RIGHT_BACK] = placement_map[Placement::BACK];

    if(placement_map[Placement::LEFT] == placement_map[Placement::BACK])
      placement_map[Placement::LEFT_BACK] = placement_map[Placement::LEFT];
    else
      placement_map[Placement::LEFT_BACK] = placement_map[Placement::BACK];
  }
}