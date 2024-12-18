#include "LinkedCellContainer.h"
#include "../../Simulation.h"
#include "io/input/cli/SimParams.h"
#include <cmath>

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

  placement_map[Placement::TOP] = boundary_conditions.top;
  placement_map[Placement::BOTTOM] = boundary_conditions.bottom;
  placement_map[Placement::LEFT] = boundary_conditions.left;
  placement_map[Placement::RIGHT] = boundary_conditions.right;
  placement_map[Placement::FRONT] = boundary_conditions.front;
  placement_map[Placement::BACK] = boundary_conditions.back;
  
  if(placement_map[Placement::TOP] == placement_map[Placement::RIGHT])
    placement_map[Placement::TOP_RIGHT_CORNER] = placement_map[Placement::RIGHT];
  else 
    placement_map[Placement::TOP_RIGHT_CORNER] = placement_map[Placement::TOP];

  if(placement_map[Placement::TOP] == placement_map[Placement::LEFT])
    placement_map[Placement::TOP_LEFT_CORNER] = placement_map[Placement::LEFT];
  else 
    placement_map[Placement::TOP_LEFT_CORNER] = placement_map[Placement::TOP];

  if(placement_map[Placement::BOTTOM] == placement_map[Placement::RIGHT])
    placement_map[Placement::BOTTOM_RIGHT_CORNER] = placement_map[Placement::RIGHT];
  else 
    placement_map[Placement::BOTTOM_RIGHT_CORNER] = placement_map[Placement::BOTTOM];

  if(placement_map[Placement::BOTTOM] == placement_map[Placement::LEFT])
    placement_map[Placement::BOTTOM_LEFT_CORNER] = placement_map[Placement::LEFT];
  else 
    placement_map[Placement::BOTTOM_LEFT_CORNER] = placement_map[Placement::BOTTOM];
  
  mark_halo_cells();
}

LinkedCellContainer::LinkedCellContainer()
    : domain_size_{0, 0, 0}, r_cutoff_{0}, left_corner_coordinates{0.0, 0.0,
                                                                   0.0},
      x{0}, y{0}, z{0}, boundary_conditions_{}, cells_map{}, particle_id{0},
      particles_left_domain{0}, is_wrapper{false}, halo_count{0},
      reflective_flag{false}, periodic_flag{false}, halo_cell_indices{},
      particles_outbound{}, placement_map{} {}

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

  if (current_index != old_index) {
    if (is_within_domain(old_position)) {
      cells[old_index].remove(particle_id);
    }
    if (is_within_domain(cells_map[particle_id]->getX())) {
      cells[current_index].insert(particle_id);
    } else {
      particles_outbound.push_back(particle_id);
    }
  }
}

std::vector<ParticlePointer>
LinkedCellContainer::get_neighbours(int particle_id, bool check_periodic_neighbours) {
  std::vector<ParticlePointer> neighbours{};
  if (cells_map[particle_id]->left_domain || cells_map[particle_id]->outbound) {
    return neighbours;
  }
  std::array<double, 3> position = cells_map[particle_id]->getX();
  int cell_index = get_cell_index(position);

  // logger.info("Current cell index: " + std::to_string(index));
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

  auto& cell = cells[cell_index];

  // add periodic neighbours
  if(check_periodic_neighbours && cell.is_halo && placement_map[cell.placement] == BoundaryCondition::Periodic) {
    

  }

  return neighbours;
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
  for (size_t i = 0; i < x; i++) {
    for (size_t j = 0; j < y; j++) {
      for (size_t k = 0; k < z; k++) {
        if (i == 0 || i == x - 1 || j == 0 || j == y - 1 ||
            (z > 1 && (k == 0 || k == z - 1))) {

          size_t index = i + j * x + k * x * y;
          cells[index].is_halo = true;
          halo_cell_indices.push_back(index);
          halo_count++;

          if (index == 0) {
            cells[index].placement = Placement::BOTTOM_LEFT_CORNER;
          } else if (index == x - 1) {
            cells[index].placement = Placement::BOTTOM_RIGHT_CORNER;
          } else if (index == x * y - x) {
            cells[index].placement = Placement::TOP_LEFT_CORNER;
          } else if (index == x * y - 1) {
            cells[index].placement = Placement::TOP_RIGHT_CORNER;
          } else if (index < x) {
            cells[index].placement = Placement::BOTTOM;
          } else if (index >= x * (y - 1)) {
            cells[index].placement = Placement::TOP;
          } else if (index % x == 0) {
            cells[index].placement = Placement::LEFT;
          } else if ((index + 1) % x == 0) {
            cells[index].placement = Placement::RIGHT;
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
}