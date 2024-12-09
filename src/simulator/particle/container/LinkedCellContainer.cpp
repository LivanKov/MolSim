#include "LinkedCellContainer.h"
#include <cmath>

size_t Cell::size() const { return particles.size(); }

ParticlePointer Cell::operator[](size_t index) { return particles[index]; }

LinkedCellContainer::LinkedCellContainer(
    std::initializer_list<double> domain_size, double r_cutoff,
    const DomainBoundaryConditions &boundary_conditions)
    : domain_size_{domain_size}, r_cutoff_{r_cutoff},
      r_cutoff_x{r_cutoff}, r_cutoff_y{r_cutoff}, r_cutoff_z{r_cutoff},
      left_corner_coordinates{0.0, 0.0, 0.0}, x{0}, y{0}, z{0},
      boundary_conditions_{boundary_conditions}, particles{} {
  if (domain_size.size() != 3 && domain_size.size() != 2) {
    throw std::invalid_argument("Domain size must have 2 or 3 elements");
  }

  double remainder_x = std::fmod(domain_size_[0], r_cutoff);
  double remainder_y = std::fmod(domain_size_[1], r_cutoff);
  double remainder_z =
      domain_size.size() == 3 ? std::fmod(domain_size_[2], r_cutoff) : 0.0;
  if (std::abs(remainder_x) > DIVISION_TOLERANCE)
    r_cutoff_x += remainder_x/(std::floor(domain_size_[0]/r_cutoff));
  
  if (std::abs(remainder_y) > DIVISION_TOLERANCE)
    r_cutoff_y += remainder_y/(std::floor(domain_size_[1]/r_cutoff));

  if (std::abs(remainder_z) > DIVISION_TOLERANCE)
    r_cutoff_z += remainder_z/(std::floor(domain_size_[1]/r_cutoff));

  // TODO rework this
  x = static_cast<size_t>(domain_size_[0] / r_cutoff);
  y = static_cast<size_t>(domain_size_[1] / r_cutoff);
  z = domain_size.size() == 3
          ? (static_cast<size_t>((domain_size_[2] / r_cutoff)))
          : 1;
  cells = std::vector<Cell>(x * y * z, Cell());
}

LinkedCellContainer::LinkedCellContainer()
    : domain_size_{0, 0, 0}, r_cutoff_{0}, left_corner_coordinates{0.0, 0.0, 0.0},
      x{0}, y{0}, z{0}, boundary_conditions_{} {}

void LinkedCellContainer::insert(Particle &p, bool placement) {
  ParticlePointer p_ptr = std::make_shared<Particle>(p);
  if (placement && is_within_domain(p_ptr->getX())) {
    std::array<double, 3> position = p.getX();
    size_t i = static_cast<size_t>((position[0] - left_corner_coordinates[0]) /
                                   r_cutoff_x);
    size_t j = static_cast<size_t>((position[1] - left_corner_coordinates[1]) /
                                   r_cutoff_y);
    size_t k = domain_size_.size() == 3
                   ? static_cast<size_t>(
                         (position[2] - left_corner_coordinates[2]) / r_cutoff_z)
                   : 0;
    size_t index = i + j * x + k * x * y;
    cells[index].particles.push_back(p_ptr);
  } else if(!is_within_domain(p_ptr->getX())) {
    p_ptr->left_domain = true;
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
    ParticlePointer p, const std::array<double, 3> &old_position) {
  if (is_within_domain(old_position)) {
    size_t i = static_cast<size_t>(
        (old_position[0] - left_corner_coordinates[0]) / r_cutoff_);
    size_t j = static_cast<size_t>(
        (old_position[1] - left_corner_coordinates[1]) / r_cutoff_);
    size_t k =
        domain_size_.size() == 3
            ? static_cast<size_t>(
                  (old_position[2] - left_corner_coordinates[2]) / r_cutoff_)
            : 0;
    size_t old_index = i + j * x + k * x * y;
    cells[old_index].particles.erase(
        std::remove(cells[old_index].particles.begin(),
                    cells[old_index].particles.end(), p),
        cells[old_index].particles.end());
  }
  if (is_within_domain(p->getX())) {
    size_t i = static_cast<size_t>((p->getX()[0] - left_corner_coordinates[0]) /
                                   r_cutoff_);
    size_t j = static_cast<size_t>((p->getX()[1] - left_corner_coordinates[1]) /
                                   r_cutoff_);
    size_t k =
        domain_size_.size() == 3
            ? static_cast<size_t>((p->getX()[2] - left_corner_coordinates[2]) /
                                  r_cutoff_)
            : 0;
    size_t index = i + j * x + k * x * y;
    cells[index].particles.push_back(p);
  } else {
    p->left_domain = true;
  }
}

Cell &LinkedCellContainer::get_cell(size_t index) { return cells[index]; }

std::vector<ParticlePointer> LinkedCellContainer::get_neighbours(Particle &p) {
  std::array<double, 3> position = p.getX();
  size_t i = static_cast<size_t>((position[0] - left_corner_coordinates[0]) /
                                 r_cutoff_);
  size_t j = static_cast<size_t>((position[1] - left_corner_coordinates[1]) /
                                 r_cutoff_);
  size_t k = domain_size_.size() == 3
                 ? static_cast<size_t>(
                       (position[2] - left_corner_coordinates[2]) / r_cutoff_)
                 : 0;

  int index = i + j * x + k * x * y;
  std::vector<ParticlePointer> neighbours(cells[index].particles);

  int X = static_cast<int>(x);
  int Y = static_cast<int>(y);
  int Z = static_cast<int>(z);
  // handle 2d case
  // Convert 1D index to 3D coordinates
  size_t i_ = index / (Y * Z);
  size_t j_ = (index % (Y * Z)) / Z;
  size_t k_ = index % Z;

  for (int di = -1; di <= 1; ++di) {
    for (int dj = -1; dj <= 1; ++dj) {
      for (int dk = -1; dk <= 1; ++dk) {
        if (di == 0 && dj == 0 && dk == 0) {
          continue;
        }
        int ni = i_ + di;
        int nj = j_ + dj;
        int nk = k_ + dk;

        if (ni >= 0 && ni < X && nj >= 0 && nj < Y && nk >= 0 && nk < Z) {
          int neighborIndex = ni * (Y * Z) + nj * Z + nk;
          neighbours.insert(neighbours.end(),
                            cells[neighborIndex].particles.begin(),
                            cells[neighborIndex].particles.end());
        }
      }
    }
  }

  return neighbours;
}

void LinkedCellContainer::clear() {
  for (size_t i = 0; i < cells.size(); ++i) {
    cells[i].particles.clear();
  }
  particles.clear();
}

// needs to find the leftmost and rightmost corner
void LinkedCellContainer::reinitialize(DirectSumContainer &container) {
  clear();
  auto current_low_left = container[0].getX();
  auto current_up_right = container[0].getX();
  for (auto &p : container) {
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
  for (auto &p : container) {
    insert(p, true);
  }
}

void LinkedCellContainer::reinitialize(std::vector<Particle> &particles) {
  clear();
  auto current_low_left = particles[0].getX();
  auto current_up_right = particles[0].getX();
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
  for (auto &p : particles) {
    insert(p, true);
  }
}

void LinkedCellContainer::reinitialize(
    std::vector<ParticlePointer> &particles) {
  clear();
  auto current_low_left = particles[0]->getX();
  auto current_up_right = particles[0]->getX();
  for (auto &p : particles) {
    if (p->getX()[0] < current_low_left[0] ||
        p->getX()[1] < current_low_left[1] ||
        p->getX()[2] < current_low_left[2]) {
      current_low_left = p->getX();
    }
    if (p->getX()[0] > current_up_right[0] ||
        p->getX()[1] > current_up_right[1] ||
        p->getX()[2] > current_up_right[2]) {
      current_up_right = p->getX();
    }
  }
  readjust_coordinates(current_low_left, current_up_right);
  for (auto &p : particles) {
    insert(*p, true);
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
  logger.debug("Current low left: " + std::to_string(current_low_left[0]) +
               " " + std::to_string(current_low_left[1]) + " " +
               std::to_string(current_low_left[2]));
  logger.debug("Current up right: " + std::to_string(current_up_right[0]) +
               " " + std::to_string(current_up_right[1]) + " " +
               std::to_string(current_up_right[2]));
  readjust_coordinates(current_low_left, current_up_right);
  auto particles_ = particles;
  clear();
  for (auto &p : particles_) {
    insert(p, true);
  }
  logger.debug("New low left: " + std::to_string(left_corner_coordinates[0]) +
               " " + std::to_string(left_corner_coordinates[1]) + " " +
               std::to_string(left_corner_coordinates[2]));
}

void LinkedCellContainer::handleBoundaryConditions(Particle &p) {
  auto position = p.getX();
  auto velocity = p.getV();

  // left boundary
  if (position[0] < left_corner_coordinates[0]) {
    if (boundary_conditions_.left == BoundaryCondition::Reflecting) {
      position[0] = 2 * left_corner_coordinates[0] - position[0];
      velocity[0] = -velocity[0];
    } else if (boundary_conditions_.left == BoundaryCondition::Outflow) {
      p.markForRemoval();
    }
  }

  // Right boundary
  if (position[0] > left_corner_coordinates[0] + domain_size_[0]) {
    if (boundary_conditions_.right == BoundaryCondition::Reflecting) {
      position[0] =
          2 * (left_corner_coordinates[0] + domain_size_[0]) - position[0];
      velocity[0] = -velocity[0];
    } else if (boundary_conditions_.right == BoundaryCondition::Outflow) {
      p.markForRemoval();
    }
  }

  // Bottom boundary
  if (position[1] < left_corner_coordinates[1]) {
    if (boundary_conditions_.bottom == BoundaryCondition::Reflecting) {
      position[1] = 2 * left_corner_coordinates[1] - position[1];
      velocity[1] = -velocity[1];
    } else if (boundary_conditions_.bottom == BoundaryCondition::Outflow) {
      p.markForRemoval();
    }
  }

  // Top boundary
  if (position[1] > left_corner_coordinates[1] + domain_size_[1]) {
    if (boundary_conditions_.top == BoundaryCondition::Reflecting) {
      position[1] =
          2 * (left_corner_coordinates[1] + domain_size_[1]) - position[1];
      velocity[1] = -velocity[1];
    } else if (boundary_conditions_.top == BoundaryCondition::Outflow) {
      p.markForRemoval();
    }
  }

  if (domain_size_.size() > 2) {
    // Front boundary
    if (position[2] < left_corner_coordinates[2]) {
      if (boundary_conditions_.front == BoundaryCondition::Reflecting) {
        position[2] = 2 * left_corner_coordinates[2] - position[2];
        velocity[2] = -velocity[2];
      } else if (boundary_conditions_.front == BoundaryCondition::Outflow) {
        p.markForRemoval();
      }
    }
    // Back boundary
    if (position[2] > left_corner_coordinates[2] + domain_size_[2]) {
      if (boundary_conditions_.back == BoundaryCondition::Reflecting) {
        position[2] =
            2 * (left_corner_coordinates[2] + domain_size_[2]) - position[2];
        velocity[2] = -velocity[2];
      } else if (boundary_conditions_.back == BoundaryCondition::Outflow) {
        p.markForRemoval();
      }
    }
  }

  // TODO optimize later
  p.updateX(position[0], position[1], position[2]);
  p.updateV(velocity[0], velocity[1], velocity[2]);
}

void LinkedCellContainer::removeOutflowParticles() {
  for (auto &cell : cells) {
    cell.particles.erase(std::remove_if(cell.particles.begin(),
                                        cell.particles.end(),
                                        [](const ParticlePointer &p) {
                                          return p->isMarkedForRemoval();
                                        }),
                         cell.particles.end());
  }
}

void LinkedCellContainer::updateParticles() {
  for (auto &cell : cells) {
    for (auto &p : cell.particles) {
      handleBoundaryConditions(*p);
    }
  }

  removeOutflowParticles();
}

size_t LinkedCellContainer::size() { return particles.size(); }

Particle &LinkedCellContainer::operator[](size_t index) {
  return particles[index];
}
