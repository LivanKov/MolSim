#include "LinkedCellContainer.h"
#include <cmath>

void LinkedCellContainer::insert(Particle &p) {
  std::array<double, 3> position = p.getX();
  size_t i = static_cast<size_t>((position[0] - left_corner_coordinates[0]) /
                                 r_cutoff_);
  size_t j = static_cast<size_t>((position[1] - left_corner_coordinates[1]) /
                                 r_cutoff_);
  size_t k = domain_size_.size() == 3
                 ? static_cast<size_t>(
                       (position[2] - left_corner_coordinates[2]) / r_cutoff_)
                 : 0;
  size_t index = i + j * x + k * x * y;
  unwrapped_cells_[index].particles.push_back(std::make_shared<Particle>(p));
}

void LinkedCellContainer::update_particle_location(
    Particle &p, std::array<double, 3> &old_position) {
  size_t i = static_cast<size_t>(
      (old_position[0] - left_corner_coordinates[0]) / r_cutoff_);
  size_t j = static_cast<size_t>(
      (old_position[1] - left_corner_coordinates[1]) / r_cutoff_);
  size_t k =
      domain_size_.size() == 3
          ? static_cast<size_t>((old_position[2] - left_corner_coordinates[2]) /
                                r_cutoff_)
          : 0;
  size_t old_index = i + j * x + k * x * y;
  unwrapped_cells_[old_index].particles.erase(
      std::remove_if(
          unwrapped_cells_[old_index].particles.begin(),
          unwrapped_cells_[old_index].particles.end(),
          [&p](const ParticlePointer &particle) { return *particle == p; }),
      unwrapped_cells_[old_index].particles.end());
  insert(p);
}

Cell &LinkedCellContainer::get_cell(size_t index) {
  return unwrapped_cells_[index];
}

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
  size_t index = i + j * x + k * x * y;
  std::vector<ParticlePointer> neighbours;

  // handle 2d case
  if (domain_size_.size() == 2) {
    if (i == x - 1 && j == y - 1) {
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - 1].particles.begin(),
                        unwrapped_cells_[index - 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x].particles.begin(),
                        unwrapped_cells_[index - x].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x - 1].particles.begin(),
                        unwrapped_cells_[index - x - 1].particles.end());
    } else if (i == x - 1 && j == 0) {
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - 1].particles.begin(),
                        unwrapped_cells_[index - 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + x].particles.begin(),
                        unwrapped_cells_[index + x].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + x - 1].particles.begin(),
                        unwrapped_cells_[index + x - 1].particles.end());
    } else if (i == 0 && j == y - 1) {
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + 1].particles.begin(),
                        unwrapped_cells_[index + 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x].particles.begin(),
                        unwrapped_cells_[index - x].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x + 1].particles.begin(),
                        unwrapped_cells_[index - x + 1].particles.end());
    } else if (i == 0 && j == 0) {
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + 1].particles.begin(),
                        unwrapped_cells_[index + 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + x].particles.begin(),
                        unwrapped_cells_[index + x].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + x + 1].particles.begin(),
                        unwrapped_cells_[index + x + 1].particles.end());
    } else if (i == x - 1) {
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - 1].particles.begin(),
                        unwrapped_cells_[index - 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + x].particles.begin(),
                        unwrapped_cells_[index + x].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x].particles.begin(),
                        unwrapped_cells_[index - x].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + x - 1].particles.begin(),
                        unwrapped_cells_[index + x - 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x - 1].particles.begin(),
                        unwrapped_cells_[index - x - 1].particles.end());
    } else if (j == y - 1) {
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + 1].particles.begin(),
                        unwrapped_cells_[index + 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x].particles.begin(),
                        unwrapped_cells_[index - x].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - 1].particles.begin(),
                        unwrapped_cells_[index - 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x + 1].particles.begin(),
                        unwrapped_cells_[index - x + 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x - 1].particles.begin(),
                        unwrapped_cells_[index - x - 1].particles.end());
    } else if (i == 0) {
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + 1].particles.begin(),
                        unwrapped_cells_[index + 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + x].particles.begin(),
                        unwrapped_cells_[index + x].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x].particles.begin(),
                        unwrapped_cells_[index - x].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + x + 1].particles.begin(),
                        unwrapped_cells_[index + x + 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x + 1].particles.begin(),
                        unwrapped_cells_[index - x + 1].particles.end());
    } else if (j == 0) {
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + 1].particles.begin(),
                        unwrapped_cells_[index + 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x].particles.begin(),
                        unwrapped_cells_[index - x].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - 1].particles.begin(),
                        unwrapped_cells_[index - 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + x - 1].particles.begin(),
                        unwrapped_cells_[index + x - 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x - 1].particles.begin(),
                        unwrapped_cells_[index - x - 1].particles.end());
    } else {
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + 1].particles.begin(),
                        unwrapped_cells_[index + 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - 1].particles.begin(),
                        unwrapped_cells_[index - 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + x].particles.begin(),
                        unwrapped_cells_[index + x].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x].particles.begin(),
                        unwrapped_cells_[index - x].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + x + 1].particles.begin(),
                        unwrapped_cells_[index + x + 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index + x - 1].particles.begin(),
                        unwrapped_cells_[index + x - 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x + 1].particles.begin(),
                        unwrapped_cells_[index - x + 1].particles.end());
      neighbours.insert(neighbours.end(),
                        unwrapped_cells_[index - x - 1].particles.begin(),
                        unwrapped_cells_[index - x - 1].particles.end());
    }
  } else {
  }

  return neighbours;
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

  p.updateX(position[0], position[1], position[2]);
  p.updateV(velocity[0], velocity[1], velocity[2]);
}

void LinkedCellContainer::removeOutflowParticles() {
  for (auto &cell : unwrapped_cells_) {
    cell.particles.erase(std::remove_if(cell.particles.begin(),
                                        cell.particles.end(),
                                        [](const ParticlePointer &p) {
                                          return p->isMarkedForRemoval();
                                        }),
                         cell.particles.end());
  }
}

void LinkedCellContainer::updateParticles() {
    for (auto &cell : unwrapped_cells_) {
        for (auto &p : cell.particles) {
            handleBoundaryConditions(*p);
    
        }
    }

    removeOutflowParticles();
}