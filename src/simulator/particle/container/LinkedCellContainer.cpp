#include "LinkedCellContainer.h"
#include <cmath>
#include "../../Simulation.h"

size_t LinkedCellContainer::Cell::size() const { return particle_ids.size(); }

void LinkedCellContainer::Cell::insert(int id) { particle_ids.insert(id); }

void LinkedCellContainer::Cell::remove(int id) { particle_ids.erase(id); }

LinkedCellContainer::LinkedCellContainer(
    std::initializer_list<double> domain_size, double r_cutoff,
    const DomainBoundaryConditions &boundary_conditions)
    : domain_size_{domain_size}, r_cutoff_{r_cutoff},
      r_cutoff_x{r_cutoff}, r_cutoff_y{r_cutoff}, r_cutoff_z{r_cutoff},
      left_corner_coordinates{0.0, 0.0, 0.0}, x{0}, y{0}, z{0},
      boundary_conditions_{boundary_conditions}, particles{}, cells_map{}, particle_id{0}, particles_left_domain{0}, is_wrapper{false} {
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

  if(domain_size.size() == 2){
    domain_size_.push_back(r_cutoff_z);
  }

  // TODO rework this
  x = static_cast<size_t>(domain_size_[0] / r_cutoff);
  y = static_cast<size_t>(domain_size_[1] / r_cutoff);
  z = domain_size.size() == 3
          ? (static_cast<size_t>((domain_size_[2] / r_cutoff)))
          : 1;
  cells = std::vector<Cell>(x * y * z, Cell());

  for(auto& p : particles.get_all_particles()){
    cells_map[p->getType()] = p;
  }

}

LinkedCellContainer::LinkedCellContainer()
    : domain_size_{0, 0, 0}, r_cutoff_{0}, left_corner_coordinates{0.0, 0.0, 0.0},
      x{0}, y{0}, z{0}, boundary_conditions_{}, cells_map{}, particle_id{0}, particles_left_domain{0}, is_wrapper{false} {}

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
    cells[index].insert(p_ptr->getType());
  } else if(!is_within_domain(p_ptr->getX())) {
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

    //compute previous index
      size_t i = static_cast<size_t>(
        (old_position[0] - left_corner_coordinates[0]) / r_cutoff_x);
    size_t j = static_cast<size_t>(
        (old_position[1] - left_corner_coordinates[1]) / r_cutoff_y);
    size_t k =
        domain_size_.size() == 3
            ? static_cast<size_t>(
                  (old_position[2] - left_corner_coordinates[2]) / r_cutoff_z)
            : 0;
    size_t old_index = i + j * x + k * x * y;


    size_t q = static_cast<size_t>((cells_map[particle_id]->getX()[0] - left_corner_coordinates[0]) /
                                   r_cutoff_x);
    size_t v = static_cast<size_t>((cells_map[particle_id]->getX()[1] - left_corner_coordinates[1]) /
                                   r_cutoff_y);
    size_t w =
        domain_size_.size() == 3
            ? static_cast<size_t>((cells_map[particle_id]->getX()[2] - left_corner_coordinates[2]) /
                                  r_cutoff_z)
            : 0;
    size_t current_index = q + v * x + w * x * y;
    //compute current index
  if(current_index != old_index){
    
  if (is_within_domain(old_position)) {
     cells[old_index].remove(cells_map[particle_id]->getType());
  } 
  if (is_within_domain(cells_map[particle_id]->getX())) {
    cells[current_index].insert(cells_map[particle_id]->getType());
  } else {
    cells_map[particle_id]->left_domain = true;
    particles_left_domain++;
  }
  }
}

LinkedCellContainer::Cell &LinkedCellContainer::get_cell(size_t index) { return cells[index]; }

std::vector<ParticlePointer> LinkedCellContainer::get_neighbours(int particle_id) {
  /*logger.info("Amount of Cells: " + std::to_string(cells.size()));
  logger.info("Available particles: " + std::to_string(particles.size()));
  logger.info("left domain particles: " + std::to_string(particles_left_domain));
  logger.info("Current particle id: " + std::to_string(particle_id));*/


  /*for(size_t i = 0; i < cells.size(); i++){
    if(cells[i].size() > 0){
      logger.info("New cell: "+ std::to_string(i));
      for(auto& p : cells[i].particle_ids){
        logger.info("Particle ID: " + std::to_string(p));
      }
    }
  }*/
  std::vector<ParticlePointer> neighbours{};
  if(cells_map[particle_id]->left_domain){
    return neighbours;
  }
  std::array<double, 3> position = cells_map[particle_id]->getX();
  size_t i = static_cast<size_t>((position[0] - left_corner_coordinates[0]) /
                                 r_cutoff_x);                         
  size_t j = static_cast<size_t>((position[1] - left_corner_coordinates[1]) /
                                 r_cutoff_y);
  size_t k = domain_size_.size() == 3
                 ? static_cast<size_t>(
                       (position[2] - left_corner_coordinates[2]) / r_cutoff_z)
                 : 0;

  int index = i + j * x + k * x * y;

  //logger.info("Current cell index: " + std::to_string(index));
  for(auto &i : cells[index].particle_ids){
      neighbours.push_back(cells_map[i]);
  }

  //logger.info("X: " + std::to_string(x));
  //logger.info("Y: " + std::to_string(y));
  //logger.info("Z: " + std::to_string(z));


  //logger.info("Current cell index: " + std::to_string(i) + " " + std::to_string(j) + " " + std::to_string(k));

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

          //logger.info("Neighbour cell index: " + std::to_string(neighborIndex));
          for(auto &s : cells[neighborIndex].particle_ids){
            neighbours.push_back(cells_map[s]);
          }
        }
      }
    }
  }
  //logger.info("Neighbours size: " + std::to_string(neighbours.size()));
  //logger.info("Cells size" + std::to_string(cells.size()));
  
  /*for(auto &n : neighbours){
    if(n -> getType() != particle_id)
      logger.info("Neighbour: " + n->toString());
  }
  for(auto &n : particles.get_all_particles()){
    if(n->getType() != particle_id)
      logger.info("Particle: " + n->toString());
  }*/

  //logger.info("Checkpoint");
  return neighbours;
}

void LinkedCellContainer::clear() {
  for (size_t i = 0; i < cells.size(); ++i) {
    cells[i].particle_ids.clear();
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
      p.left_domain = true;
      particles_left_domain++;
    }
  }

  // Right boundary
  if (position[0] > left_corner_coordinates[0] + domain_size_[0]) {
    if (boundary_conditions_.right == BoundaryCondition::Reflecting) {
      position[0] =
          2 * (left_corner_coordinates[0] + domain_size_[0]) - position[0];
      velocity[0] = -velocity[0];
    } else if (boundary_conditions_.right == BoundaryCondition::Outflow) {
      p.left_domain = true;
      particles_left_domain++;
    }
  }

  // Bottom boundary
  if (position[1] < left_corner_coordinates[1]) {
    if (boundary_conditions_.bottom == BoundaryCondition::Reflecting) {
      position[1] = 2 * left_corner_coordinates[1] - position[1];
      velocity[1] = -velocity[1];
    } else if (boundary_conditions_.bottom == BoundaryCondition::Outflow) {
      p.left_domain = true;
      particles_left_domain++;
    }
  }

  // Top boundary
  if (position[1] > left_corner_coordinates[1] + domain_size_[1]) {
    if (boundary_conditions_.top == BoundaryCondition::Reflecting) {
      position[1] =
          2 * (left_corner_coordinates[1] + domain_size_[1]) - position[1];
      velocity[1] = -velocity[1];
    } else if (boundary_conditions_.top == BoundaryCondition::Outflow) {
      p.left_domain = true;
      particles_left_domain++;
    }
  }

  if (domain_size_.size() > 2) {
    // Front boundary
    if (position[2] < left_corner_coordinates[2]) {
      if (boundary_conditions_.front == BoundaryCondition::Reflecting) {
        position[2] = 2 * left_corner_coordinates[2] - position[2];
        velocity[2] = -velocity[2];
      } else if (boundary_conditions_.front == BoundaryCondition::Outflow) {
        p.left_domain = true;
        particles_left_domain++;
      }
    }
    // Back boundary
    if (position[2] > left_corner_coordinates[2] + domain_size_[2]) {
      if (boundary_conditions_.back == BoundaryCondition::Reflecting) {
        position[2] =
            2 * (left_corner_coordinates[2] + domain_size_[2]) - position[2];
        velocity[2] = -velocity[2];
      } else if (boundary_conditions_.back == BoundaryCondition::Outflow) {
        p.left_domain = true;
        particles_left_domain++;
      }
    }
  }

  // TODO optimize later
  p.updateX(position[0], position[1], position[2]);
  p.updateV(velocity[0], velocity[1], velocity[2]);
}

void LinkedCellContainer::removeOutflowParticles() {
  /*for (auto &cell : cells) {
    cell.particles.erase(std::remove_if(cell.particles.begin(),
                                        cell.particles.end(),
                                        [](const ParticlePointer &p) {
                                          return p->left_domain;
                                        }),
                         cell.particles.end());
  }*/
}

void LinkedCellContainer::updateParticles() {
  /*for (auto &cell : cells) {
    for (auto &p : cell.particles) {
      handleBoundaryConditions(*p);
    }
  }

  removeOutflowParticles();*/
}

size_t LinkedCellContainer::size() { return particles.size(); }

Particle &LinkedCellContainer::operator[](size_t index) {
  return particles[index];
}
