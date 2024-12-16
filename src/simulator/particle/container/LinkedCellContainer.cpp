#include "LinkedCellContainer.h"
#include "../../Simulation.h"
#include <cmath>

size_t LinkedCellContainer::Cell::size() const { return particle_ids.size(); }

void LinkedCellContainer::Cell::insert(int id) { particle_ids.insert(id); }

void LinkedCellContainer::Cell::remove(int id) { particle_ids.erase(id); }

LinkedCellContainer::LinkedCellContainer(
    std::initializer_list<double> domain_size, double r_cutoff,
    const DomainBoundaryConditions &boundary_conditions)
    : domain_size_{domain_size}, r_cutoff_{r_cutoff}, r_cutoff_x{r_cutoff},
      r_cutoff_y{r_cutoff}, r_cutoff_z{r_cutoff},
      left_corner_coordinates{0.0, 0.0, 0.0}, x{0}, y{0}, z{0},
      boundary_conditions_{boundary_conditions}, particles{}, cells_map{},
      particle_id{0}, particles_left_domain{0}, is_wrapper{false},
      halo_count{0}, reflective_flag{false}, periodic_flag{false} {
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

  mark_halo_cells();
}

LinkedCellContainer::LinkedCellContainer()
    : domain_size_{0, 0, 0}, r_cutoff_{0}, left_corner_coordinates{0.0, 0.0,
                                                                   0.0},
      x{0}, y{0}, z{0}, boundary_conditions_{}, cells_map{}, particle_id{0},
      particles_left_domain{0}, is_wrapper{false}, halo_count{0},
      reflective_flag{false}, periodic_flag{false} {}

void LinkedCellContainer::insert(Particle &p, bool placement) {
  ParticlePointer p_ptr = std::make_shared<Particle>(p);
  if (placement && is_within_domain(p_ptr->getX())) {
    std::array<double, 3> position = p.getX();
    size_t i = static_cast<size_t>((position[0] - left_corner_coordinates[0]) /
                                   r_cutoff_x);
    size_t j = static_cast<size_t>((position[1] - left_corner_coordinates[1]) /
                                   r_cutoff_y);
    size_t k =
        domain_size_.size() == 3
            ? static_cast<size_t>((position[2] - left_corner_coordinates[2]) /
                                  r_cutoff_z)
            : 0;
    size_t index = i + j * x + k * x * y;
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

  // compute previous index
  size_t i = static_cast<size_t>(
      (old_position[0] - left_corner_coordinates[0]) / r_cutoff_x);
  size_t j = static_cast<size_t>(
      (old_position[1] - left_corner_coordinates[1]) / r_cutoff_y);
  size_t k =
      domain_size_.size() == 3
          ? static_cast<size_t>((old_position[2] - left_corner_coordinates[2]) /
                                r_cutoff_z)
          : 0;
  size_t old_index = i + j * x + k * x * y;

  size_t q = static_cast<size_t>(
      (cells_map[particle_id]->getX()[0] - left_corner_coordinates[0]) /
      r_cutoff_x);
  size_t v = static_cast<size_t>(
      (cells_map[particle_id]->getX()[1] - left_corner_coordinates[1]) /
      r_cutoff_y);
  size_t w = domain_size_.size() == 3
                 ? static_cast<size_t>((cells_map[particle_id]->getX()[2] -
                                        left_corner_coordinates[2]) /
                                       r_cutoff_z)
                 : 0;
  size_t current_index = q + v * x + w * x * y;
  // compute current index
  if (current_index != old_index) {

    if (is_within_domain(old_position)) {
      cells[old_index].remove(cells_map[particle_id]->getType());
    }
    if (is_within_domain(cells_map[particle_id]->getX())) {
      cells[current_index].insert(cells_map[particle_id]->getType());
      if(cells_map[particle_id]->left_domain){
        cells_map[particle_id]->left_domain = false;
        particles_left_domain--;
      }
      if(cells_map[particle_id]->is_periodic_copy && cells_map[particle_id]->secondary_copy_flag){
        cells_map[particle_id]->is_periodic_copy = false;
        particles_left_domain--;
      }

      if (reflective_flag) {
        handle_boundary_conditions(particle_id, current_index);
      } else if(periodic_flag){
        handle_periodic_boundary_conditions(particle_id, current_index);
      }
    } else {
      if(cells_map[particle_id]->left_domain == false){
        cells_map[particle_id]->left_domain = true;
        particles_left_domain++;
      }
    }
  }
}

LinkedCellContainer::Cell &LinkedCellContainer::get_cell(size_t index) {
  return cells[index];
}

std::vector<ParticlePointer>
LinkedCellContainer::get_neighbours(int particle_id) {
  std::vector<ParticlePointer> neighbours{};
  if (cells_map[particle_id]->left_domain) {
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

  // logger.info("Current cell index: " + std::to_string(index));
  for (auto &i : cells[index].particle_ids) {
    neighbours.push_back(cells_map[i]);
  }

  // logger.info("X: " + std::to_string(x));
  // logger.info("Y: " + std::to_string(y));
  // logger.info("Z: " + std::to_string(z));

  // logger.info("Current cell index: " + std::to_string(i) + " " +
  // std::to_string(j) + " " + std::to_string(k));

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

          // logger.info("Neighbour cell index: " +
          // std::to_string(neighborIndex));
          for (auto &s : cells[neighborIndex].particle_ids) {
            neighbours.push_back(cells_map[s]);
          }
        }
      }
    }
  }
  // logger.info("Neighbours size: " + std::to_string(neighbours.size()));
  // logger.info("Cells size" + std::to_string(cells.size()));

  /*for(auto &n : neighbours){
    if(n -> getType() != particle_id)
      logger.info("Neighbour: " + n->toString());
  }
  for(auto &n : particles.get_all_particles()){
    if(n->getType() != particle_id)
      logger.info("Particle: " + n->toString());
  }*/

  // logger.info("Checkpoint");
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
  readjust_coordinates(current_low_left, current_up_right);
  auto particles_ = particles;
  clear();
  for (auto &p : particles_) {
    insert(p, true);
  }
}

void LinkedCellContainer::handle_boundary_conditions(int particle_id, int cell_index) {
  //ensure that the particle is in halo cell
  if(!cells[cell_index].is_halo){
    return;
  }

  auto& velocity = cells_map[particle_id]->getV();
  logger.debug("Checking");
  // Left boundary
  if (cell_index % x == 0) {
    if (boundary_conditions_.left == BoundaryCondition::Reflecting) {
      logger.debug("Reflecting left boundary");
      cells_map[particle_id]->updateV(-velocity[0], velocity[1], velocity[2]);
      //logger.debug("Velocity: " + std::to_string(cells_map[particle_id]->getV()[0]) + " " + std::to_string(cells_map[particle_id]->getV()[1]) + " " + std::to_string(cells_map[particle_id]->getV()[2]));
    }
  }

  // Right boundary
  if ((cell_index + 1) % x == 0) { 
    if (boundary_conditions_.right == BoundaryCondition::Reflecting) {
      logger.debug("Reflecting right boundary");
      cells_map[particle_id]->updateV(-velocity[0], velocity[1], velocity[2]);
    }
  }

  // Bottom boundary
  for(size_t i = 0; i < z; i++){
    if (x * y * i <= cell_index && cell_index < x * y * i + x) {
      if (boundary_conditions_.bottom == BoundaryCondition::Reflecting) {
        logger.debug("Reflecting bottom boundary");
        logger.debug("Particle location: " + std::to_string(cells_map[particle_id]->getX()[0]) + " " + std::to_string(cells_map[particle_id]->getX()[1]) + " " + std::to_string(cells_map[particle_id]->getX()[2]));  
        logger.debug("Velocity: " + std::to_string(cells_map[particle_id]->getV()[0]) + " " + std::to_string(cells_map[particle_id]->getV()[1]) + " " + std::to_string(cells_map[particle_id]->getV()[2]));
        cells_map[particle_id]->updateV(velocity[0], -velocity[1], velocity[2]);
        logger.debug("Velocity: " + std::to_string(cells_map[particle_id]->getV()[0]) + " " + std::to_string(cells_map[particle_id]->getV()[1]) + " " + std::to_string(cells_map[particle_id]->getV()[2]));

      }
    }
  }

  // Top boundary
  for(size_t i = 1; i <= z; i++){
    if (x * y * i - x <= cell_index && cell_index < x * y * i) {
      if (boundary_conditions_.top == BoundaryCondition::Reflecting) {
        logger.debug("Reflecting top boundary");
        cells_map[particle_id]->updateV(velocity[0], -velocity[1], velocity[2]);
      }
    }
  }

  if (z > 1) {
    // Front boundary
    if (cell_index < x * y) {
      if (boundary_conditions_.front == BoundaryCondition::Reflecting) {
        logger.debug("Reflecting front boundary");
        cells_map[particle_id]->updateV(velocity[0], velocity[1], -velocity[2]);
      }
    }
    // Back boundary
    if (cell_index >= x * y * (z - 1)) {
      if (boundary_conditions_.back == BoundaryCondition::Reflecting) {
        logger.debug("Reflecting back boundary");
        cells_map[particle_id]->updateV(velocity[0], velocity[1], -velocity[2]);
      }
    }
  }
}

void LinkedCellContainer::handle_periodic_boundary_conditions(int particle_id,
                                                              int cell_index) {

  if(cells_map[particle_id]->secondary_copy_flag && !cells[cell_index].is_halo){
    cells_map[particle_id]->secondary_copy_flag = false;
    return;
  }

  if(!cells[cell_index].is_halo || cells_map[particle_id]->secondary_copy_flag && cells[cell_index].is_halo){ 
    return;
  }

  logger.info("Handling periodic boundary conditions");

  //this is way too long and needs to be refactored
  if (z == 1) {
    logger.info("2D periodic boundary conditions");
    // handle corner case
    if (cell_index == 0) {
      logger.info("Bottom left corner");
      cells_map[particle_id]->updateX(
          cells_map[particle_id]->getX()[0] + domain_size_[0],
          cells_map[particle_id]->getX()[1] + domain_size_[1],
          cells_map[particle_id]->getX()[2]);
      cells_map[particle_id]->left_domain = true;
      particle_id++;
      Particle p_1{std::array<double, 3>{cells_map[particle_id]->getX()[0] +
                                             domain_size_[0],
                                         cells_map[particle_id]->getX()[1],
                                         cells_map[particle_id]->getX()[2]},
                   std::array<double, 3>{cells_map[particle_id]->getV()[0],
                                         -cells_map[particle_id]->getV()[1],
                                         cells_map[particle_id]->getV()[2]},
                   cells_map[particle_id]->getM(),
                   particle_id,
                   cells_map[particle_id]->getEpsilon(),
                   cells_map[particle_id]->getSigma()};
      particle_id++;
      Particle p_2{std::array<double, 3>{cells_map[particle_id]->getX()[0],
                                         cells_map[particle_id]->getX()[1] +
                                             domain_size_[1],
                                         cells_map[particle_id]->getX()[2]},
                   std::array<double, 3>{-cells_map[particle_id]->getV()[0],
                                         cells_map[particle_id]->getV()[1],
                                         cells_map[particle_id]->getV()[2]},
                   cells_map[particle_id]->getM(),
                   particle_id,
                   cells_map[particle_id]->getEpsilon(),
                   cells_map[particle_id]->getSigma()};
      insert(p_1, false);
      insert(p_2, false);

    } else if (cell_index == x - 1) {
      logger.info("Bottom right corner");
      cells_map[particle_id]->updateX(
          cells_map[particle_id]->getX()[0] - domain_size_[0],
          cells_map[particle_id]->getX()[1] + domain_size_[1],
          cells_map[particle_id]->getX()[2]);
      cells_map[particle_id]->left_domain = true;
      particle_id++;
      Particle p_1{std::array<double, 3>{cells_map[particle_id]->getX()[0] -
                                             domain_size_[0],
                                         cells_map[particle_id]->getX()[1],
                                         cells_map[particle_id]->getX()[2]},
                   std::array<double, 3>{cells_map[particle_id]->getV()[0],
                                         -cells_map[particle_id]->getV()[1],
                                         cells_map[particle_id]->getV()[2]},
                   cells_map[particle_id]->getM(),
                   particle_id,
                   cells_map[particle_id]->getEpsilon(),
                   cells_map[particle_id]->getSigma()};
      particle_id++;
      Particle p_2{std::array<double, 3>{cells_map[particle_id]->getX()[0],
                                         cells_map[particle_id]->getX()[1] +
                                             domain_size_[1],
                                         cells_map[particle_id]->getX()[2]},
                   std::array<double, 3>{-cells_map[particle_id]->getV()[0],
                                         cells_map[particle_id]->getV()[1],
                                         cells_map[particle_id]->getV()[2]},
                   cells_map[particle_id]->getM(),
                   particle_id,
                   cells_map[particle_id]->getEpsilon(),
                   cells_map[particle_id]->getSigma()};
      insert(p_1, false);
      insert(p_2, false);

    } else if (cell_index == x * y - 1) {
      logger.info("Top right corner");
      cells_map[particle_id]->updateX(
          cells_map[particle_id]->getX()[0] - domain_size_[0],
          cells_map[particle_id]->getX()[1] - domain_size_[1],
          cells_map[particle_id]->getX()[2]);
      cells_map[particle_id]->left_domain = true;

      particle_id++;
      Particle p_1{std::array<double, 3>{cells_map[particle_id]->getX()[0] -
                                             domain_size_[0],
                                         cells_map[particle_id]->getX()[1],
                                         cells_map[particle_id]->getX()[2]},
                   std::array<double, 3>{-cells_map[particle_id]->getV()[0],
                                         cells_map[particle_id]->getV()[1],
                                         cells_map[particle_id]->getV()[2]},
                   cells_map[particle_id]->getM(),
                   particle_id,
                   cells_map[particle_id]->getEpsilon(),
                   cells_map[particle_id]->getSigma()};
      particle_id++;
      Particle p_2{std::array<double, 3>{cells_map[particle_id]->getX()[0],
                                         cells_map[particle_id]->getX()[1] -
                                             domain_size_[1],
                                         cells_map[particle_id]->getX()[2]},
                   std::array<double, 3>{cells_map[particle_id]->getV()[0],
                                         -cells_map[particle_id]->getV()[1],
                                         cells_map[particle_id]->getV()[2]},
                   cells_map[particle_id]->getM(),
                   particle_id,
                   cells_map[particle_id]->getEpsilon(),
                   cells_map[particle_id]->getSigma()};
      insert(p_1, false);
      insert(p_2, false);

    } else if (cell_index == x * y - x) {
      logger.info("Top left corner");
      cells_map[particle_id]->updateX(
          cells_map[particle_id]->getX()[0] + domain_size_[0],
          cells_map[particle_id]->getX()[1] - domain_size_[1],
          cells_map[particle_id]->getX()[2]);
      cells_map[particle_id]->left_domain = true;

      particle_id++;
      Particle p_1{std::array<double, 3>{cells_map[particle_id]->getX()[0] +
                                             domain_size_[0],
                                         cells_map[particle_id]->getX()[1],
                                         cells_map[particle_id]->getX()[2]},
                   std::array<double, 3>{-cells_map[particle_id]->getV()[0],
                                         cells_map[particle_id]->getV()[1],
                                         cells_map[particle_id]->getV()[2]},
                   cells_map[particle_id]->getM(),
                   particle_id,
                   cells_map[particle_id]->getEpsilon(),
                   cells_map[particle_id]->getSigma()};
      particle_id++;
      Particle p_2{std::array<double, 3>{cells_map[particle_id]->getX()[0],
                                         cells_map[particle_id]->getX()[1] -
                                             domain_size_[1],
                                         cells_map[particle_id]->getX()[2]},
                   std::array<double, 3>{cells_map[particle_id]->getV()[0],
                                         -cells_map[particle_id]->getV()[1],
                                         cells_map[particle_id]->getV()[2]},
                   cells_map[particle_id]->getM(),
                   particle_id,
                   cells_map[particle_id]->getEpsilon(),
                   cells_map[particle_id]->getSigma()};
      insert(p_1, false);
      insert(p_2, false);
    } else if( cell_index % x == 0){
      logger.info("Left boundary");
      cells_map[particle_id]->updateX(
          cells_map[particle_id]->getX()[0] + domain_size_[0],
          cells_map[particle_id]->getX()[1],
          cells_map[particle_id]->getX()[2]);
      cells_map[particle_id]->left_domain = true;
    } else if( (cell_index + 1) % x == 0){
      logger.info("Right boundary");
      cells_map[particle_id]->updateX(
          cells_map[particle_id]->getX()[0] - domain_size_[0],
          cells_map[particle_id]->getX()[1],
          cells_map[particle_id]->getX()[2]);
      cells_map[particle_id]->left_domain = true;
    } else if( cell_index < x){
      logger.info("Bottom boundary");
      cells_map[particle_id]->updateX(
          cells_map[particle_id]->getX()[0],
          cells_map[particle_id]->getX()[1] + domain_size_[1],
          cells_map[particle_id]->getX()[2]);
      cells_map[particle_id]->is_periodic_copy = true;
      cells_map[particle_id]->secondary_copy_flag = true;
      particles_left_domain++;
    } else if( cell_index >= x * (y - 1)){
      logger.info("Top boundary");
      cells_map[particle_id]->updateX(
          cells_map[particle_id]->getX()[0],
          cells_map[particle_id]->getX()[1] - domain_size_[1],
          cells_map[particle_id]->getX()[2]);
      cells_map[particle_id]->is_periodic_copy = true;
      cells_map[particle_id]->secondary_copy_flag = true;
      particles_left_domain++;
    }
  }
  // determine if corner cel
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
          halo_count++;
        }
      }
    }
  }
}
