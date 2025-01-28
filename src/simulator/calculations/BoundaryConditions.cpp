#include "BoundaryConditions.h"
#include "../particle/container/LinkedCellContainer.h"
#include <iostream>

void BoundaryConditions::run(LinkedCellContainer &particles) {
  particles.clear_ghost_particles();
  for (auto &cell_index : particles.halo_cell_indices) {
    for (auto &particle_id : particles.cells[cell_index].particle_ids) {
      auto position = particles.cells[cell_index].placement;
      if (particles.placement_map[position] == BoundaryCondition::Reflecting) {
        handle_reflect_conditions(particle_id, cell_index, particles);
      } else if (particles.placement_map[position] ==
                 BoundaryCondition::Periodic) {
        particles.create_ghost_particles(particle_id, cell_index);
      }
    }

    for (auto &particle_id : particles.particles_outbound) {
      auto &particle = particles.cells_map[particle_id];
      if (!particle->left_domain && !particle->outbound) {
        auto cell_index = particles.get_cell_index(particle->getOldX());
        auto position = particles.cells[cell_index].placement;
        particle->outbound = true;
        if (particles.placement_map[position] == BoundaryCondition::Outflow) {
          handle_outflow_conditions(particle_id, cell_index, particles);
        } else if (particles.placement_map[position] ==
                   BoundaryCondition::Periodic) {
          handle_periodic_conditions(particle_id, cell_index, particles);
        }
      }
    }
  }
}

void BoundaryConditions::handle_reflect_conditions(
    int particle_id, int cell_index, LinkedCellContainer &particles) {
  auto &velocity = particles.cells_map[particle_id]->getV();
  auto &cell = particles.cells[cell_index];
  // Left boundary or Right boundary
  if (cell.placement == Placement::LEFT || cell.placement == Placement::RIGHT)
    particles.cells_map[particle_id]->updateV(-velocity[0], velocity[1],
                                              velocity[2]);

  // Bottom boundary
  if (cell.placement == Placement::BOTTOM || cell.placement == Placement::TOP)
    particles.cells_map[particle_id]->updateV(velocity[0], -velocity[1],
                                              velocity[2]);

  // Front or Back
  if (particles.z > 1 &&
      (cell.placement == Placement::FRONT || cell.placement == Placement::BACK))
    particles.cells_map[particle_id]->updateV(velocity[0], velocity[1],
                                              -velocity[2]);

  if (particles.domain_size_.size() == 2) {
    // handle corners
    if (cell.placement == Placement::BOTTOM_LEFT_CORNER ||
        cell.placement == Placement::TOP_RIGHT_CORNER ||
        cell.placement == Placement::TOP_LEFT_CORNER ||
        cell.placement == Placement::BOTTOM_RIGHT_CORNER)
      particles.cells_map[particle_id]->updateV(-velocity[0], -velocity[1],
                                                velocity[2]);
  } else if (particles.domain_size_.size() == 3) {
    if (cell.placement == Placement::TOP_BACK_LEFT_CORNER ||
        cell.placement == Placement::TOP_BACK_RIGHT_CORNER ||
        cell.placement == Placement::TOP_FRONT_LEFT_CORNER ||
        cell.placement == Placement::TOP_FRONT_RIGHT_CORNER ||
        cell.placement == Placement::BOTTOM_BACK_LEFT_CORNER ||
        cell.placement == Placement::BOTTOM_BACK_RIGHT_CORNER ||
        cell.placement == Placement::BOTTOM_FRONT_LEFT_CORNER ||
        cell.placement == Placement::BOTTOM_FRONT_RIGHT_CORNER ||
        cell.placement == Placement::BOTTOM_BACK_EDGE ||
        cell.placement == Placement::BOTTOM_FRONT_EDGE ||
        cell.placement == Placement::BOTTOM_LEFT_EDGE ||
        cell.placement == Placement::BOTTOM_RIGHT_EDGE ||
        cell.placement == Placement::TOP_BACK_EDGE ||
        cell.placement == Placement::TOP_FRONT_EDGE ||
        cell.placement == Placement::TOP_RIGHT_EDGE ||
        cell.placement == Placement::TOP_LEFT_EDGE) {
      particles.cells_map[particle_id]->updateV(velocity[0], -velocity[1],
                                                velocity[2]);
    }
  }
}

void BoundaryConditions::handle_periodic_conditions(
    int particle_id, int cell_index, LinkedCellContainer &particles) {
  particles.cells_map[particle_id]->outbound = false;
  particles.particles_outbound.erase(
      std::remove(particles.particles_outbound.begin(),
                  particles.particles_outbound.end(), particle_id),
      particles.particles_outbound.end());

  auto old_x = particles.cells_map[particle_id]->getOldX();

  std::array<double, 3> location = particles.cells_map[particle_id]->getX();

  if (location[0] < particles.left_corner_coordinates[0])
    particles.cells_map[particle_id]->updateX(
        location[0] + particles.domain_size_[0], location[1], location[2]);
  if (location[0] >
      particles.left_corner_coordinates[0] + particles.domain_size_[0])
    particles.cells_map[particle_id]->updateX(
        location[0] - particles.domain_size_[0], location[1], location[2]);
  if (location[1] < particles.left_corner_coordinates[1])
    particles.cells_map[particle_id]->updateX(
        location[0], location[1] + particles.domain_size_[1], location[2]);
  if (location[1] >
      particles.left_corner_coordinates[1] + particles.domain_size_[1])
    particles.cells_map[particle_id]->updateX(
        location[0], location[1] - particles.domain_size_[1], location[2]);
  if (particles.domain_size_.size() == 3) {
    if (location[2] < particles.left_corner_coordinates[2])
      particles.cells_map[particle_id]->updateX(
          location[0], location[1], location[2] + particles.domain_size_[2]);
    if (location[2] >
        particles.left_corner_coordinates[2] + particles.domain_size_[2])
      particles.cells_map[particle_id]->updateX(
          location[0], location[1], location[2] - particles.domain_size_[2]);
  }

  particles.cells_map[particle_id]->updateOldX(location[0], location[1],
                                               location[2]);
  particles.update_particle_location(particle_id, location);
  if (!particles.is_within_domain(particles.cells_map[particle_id]->getX())) {
    std::cout << "Out of domain" << std::endl;
    std::cout << particles.cells_map[particle_id]->toString() << std::endl;
    std::cout << particles.is_within_domain(
                     particles.cells_map[particle_id]->getOldX())
              << std::endl;
    std::cout << "older location: " << old_x[0] << " " << old_x[1] << " "
              << old_x[2] << std::endl;
    std::cout << particles.is_within_domain(old_x) << std::endl;
  }
}

void BoundaryConditions::handle_outflow_conditions(
    int particle_id, int cell_index, LinkedCellContainer &particles) {

  auto &particle = particles.cells_map[particle_id];
  particle->left_domain = true;
  particles.particles_left_domain++;
}