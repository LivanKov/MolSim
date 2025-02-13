#include "BoundaryConditions.h"
#include "../particle/container/LinkedCellContainer.h"
#include <iostream>

void BoundaryConditions::run(LinkedCellContainer &particles) {
  particles.clear_ghost_particles();
  for (auto &cell_index : particles.halo_cell_indices) {
    for (auto &particle_id : particles.cells[cell_index].particle_ids) {
      auto position = particles.cells[cell_index].placement;
      if (particles.placement_map[position] == BoundaryCondition::Periodic) {
        particles.create_ghost_particles(particle_id, cell_index);
      }
    }

    for (size_t i = 0; i < particles.particles_outbound.size(); ++i) {
      auto &particle_id = particles.particles_outbound[i];
      // for (auto &particle_id : particles.particles_outbound) {
      auto &particle = particles.at(particle_id);
      if (!particle->left_domain && !particle->outbound) {
        auto cell_index = particles.get_cell_index(particle->getOldX());
        auto position = particles.cells[cell_index].placement;
        particle->outbound = true;
        if (particles.placement_map[position] == BoundaryCondition::Outflow) {
          handle_outflow_conditions(particle_id, cell_index, particles);
        } else if (particles.placement_map[position] ==
                   BoundaryCondition::Periodic) {
          handle_periodic_conditions(particle_id, cell_index, particles);
        } else if (particles.placement_map[position] ==
                   BoundaryCondition::Reflecting) {
          handle_reflect_conditions(particle_id, cell_index, particles);
        }
      }
    }
  }
}

void BoundaryConditions::handle_reflect_conditions(
    int particle_id, int cell_index, LinkedCellContainer &particles) {
  auto &velocity = particles.at(particle_id)->getV();
  particles.at(particle_id)->outbound = false;
  particles.particles_outbound.erase(
      std::remove(particles.particles_outbound.begin(),
                  particles.particles_outbound.end(), particle_id),
      particles.particles_outbound.end());

  std::array<double, 3> location = particles.at(particle_id)->getX();

  auto newX = location[0];
  auto newY = location[1];
  auto newZ = location[2];
  auto newV_x = velocity[0];
  auto newV_y = velocity[1];
  auto newV_z = velocity[2];

  if (location[0] < particles.left_corner_coordinates[0]) {
    newX = 2 * particles.left_corner_coordinates[0] - location[0];
    newV_x = -velocity[0];
  }

  if (location[0] >
      particles.left_corner_coordinates[0] + particles.domain_size_[0]) {
    newX =
        2 * (particles.left_corner_coordinates[0] + particles.domain_size_[0]) -
        location[0];
    newV_x = -velocity[0];
  }

  if (location[1] < particles.left_corner_coordinates[1]) {
    newY = 2 * particles.left_corner_coordinates[1] - location[1];
    newV_y = -velocity[1];
  }

  if (location[1] >
      particles.left_corner_coordinates[1] + particles.domain_size_[1]) {
    newY =
        2 * (particles.left_corner_coordinates[1] + particles.domain_size_[1]) -
        location[1];
    newV_y = -velocity[1];
  }

  if (particles.domain_size_.size() == 3) {
    if (location[2] < particles.left_corner_coordinates[2]) {
      newZ = 2 * particles.left_corner_coordinates[2] - location[2];
      newV_z = -velocity[2];
    }

    if (location[2] >
        particles.left_corner_coordinates[2] + particles.domain_size_[2]) {
      newZ = 2 * (particles.left_corner_coordinates[2] +
                  particles.domain_size_[2]) -
             location[2];
      newV_z = -velocity[2];
    }
  }

  particles.at(particle_id)->updateX(newX, newY, newZ);
  particles.at(particle_id)->updateV(newV_x, newV_y, newV_z);
  particles.at(particle_id)->updateOldX(location[0], location[1], location[2]);
  particles.update_particle_location(particle_id, location);
}

void BoundaryConditions::handle_periodic_conditions(
    int particle_id, int cell_index, LinkedCellContainer &particles) {
  particles.at(particle_id)->outbound = false;
  particles.particles_outbound.erase(
      std::remove(particles.particles_outbound.begin(),
                  particles.particles_outbound.end(), particle_id),
      particles.particles_outbound.end());

  auto old_x = particles.at(particle_id)->getOldX();

  std::array<double, 3> location = particles.at(particle_id)->getX();

  auto newX = location[0];
  auto newY = location[1];
  auto newZ = location[2];

  if (location[0] < particles.left_corner_coordinates[0]) {
    newX = location[0] + particles.domain_size_[0];
  }

  if (location[0] >
      particles.left_corner_coordinates[0] + particles.domain_size_[0]) {
    newX = location[0] - particles.domain_size_[0];
  }

  if (location[1] < particles.left_corner_coordinates[1]) {
    newY = location[1] + particles.domain_size_[1];
  }

  if (location[1] >
      particles.left_corner_coordinates[1] + particles.domain_size_[1]) {
    newY = location[1] - particles.domain_size_[1];
  }

  if (particles.domain_size_.size() == 3) {
    if (location[2] < particles.left_corner_coordinates[2]) {
      newZ = location[2] + particles.domain_size_[2];
    }

    if (location[2] >
        particles.left_corner_coordinates[2] + particles.domain_size_[2]) {
      newZ = location[2] - particles.domain_size_[2];
    }
  }
  particles.at(particle_id)->updateX(newX, newY, newZ);
  particles.at(particle_id)->updateOldX(location[0], location[1], location[2]);
  particles.update_particle_location(particle_id, location);
  if (!particles.is_within_domain(particles.at(particle_id)->getX())) {
    std::cout << "Out of domain" << std::endl;
    std::cout << particles.at(particle_id)->toString() << std::endl;
    std::cout << particles.is_within_domain(
                     particles.at(particle_id)->getOldX())
              << std::endl;
    std::cout << "older location: " << old_x[0] << " " << old_x[1] << " "
              << old_x[2] << std::endl;
    std::cout << particles.is_within_domain(old_x) << std::endl;
  }
}

void BoundaryConditions::handle_outflow_conditions(
    int particle_id, int cell_index, LinkedCellContainer &particles) {

  auto &particle = particles.at(particle_id);
  particle->left_domain = true;
  particles.particles_left_domain++;
}