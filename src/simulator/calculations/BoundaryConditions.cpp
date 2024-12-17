#include "BoundaryConditions.h"
#include "../particle/container/LinkedCellContainer.h"


void BoundaryConditions::run(LinkedCellContainer &particles) {

    for(auto& cell_index : particles.halo_cell_indices){
        for(auto& particle_id : particles.cells[cell_index].particle_ids){
            auto position = particles.cells[cell_index].placement;
            if(particles.placement_map[position] == BoundaryCondition::Reflecting){
                handle_reflect_conditions(particle_id, cell_index, particles);
            }
        }
    }

    for(auto& particle_id : particles.particles_outbound){
        auto& particle = particles.cells_map[particle_id];
        if(!particle->left_domain){
            auto cell_index = particles.get_cell_index(particle->getOldX());
            handle_outflow_conditions(particle_id, cell_index, particles);
        }
    }
    
}


void BoundaryConditions::handle_reflect_conditions(int particle_id, int cell_index, LinkedCellContainer &particles) {
  auto& velocity = particles.cells_map[particle_id]->getV();
  auto& cell = particles.cells[cell_index];
  // Left boundary or Right boundary
  if (cell.placement == Placement::LEFT || cell.placement == Placement::RIGHT)
    particles.cells_map[particle_id]->updateV(-velocity[0], velocity[1], velocity[2]);

  // Bottom boundary
  if (cell.placement == Placement::BOTTOM || cell.placement == Placement::TOP)      
    particles.cells_map[particle_id]->updateV(velocity[0], -velocity[1], velocity[2]);

  // Front or Back
  if (particles.z > 1 && (cell.placement == Placement::FRONT || cell.placement == Placement::BACK))
    particles.cells_map[particle_id]->updateV(velocity[0], velocity[1], -velocity[2]);

  //handle corners
  if(cell.placement == Placement::BOTTOM_LEFT_CORNER || cell.placement == Placement::TOP_RIGHT_CORNER || cell.placement == Placement::TOP_LEFT_CORNER || cell.placement == Placement::BOTTOM_RIGHT_CORNER)
    particles.cells_map[particle_id]->updateV(-velocity[0], -velocity[1], velocity[2]);
}

void BoundaryConditions::handle_periodic_conditions(LinkedCellContainer &particles) {
  
}

void BoundaryConditions::handle_outflow_conditions(int particle_id, int cell_index, LinkedCellContainer &particles) {
    auto& cell = particles.cells[cell_index];
    if(particles.placement_map[cell.placement] == BoundaryCondition::Outflow){
        auto& particle = particles.cells_map[particle_id];
        particle->left_domain = true;
        particles.particles_left_domain++;
    }
}