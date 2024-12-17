#include "BoundaryConditions.h"
#include "../particle/container/LinkedCellContainer.h"


void BoundaryConditions::run(LinkedCellContainer &particles) {

    for(auto& index : particles.halo_cell_indices){
        for(auto& id : particles.cells[index].particle_ids){
               
        }
    }
    
}