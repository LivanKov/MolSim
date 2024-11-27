#include "LinkedCells.h"
#include <array>

LinkedCells<2>::LinkedCells(std::array<double, 2> &domain_size, double r_cutoff) : r_cutoff_{r_cutoff}, width{domain_size[0]}, height{domain_size[1]}, depth{r_cutoff} {
  cells_ = std::vector<std::vector<Cell>>(width, std::vector<Cell>(height));
}

void LinkedCells<2>::insert_particles(std::vector<ParticlePointer> &particles) {
    for(auto &particle : particles){
        auto coord = particle->getX();
        int x = coord[0]/r_cutoff_;
        int y = coord[1]/r_cutoff_;
        cells_[x][y].particles.push_back(particle);
    }
  // ...
}


LinkedCells<3>::LinkedCells(std::array<double, 3> &domain_size, double r_cutoff) : r_cutoff_{r_cutoff}, width{domain_size[0]}, height{domain_size[1]}, depth{domain_size[2]} {
  cells_ = std::vector<std::vector<std::vector<Cell>>>(width, std::vector<std::vector<Cell>>(height, std::vector<Cell>(depth)));
}

void LinkedCells<3>::insert_particles(std::vector<ParticlePointer> &particles) {
  // ...
}
