#include "LinkedCells.h"
#include <array>

LinkedCells<2>::LinkedCells(std::array<double,3>& left_corner_coordinates, std::array<double, 2> &domain_size, double r_cutoff) : r_cutoff_{r_cutoff}, width{domain_size[0]}, height{domain_size[1]}, depth{r_cutoff} {
  cells_ = std::vector<std::vector<Cell>>(width, std::vector<Cell>(height));
  for(size_t i = cells_.size() - 1; i >= 0; --i){
    for(size_t j = 0; j < cells_[i].size(); ++j){
      if(i == 0 || j == 0 || i == cells_.size() - 1 || j == cells_[i].size() - 1){
        cells_[i][j].is_boundary = true;
      }
      cells_[i][j].left_corner_coordinates[0] = left_corner_coordinates[0] * i * r_cutoff;
      cells_[i][j].left_corner_coordinates[1] = left_corner_coordinates[1] * j * r_cutoff;
    }
  }
}

void LinkedCells<2>::insert_particles(std::vector<ParticlePointer> &particles) {
  for(auto &particle : particles){
    Cell &cell = get_corresponding_cell(particle);
    cell.particles.insert(particle);
  }
}

Cell& LinkedCells<2>::get_corresponding_cell(ParticlePointer &particle) {
  auto coordinates = particle->getX();
  auto x = (coordinates[0] - left_corner_coordinates[0]) / r_cutoff_;
  auto y = (coordinates[1] - left_corner_coordinates[1]) / r_cutoff_;
  return cells_[x][y];
}



LinkedCells<3>::LinkedCells(std::array<double,3>& left_corner_coordinates, std::array<double, 3> &domain_size, double r_cutoff) : r_cutoff_{r_cutoff}, width{domain_size[0]}, height{domain_size[1]}, depth{domain_size[2]} {
  cells_ = std::vector<std::vector<std::vector<Cell>>>(width, std::vector<std::vector<Cell>>(height, std::vector<Cell>(depth)));
  for(size_t i = cells_.size() - 1; i >= 0; --i){
    for(size_t j = 0; j < cells_[i].size(); ++j){
      for(size_t k = 0; k < cells_[i][j].size(); ++k){
        if(i == 0 || j == 0 || k == 0 || i == cells_.size() - 1 || j == cells_[i].size() - 1 || k == cells_[i][j].size() - 1){
          cells_[i][j][k].is_boundary = true;
        }
        cells_[i][j][k].left_corner_coordinates[0] = left_corner_coordinates[0] * i * r_cutoff;
        cells_[i][j][k].left_corner_coordinates[1] = left_corner_coordinates[1] * j * r_cutoff;
        cells_[i][j][k].left_corner_coordinates[2] = left_corner_coordinates[2] * k * r_cutoff;
      }
    }
  }
}

Cell& LinkedCells<3>::get_corresponding_cell(ParticlePointer &particle) {
  auto coordinates = particle->getX();
  auto x = (coordinates[0] - left_corner_coordinates[0]) / r_cutoff_;
  auto y = (coordinates[1] - left_corner_coordinates[1]) / r_cutoff_;
  auto z = (coordinates[2] - left_corner_coordinates[2]) / r_cutoff_;
  return cells_[x][y][z];
}


void LinkedCells<3>::insert_particles(std::vector<ParticlePointer> &particles) {
  for(auto &particle : particles){
    Cell &cell = get_corresponding_cell(particle);
    cell.particles.insert(particle);
  }
}
