#include "LinkedCellContainer.h"
#include <cmath>

LinkedCellContainer::LinkedCellContainer(
    std::initializer_list<double> domain_size, double r_cutoff,
    std::initializer_list<double> left_corner_coordinates)
    : domain_size_{domain_size}, r_cutoff_{r_cutoff},
      left_corner_coordinates{left_corner_coordinates} {
        if(domain_size.size() != 3 && domain_size.size() != 2) {
            throw std::invalid_argument("Domain size must have 2 or 3 elements");
        }
        x = static_cast<size_t>(domain_size_[0] / r_cutoff);
        y = static_cast<size_t>(domain_size_[1] / r_cutoff);
        z = domain_size.size() == 3 ? static_cast<size_t>((domain_size_[2] / r_cutoff)) : 1;
        unwrapped_cells_ = std::vector<Cell>(x * y * z, Cell());     
      }

void LinkedCellContainer::insert(Particle &p) {
    std::array<double, 3> position = p.getX();
    size_t i = static_cast<size_t>((position[0] - left_corner_coordinates[0]) / r_cutoff_);
    size_t j = static_cast<size_t>((position[1] - left_corner_coordinates[1]) / r_cutoff_);
    size_t k = domain_size_.size() == 3 ? static_cast<size_t>((position[2] - left_corner_coordinates[2]) / r_cutoff_) : 0;
    size_t index = i + j * x + k * x * y;
    unwrapped_cells_[index].particles.push_back(std::make_shared<Particle>(p));
}


void LinkedCellContainer::update_particle_location(Particle &p, std::array<double, 3> &old_position) {
    size_t i = static_cast<size_t>((old_position[0] - left_corner_coordinates[0]) / r_cutoff_);
    size_t j = static_cast<size_t>((old_position[1] - left_corner_coordinates[1]) / r_cutoff_);
    size_t k = domain_size_.size() == 3 ? static_cast<size_t>((old_position[2] - left_corner_coordinates[2]) / r_cutoff_) : 0;
    size_t old_index = i + j * x + k * x * y;
    unwrapped_cells_[old_index].particles.erase(std::remove_if(unwrapped_cells_[old_index].particles.begin(), unwrapped_cells_[old_index].particles.end(), [&p](const ParticlePointer &particle) { return *particle == p; }), unwrapped_cells_[old_index].particles.end());
    insert(p);
}

Cell &LinkedCellContainer::get_cell(size_t index) {
    return unwrapped_cells_[index];
}

std::vector<ParticlePointer> LinkedCellContainer::get_neighbours(Particle &p) {
    std::array<double, 3> position = p.getX();
    size_t i = static_cast<size_t>((position[0] - left_corner_coordinates[0]) / r_cutoff_);
    size_t j = static_cast<size_t>((position[1] - left_corner_coordinates[1]) / r_cutoff_);
    size_t k = domain_size_.size() == 3 ? static_cast<size_t>((position[2] - left_corner_coordinates[2]) / r_cutoff_) : 0;
    size_t index = i + j * x + k * x * y;
    std::vector<ParticlePointer> neighbours;
    
    // handle 2d case
    if(domain_size_.size() == 2) {
        if(i == x - 1 && j == y -1){
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - 1].particles.begin(), unwrapped_cells_[index - 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x].particles.begin(), unwrapped_cells_[index - x].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x - 1].particles.begin(), unwrapped_cells_[index - x - 1].particles.end());
        } else if(i == x - 1 && j == 0) {
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - 1].particles.begin(), unwrapped_cells_[index - 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + x].particles.begin(), unwrapped_cells_[index + x].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + x - 1].particles.begin(), unwrapped_cells_[index + x - 1].particles.end());
        } else if(i == 0 && j == y - 1) {
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + 1].particles.begin(), unwrapped_cells_[index + 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x].particles.begin(), unwrapped_cells_[index - x].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x + 1].particles.begin(), unwrapped_cells_[index - x + 1].particles.end());
        } else if(i == 0 && j == 0) {
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + 1].particles.begin(), unwrapped_cells_[index + 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + x].particles.begin(), unwrapped_cells_[index + x].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + x + 1].particles.begin(), unwrapped_cells_[index + x + 1].particles.end());
        } else if(i == x - 1) {
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - 1].particles.begin(), unwrapped_cells_[index - 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + x].particles.begin(), unwrapped_cells_[index + x].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x].particles.begin(), unwrapped_cells_[index - x].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + x - 1].particles.begin(), unwrapped_cells_[index + x - 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x - 1].particles.begin(), unwrapped_cells_[index - x - 1].particles.end());
        } else if(j == y - 1) {
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + 1].particles.begin(), unwrapped_cells_[index + 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x].particles.begin(), unwrapped_cells_[index - x].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - 1].particles.begin(), unwrapped_cells_[index - 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x + 1].particles.begin(), unwrapped_cells_[index - x + 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x - 1].particles.begin(), unwrapped_cells_[index - x - 1].particles.end());
        } else if(i == 0) {
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + 1].particles.begin(), unwrapped_cells_[index + 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + x].particles.begin(), unwrapped_cells_[index + x].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x].particles.begin(), unwrapped_cells_[index - x].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + x + 1].particles.begin(), unwrapped_cells_[index + x + 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x + 1].particles.begin(), unwrapped_cells_[index - x + 1].particles.end());
        } else if(j == 0) {
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + 1].particles.begin(), unwrapped_cells_[index + 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x].particles.begin(), unwrapped_cells_[index - x].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - 1].particles.begin(), unwrapped_cells_[index - 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + x - 1].particles.begin(), unwrapped_cells_[index + x - 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x - 1].particles.begin(), unwrapped_cells_[index - x -1].particles.end());
        } else {
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + 1].particles.begin(), unwrapped_cells_[index + 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - 1].particles.begin(), unwrapped_cells_[index - 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + x].particles.begin(), unwrapped_cells_[index + x].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x].particles.begin(), unwrapped_cells_[index - x].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + x + 1].particles.begin(), unwrapped_cells_[index + x + 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index + x - 1].particles.begin(), unwrapped_cells_[index + x - 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x + 1].particles.begin(), unwrapped_cells_[index - x + 1].particles.end());
            neighbours.insert(neighbours.end(), unwrapped_cells_[index - x - 1].particles.begin(), unwrapped_cells_[index - x - 1].particles.end());
        }
    } else {

    }

    return neighbours;
}









