#include "LinkedCellContainer.h"
#include <cmath>

size_t Cell::size() const {
    return particles.size();
}

ParticlePointer Cell::operator[](size_t index) {
    return particles[index];
}


LinkedCellContainer::LinkedCellContainer(
    std::initializer_list<double> domain_size, double r_cutoff,
    std::initializer_list<double> left_corner_coordinates)
    : domain_size_{domain_size}, r_cutoff_{r_cutoff},
      left_corner_coordinates{left_corner_coordinates}, x{0}, y{0}, z{0}, ParticleContainer{} {
        if(domain_size.size() != 3 && domain_size.size() != 2) {
            throw std::invalid_argument("Domain size must have 2 or 3 elements");
        }
        x = static_cast<size_t>(domain_size_[0] / r_cutoff);
        y = static_cast<size_t>(domain_size_[1] / r_cutoff);
        z = domain_size.size() == 3 ? static_cast<size_t>((domain_size_[2] / r_cutoff)) : 1;
        cells = std::vector<Cell>(x * y * z, Cell());     
      }

void LinkedCellContainer::insert(Particle &p) {
    ParticlePointer p_ptr = std::make_shared<Particle>(p);
    std::array<double, 3> position = p.getX();
    size_t i = static_cast<size_t>((position[0] - left_corner_coordinates[0]) / r_cutoff_);
    size_t j = static_cast<size_t>((position[1] - left_corner_coordinates[1]) / r_cutoff_);
    size_t k = domain_size_.size() == 3 ? static_cast<size_t>((position[2] - left_corner_coordinates[2]) / r_cutoff_) : 0;
    size_t index = i + j * x + k * x * y;
    cells[index].particles.push_back(p_ptr);
    create_pairs(p_ptr);
    _particle_container.push_back(p_ptr);
}

bool LinkedCellContainer::is_within_domain(const std::array<double,3>& position){
    return position[0] >= left_corner_coordinates[0] && position[0] <= left_corner_coordinates[0] + domain_size_[0] &&
           position[1] >= left_corner_coordinates[1] && position[1] <= left_corner_coordinates[1] + domain_size_[1] &&
           position[2] >= left_corner_coordinates[2] && position[2] <= left_corner_coordinates[2] + domain_size_[2];
}


void LinkedCellContainer::update_particle_location(ParticlePointer p, std::array<double, 3> &old_position) {
    if(is_within_domain(old_position)) {
        size_t i = static_cast<size_t>((old_position[0] - left_corner_coordinates[0]) / r_cutoff_);
        size_t j = static_cast<size_t>((old_position[1] - left_corner_coordinates[1]) / r_cutoff_);
        size_t k = domain_size_.size() == 3 ? static_cast<size_t>((old_position[2] - left_corner_coordinates[2]) / r_cutoff_) : 0;
        size_t old_index = i + j * x + k * x * y;
        cells[old_index].particles.erase(std::remove(cells[old_index].particles.begin(), cells[old_index].particles.end(), p), cells[old_index].particles.end());
    }
    if(is_within_domain(p->getX())) {
        size_t i = static_cast<size_t>((p->getX()[0] - left_corner_coordinates[0]) / r_cutoff_);
        size_t j = static_cast<size_t>((p->getX()[1] - left_corner_coordinates[1]) / r_cutoff_);
        size_t k = domain_size_.size() == 3 ? static_cast<size_t>((p->getX()[2] - left_corner_coordinates[2]) / r_cutoff_) : 0;
        size_t index = i + j * x + k * x * y;
        cells[index].particles.push_back(p);
    }else{
        p->left_domain = true;
    }

}

Cell &LinkedCellContainer::get_cell(size_t index) {
    return cells[index];
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
            neighbours.insert(neighbours.end(), cells[index - 1].particles.begin(), cells[index - 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x].particles.begin(), cells[index - x].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x - 1].particles.begin(), cells[index - x - 1].particles.end());
        } else if(i == x - 1 && j == 0) {
            neighbours.insert(neighbours.end(), cells[index - 1].particles.begin(), cells[index - 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index + x].particles.begin(), cells[index + x].particles.end());
            neighbours.insert(neighbours.end(), cells[index + x - 1].particles.begin(), cells[index + x - 1].particles.end());
        } else if(i == 0 && j == y - 1) {
            neighbours.insert(neighbours.end(), cells[index + 1].particles.begin(), cells[index + 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x].particles.begin(), cells[index - x].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x + 1].particles.begin(), cells[index - x + 1].particles.end());
        } else if(i == 0 && j == 0) {
            neighbours.insert(neighbours.end(), cells[index + 1].particles.begin(), cells[index + 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index + x].particles.begin(), cells[index + x].particles.end());
            neighbours.insert(neighbours.end(), cells[index + x + 1].particles.begin(), cells[index + x + 1].particles.end());
        } else if(i == x - 1) {
            neighbours.insert(neighbours.end(), cells[index - 1].particles.begin(), cells[index - 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index + x].particles.begin(), cells[index + x].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x].particles.begin(), cells[index - x].particles.end());
            neighbours.insert(neighbours.end(), cells[index + x - 1].particles.begin(), cells[index + x - 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x - 1].particles.begin(), cells[index - x - 1].particles.end());
        } else if(j == y - 1) {
            neighbours.insert(neighbours.end(), cells[index + 1].particles.begin(), cells[index + 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x].particles.begin(), cells[index - x].particles.end());
            neighbours.insert(neighbours.end(), cells[index - 1].particles.begin(), cells[index - 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x + 1].particles.begin(), cells[index - x + 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x - 1].particles.begin(), cells[index - x - 1].particles.end());
        } else if(i == 0) {
            neighbours.insert(neighbours.end(), cells[index + 1].particles.begin(), cells[index + 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index + x].particles.begin(), cells[index + x].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x].particles.begin(), cells[index - x].particles.end());
            neighbours.insert(neighbours.end(), cells[index + x + 1].particles.begin(), cells[index + x + 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x + 1].particles.begin(), cells[index - x + 1].particles.end());
        } else if(j == 0) {
            neighbours.insert(neighbours.end(), cells[index + 1].particles.begin(), cells[index + 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x].particles.begin(), cells[index - x].particles.end());
            neighbours.insert(neighbours.end(), cells[index - 1].particles.begin(), cells[index - 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index + x - 1].particles.begin(), cells[index + x - 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x - 1].particles.begin(), cells[index - x -1].particles.end());
        } else {
            neighbours.insert(neighbours.end(), cells[index + 1].particles.begin(), cells[index + 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index - 1].particles.begin(), cells[index - 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index + x].particles.begin(), cells[index + x].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x].particles.begin(), cells[index - x].particles.end());
            neighbours.insert(neighbours.end(), cells[index + x + 1].particles.begin(), cells[index + x + 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index + x - 1].particles.begin(), cells[index + x - 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x + 1].particles.begin(), cells[index - x + 1].particles.end());
            neighbours.insert(neighbours.end(), cells[index - x - 1].particles.begin(), cells[index - x - 1].particles.end());
        }
    } else {

    }

    return neighbours;
}









