#include "LinkedCellContainer.h"
#include <cmath>

LinkedCellContainer::LinkedCellContainer(
    const std::vector<double> &domain_size, double r_cutoff,
    std::array<double, 3> &left_corner_coordinates)
    : domain_size_{domain_size}, r_cutoff_{r_cutoff},
      left_corner_coordinates{left_corner_coordinates} {
        if(domain_size.size() != 3 && domain_size.size() != 2) {
            throw std::invalid_argument("Domain size must have 2 or 3 elements");
        }
        x = static_cast<size_t>(std::ceil(domain_size[0] / r_cutoff));
        y = static_cast<size_t>(std::ceil(domain_size[1] / r_cutoff));
        z = domain_size.size() == 3 ? static_cast<size_t>(std::ceil(domain_size[2] / r_cutoff)) : 1;
        unwrapped_cells_.resize(x * y * z);      
      }