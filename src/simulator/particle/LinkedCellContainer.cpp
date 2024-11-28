#include "LinkedCellContainer.h"

template <size_t N>
LinkedCellContainer<N>::LinkedCellContainer(
    std::array<double, N> &domain_size, double r_cutoff,
    std::array<double, 3> &left_corner_coordinates)
    : domain_size_{domain_size}, r_cutoff_{r_cutoff},
      left_corner_coordinates{left_corner_coordinates},
      cells_{left_corner_coordinates, domain_size, r_cutoff} {}