#include "LinkedCells.h"
#include <array>

LinkedCells<2>::LinkedCells(std::array<double, 2> &domain_size, double r_cutoff) : r_cutoff_{r_cutoff} {
  cells_ = std::vector<std::vector<Cell>>(domain_size[0], std::vector<Cell>(domain_size[1]));
}



LinkedCells<3>::LinkedCells(std::array<double, 3> &domain_size, double r_cutoff) {
cells_ = std::vector<std::vector<std::vector<Cell>>>(domain_size[0], 
                 std::vector<std::vector<Cell>>(domain_size[1],
                 std::vector<Cell>(domain_size[2])));
  // ...
}