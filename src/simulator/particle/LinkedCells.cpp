#include "LinkedCells.h"
#include <array>

LinkedCells<2>::LinkedCells(std::array<double, 2> &domain_size, double r_cutoff) : r_cutoff_{r_cutoff}, width{domain_size[0]}, height{domain_size[1]}, depth{r_cutoff} {
  cells_ = std::vector<std::vector<Cell>>(domain_size[0]/r_cutoff, std::vector<Cell>(domain_size[1]/r_cutoff));
}



LinkedCells<3>::LinkedCells(std::array<double, 3> &domain_size, double r_cutoff) : r_cutoff_{r_cutoff}, width{domain_size[0]}, height{domain_size[1]}, depth{domain_size[2]} {
cells_ = std::vector<std::vector<std::vector<Cell>>>(domain_size[0]/r_cutoff, 
                 std::vector<std::vector<Cell>>(domain_size[1]/r_cutoff,
                 std::vector<Cell>(domain_size[2]/r_cutoff)));
  // ...
}