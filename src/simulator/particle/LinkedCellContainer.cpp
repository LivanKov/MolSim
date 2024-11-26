#include "LinkedCellContainer.h"

template<size_t N>
LinkedCellContainer<N>::LinkedCellContainer(std::array<double, N>& domain_size, double r_cutoff) : domain_size_{domain_size}, r_cutoff_{r_cutoff} {
    
}