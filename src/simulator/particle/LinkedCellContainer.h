#include "ParticleContainer.h"
#include <array>

#pragma once

template <size_t N>
class LinkedCellContainer : public ParticleContainer {

    LinkedCellContainer(std::array<double,N>& domain_size, double r_cutoff);

    public:
        std::array<double,N> domain_size_;
        double r_cutoff_;
        
};