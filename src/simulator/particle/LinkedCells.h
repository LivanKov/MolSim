#include "ParticleContainer.h"
#include <array>

#pragma once

// 2 or 3 dimensions ???
template <size_t N>
struct Cell{
    std::vector<ParticlePointer> particles;
};


template <size_t N>
class LinkedCells {};


template<>
class LinkedCells<2> {
    public:
        LinkedCells(std::array<double,2>& domain_size);
    
};

template<>
class LinkedCells<3> {
    public:
        LinkedCells(std::array<double,3>& domain_size);
};