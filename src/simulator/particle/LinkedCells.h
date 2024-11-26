#include "ParticleContainer.h"
#include <array>

#pragma once

// 2 or 3 dimensions ???
struct Cell{
    std::vector<ParticlePointer> particles;
};


template <size_t N>
class LinkedCells {};


template<>
class LinkedCells<2> {
    public:
        LinkedCells(std::array<double,2>& domain_size, double r_cutoff);
    private: 
        std::vector<std::vector<Cell>> cells_;
        double r_cutoff_;
    
};

template<>
class LinkedCells<3> {
    public:
        LinkedCells(std::array<double,3>& domain_size, double r_cutoff);
    private:
        std::vector<std::vector<std::vector<Cell>>> cells_;
        double r_cutoff_; 
};