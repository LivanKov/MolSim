#include "ParticleContainer.h"
#include <array>

#pragma once

/**
 * @struct Cell
 * @brief Struct that represents a cell in the linked cells container.
 */
struct Cell{
    /**
     * @brief Particles contained within the cell
     */
    std::vector<ParticlePointer> particles;
};

/**
 * @class LinkedCells
 * @brief Container class for LinkedCells.
 * @tparam N Dimension of the container, only supports 2 and 3 dimensions.
 */

template <size_t N>
class LinkedCells {};

/**
 * @brief Specialization of the LinkedCells class for 2 dimensions.
 */

template<>
class LinkedCells<2> {
    public:
        /**
         * @brief Constructor.
         * @param domain_size Domain size of the container.
         * @param r_cutoff Cutoff radius.
         */
        LinkedCells(std::array<double,2>& domain_size, double r_cutoff);
        /**
         * @brief Insert particles into the container.
         * @param particles Vector of ParticlePointer objects.
         */
        void insert_particles(std::vector<ParticlePointer>& particles);
    private: 
        std::vector<std::vector<Cell>> cells_;
        double r_cutoff_;
        double width;
        double height;
        double depth;
    
};

/**
 * @brief Specialization of the LinkedCells class for 3 dimensions.
 */
template<>
class LinkedCells<3> {
    public:
        /**
         * @brief Constructor.
         * @param domain_size Domain size of the container.
         * @param r_cutoff Cutoff radius.
         */
        LinkedCells(std::array<double,3>& domain_size, double r_cutoff);
        /**
         * @brief Insert particles into the container.
         * @param particles Vector of ParticlePointer objects.
         */
        void insert_particles(std::vector<ParticlePointer>& particles);
    private:
        std::vector<std::vector<std::vector<Cell>>> cells_;
        double r_cutoff_;
        double width;
        double height;
        double depth;
};