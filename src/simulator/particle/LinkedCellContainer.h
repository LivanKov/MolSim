#include "LinkedCells.h"
#include "ParticleContainer.h"
#include <array>

#pragma once

/**
 * @class LinkedCellContainer
 * @brief Class that provides a container for particles that uses linked cells
 * to speed up the computation.Inherits from ParticleContainer.
 * @tparam N Dimension of the container.
 */
template <size_t N> class LinkedCellContainer : public ParticleContainer {
  /**
   * @brief Constructor.
   * @param domain_size Domain size of the container.
   * @param r_cutoff Cutoff radius.
   */
  LinkedCellContainer(std::array<double, N> &domain_size, double r_cutoff,
                      std::array<double, 3> &left_corner_coordinates);

  std::array<double, N> domain_size_;
  std::array<double, 3> left_corner_coordinates;
  double r_cutoff_;
  LinkedCells<N> cells_;
};