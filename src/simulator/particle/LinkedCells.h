#pragma once

template <size_t N>
class LinkedCells {
public:
  LinkedCells(std::array<double, N> &domain_size, double r_cutoff);
};