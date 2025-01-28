//
// Created by sebastianpse on 1/27/25.
//

#pragma once
#include "../simulator/particle/container/LinkedCellContainer.h"
#include <iostream>
#include <fstream>




class ParticleProfiler {
    public:
    ParticleProfiler(LinkedCellContainer &particles,
        size_t number_of_bins,
        double x_min,
        double x_max,
        const std::string& output_file);

    void apply_profiler();


    private:
    LinkedCellContainer &particles_;
    size_t number_of_bins_;
    double x_min_;
    double x_max_;
    double bin_width_;
    std::string output_file_;


};

