#include "FileWriter.h"

output::FileWriter::FileWriter(LinkedCellContainer &particles)
    : particles(particles){};

void output::FileWriter::plot_particles(const std::string &filepath,
                                        int iteration){};
