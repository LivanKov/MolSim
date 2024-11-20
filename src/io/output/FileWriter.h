#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <string>
#include <particleSim/ParticleContainer.h>

namespace output{

class FileWriter {
public:
    FileWriter(ParticleContainer& particles);

    FileWriter() = default;

    virtual ~FileWriter() = default;

    virtual void write_file(const std::string& filepath, int iteration) = 0;
    
    virtual void plot_particles() = 0;

    ParticleContainer particles;
};
}
#endif