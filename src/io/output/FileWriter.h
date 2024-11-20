#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <string>
#include <simulator/particle/ParticleContainer.h>
#include <memory>

namespace output{

class FileWriter {
public:
    FileWriter(ParticleContainer& particles);

    virtual ~FileWriter() = default;
    
    virtual void plot_particles(const std::string& filepath, int iteration) = 0;

    ParticleContainer& particles;
};
}
#endif