#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <string>
#include <simulator/ParticleContainer.h>
#include <memory>

namespace output{

class FileWriter {
public:
    FileWriter(std::shared_ptr<ParticleContainer>& particles);

    virtual ~FileWriter() = default;
    
    virtual void plot_particles(const std::string& filepath, int iteration) = 0;

    std::shared_ptr<ParticleContainer> particles;
};
}
#endif