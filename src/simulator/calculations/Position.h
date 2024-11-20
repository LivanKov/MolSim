#include "../particle/ParticleContainer.h"

class Position {
    static void run(ParticleContainer& particles, double time_delta);
    virtual ~Position() = 0;    
};