#include "../particle/ParticleContainer.h"

class Position {
    static void calculate(ParticleContainer& particles, double time_delta);
    virtual ~Position() = 0;    
};