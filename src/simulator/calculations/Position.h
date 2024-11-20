#include "../particle/ParticleContainer.h"
#include "Calculation.h"

class Position : AbstractPolicy {
    static void run(ParticleContainer& particles, double time_delta);
    virtual ~Position() = 0;    
};