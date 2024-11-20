#include "../particle/ParticleContainer.h"
#include "Calculation.h"

class Velocity : AbstractPolicy {   
    static void run(ParticleContainer& particles, double time_delta);
    virtual ~Velocity() = 0;
};