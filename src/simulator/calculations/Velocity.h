#include "../particle/ParticleContainer.h"


class Velocity{
    static void run(ParticleContainer& particles, double time_delta);
    virtual ~Velocity() = 0;
};