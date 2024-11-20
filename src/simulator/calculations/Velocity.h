#include "../particle/ParticleContainer.h"


class Velocity{
    static void calculate(ParticleContainer& particles, double time_delta);
    virtual ~Velocity() = 0;
};