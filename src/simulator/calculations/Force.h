#include "../particle/ParticleContainer.h"

#define EPSILON 5.0
#define SIGMA 1.0

class Force{
    static void run_verlet(ParticleContainer &particles);
    static void run_lennard_jones(ParticleContainer &particles);
    virtual ~Force() = 0;    
};