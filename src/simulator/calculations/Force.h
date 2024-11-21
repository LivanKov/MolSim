#include "../particle/ParticleContainer.h"
#include "Calculation.h"

#define EPSILON 5.0
#define SIGMA 1.0

struct Force : AbstractPolicy{
    static void run_verlet(ParticleContainer &particles);
    static void run_lennard_jones(ParticleContainer &particles);
};