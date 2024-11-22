#include "../particle/ParticleContainer.h"
#include "Calculation.h"

#define EPSILON 5.0
#define SIGMA 1.0

enum Type {
    LENNARD_JONES,
    VERLET
};

struct Force : AbstractPolicy{
    static void run(ParticleContainer &particles, Type type);
private:
    static void lennard_jones(ParticleContainer &particles);
    static void verlet(ParticleContainer &particles);
};