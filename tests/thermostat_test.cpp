//
// Created by sebastianpse on 12/8/24.
//

# include "simulator/Thermostat.h"
#include <gtest/gtest.h>
#include <iostream>
#include <numeric>
#include <array>
#include <memory>

#include "simulator/particle/ParticleGenerator.h"

class ThermostatTest : public testing::Test {
protected:

    void SetUp() override {
        particles = ParticleContainer();
        thermostat = std::make_unique<Thermostat>(300.0, 100, 350.0, 0.5);
    }

    std::unique_ptr<Thermostat> thermostat;
    ParticleContainer particles;
};

// tests the kinetic energy for a simple cuboid with simple parameters for mass = 1.0 and velocity = 1.0
TEST_F(ThermostatTest, SimpleKineticEnergyTestCase) {
    ParticleGenerator::insertCuboid(
    std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{3, 3, 3}, 1.0, 1.0,
    std::array<double, 3>{1.0, 1.0, 1.0}, 0, particles);

    // since we have 3x3x3 = 27 particles the kinetic energy must be according to the formula 27 * ((1*3)/2) = 40.5
    ASSERT_EQ(thermostat->calculate_kinetic_energy(particles), 40.5);
}

// tests the output for zero particles
TEST_F(ThermostatTest, EdgeKineticEnergyTestCase) {
    ASSERT_EQ(thermostat->calculate_kinetic_energy(particles), 0.0);
}

// here we test with an arbitrary cuboid with arbitrary more complex values
TEST_F(ThermostatTest, StandardKineticEnergyTest) {
    ParticleGenerator::insertCuboid(
    std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{12, 5, 7}, 1.0, 3.7,
    std::array<double, 3>{-4.9, 9.8, 0.2}, 0, particles);
    //                                       our "sum"   mass     dot product of velocity vector
    constexpr auto expected_result = (12 * 5 * 7) * 3.7 * (4.9 * 4.9 + 9.8 * 9.8 + 0.2 * 0.2) / 2;
    ASSERT_NEAR(thermostat->calculate_kinetic_energy(particles), expected_result, 0.00000001);
}



