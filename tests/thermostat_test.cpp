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
#include "utils/logger/Logger.h"

class ThermostatTest : public testing::Test {
protected:

    std::unique_ptr<Thermostat> thermostat;
    ParticleContainer particles;

    void SetUp() override {
        particles = ParticleContainer();
        thermostat = std::make_unique<Thermostat>(particles,300, 350, 0.5, false, 3);
    }


};


TEST_F(ThermostatTest, ConstructorTest) {

    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 110, 0.5, true, 3);


    // this is the temp of the particle with its initial values
    double current_temp = 1.5 * 2 / 3 * 1;

    double scaling_factor = std::sqrt(100 / current_temp);

    std::array<double, 3> new_v{1 * scaling_factor, 1 * scaling_factor, 1 * scaling_factor};

    ASSERT_EQ(particles[0].getV(), new_v);
}

// ---------------------------- tests for calculate_kinetic_energy() ---------------------------------------------------


// tests the kinetic energy for a simple cuboid with simple parameters for mass = 1.0 and velocity = 1.0
TEST_F(ThermostatTest, KineticEnergySimpleTest) {
    ParticleGenerator::insertCuboid(
    std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{3, 3, 3}, 1.0, 1.0,
    std::array<double, 3>{1.0, 1.0, 1.0}, 0, particles);
    // since we have 3x3x3 = 27 particles the kinetic energy must be according to the formula 27 * ((1*3)/2) = 40.5
    ASSERT_EQ(thermostat->calculate_kinetic_energy(), 40.5);
}

// tests the output for zero particles
TEST_F(ThermostatTest, KineticEnergyZeroParticlesTest) {
    ASSERT_EQ(thermostat->calculate_kinetic_energy(), 0.0);
}

// tests if the correct E_kin gets calculated with one single particle
TEST_F(ThermostatTest, KineticEnergySingleParticleTest) {
    // insert particle with
    Particle particle{{0,0,0},{2.0, 0, 0},5.0};
    particles.insert(particle);
    // E_kin = 1 * 5 * (2^2 +0+0) / 2 = 10
    ASSERT_EQ(thermostat->calculate_kinetic_energy(), 10.0);
}

// checks if E_kin is 0 when mass 0
TEST_F(ThermostatTest, KineticEnergyZeroMassTest) {
    // insert particle with
    Particle particle{{10,10,10},{2.0, 3.0, -1.0},0.0};
    particles.insert(particle);
    // must be 0
    ASSERT_EQ(thermostat->calculate_kinetic_energy(), 0.0);
}

TEST_F(ThermostatTest, KineticEnergyZeroVelocityTest) {
    // insert particle with
    Particle particle{{-10,5,30},{0.0, 0.0, 0.0},2.7};
    particles.insert(particle);
    // must be 0
    ASSERT_EQ(thermostat->calculate_kinetic_energy(), 0.0);
}

TEST_F(ThermostatTest, KineticEnergyTwoDimensionsTest) {
    Particle particle1{{0.0, 0.0}, {1.5, 2.0}, 2.0};
    particles.insert(particle1);

    // we take our formula: E_kin =
    double e_kin = 1 * 0.5 * 2 * ((1.5 * 1.5) + (2 * 2));

    ASSERT_EQ(thermostat->calculate_kinetic_energy(), e_kin);
}

TEST_F(ThermostatTest, KineticEnergyMixedParticlesTest) {
    Particle particle1{{0.0, 0.0, 0.0}, {1.5, 2.0, 0.5}, 2.0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {3.0, 1.0, -2.0}, 1.0};
    particles.insert(particle2);

    // Particle 1: E_kin = 2 * (1.5^2 + 2^2 + 0.5^2) / 2 = 6.5
    // Particle 2: E_kin = 1 * (3^2 + 1^2 + 2^2) / 2 = 7
    constexpr double expected_energy = 6.5 + 7;
    ASSERT_NEAR(thermostat->calculate_kinetic_energy(), expected_energy, 1e-6);
}


// tests for a large amount of particles 10*10*10 = 1_000
TEST_F(ThermostatTest, KineticEnergyLargeNumberOfParticlesTest) {
    ParticleGenerator::insertCuboid(
    std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 10}, 1.0, 1.0,
    std::array<double, 3>{1.0, 1.0, 1.0}, 0, particles);
    thermostat->calculate_kinetic_energy();
    // E_kin must be 1_000 * 3 / 2
    ASSERT_EQ(thermostat->calculate_kinetic_energy(), 1000 * 3 / 2);
}

// here we test with an arbitrary cuboid with arbitrary more complex values
TEST_F(ThermostatTest, KineticEnergyMoreComplexParticleTest) {
    ParticleGenerator::insertCuboid(
    std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{12, 5, 7}, 1.0, 3.7,
    std::array<double, 3>{-4.9, 9.8, 0.2}, 0, particles);
    //                                       our "sum"   mass     dot product of velocity vector
    constexpr auto expected_result = (12 * 5 * 7) * 3.7 * (4.9 * 4.9 + 9.8 * 9.8 + 0.2 * 0.2) / 2;
    ASSERT_NEAR(thermostat->calculate_kinetic_energy(), expected_result, 1e-6);
}

TEST_F(ThermostatTest, KineticEnergyFiveParticlesTest) {
    // we insert 5 particles with all of them having different values
    Particle particle1{{0.5, 0.5, 0.5}, {2.0, -1.5, 0.5}, 1.5};
    particles.insert(particle1);

    Particle particle2{{-1.0, 0.0, 1.0}, {0.5, 3.0, 1.0}, 2.0};
    particles.insert(particle2);

    Particle particle3{{1.0, 2.0, 3.0}, {1.0, 1.0, -2.0}, 2.5};
    particles.insert(particle3);

    Particle particle4{{-1.5, -1.5, -1.5}, {3.5, 2.0, -1.0}, 1.0};
    particles.insert(particle4);

    Particle particle5{{2.5, -2.5, 0.0}, {-2.0, 0.0, 2.0}, 1.8};
    particles.insert(particle5);

    double expected_kinetic_energy = 0;

    // we use our E_kin formula for every particle
    expected_kinetic_energy += 0.5 * 1.5 * ((2 * 2) + (1.5 * 1.5) + (0.5 * 0.5));
    expected_kinetic_energy += 0.5 * 2 * ((0.5 * 0.5) + (3 * 3) + (1 * 1));
    expected_kinetic_energy += 0.5 * 2.5 * ((1 * 1) + (1 * 1) + (2 * 2));
    expected_kinetic_energy += 0.5 * 1 * ((3.5 * 3.5) + (2 * 2) + (1 * 1));
    expected_kinetic_energy += 0.5 * 1.8 * ((2 * 2) + 0 + (2 * 2));

    // Assert that the computed kinetic energy matches the expected value
    ASSERT_NEAR(thermostat->calculate_kinetic_energy(), expected_kinetic_energy, 1e-8);
}





// ---------------------------- calculate_current_temperature() --------------------------------------------------------

// tests calculate_current_temperature() for zero particles -> edge case
TEST_F(ThermostatTest, CurrentTemperatureNoParticlesTest) {
 //   thermostat->calculate_current_temperature(3);

    ASSERT_EQ(thermostat->get_current_temperature(), 100);  // must be 300 with using 300 as initial temperature
    // the warning method gets also displayed from the logger
}

// tests the method for one particle
TEST_F(ThermostatTest, CurrentTemperatureOneParticleTest) {
    ParticleGenerator::insertCuboid(
        std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{1, 1, 1}, 1.0, 1.0,
        std::array<double, 3>{1, 1, 1}, 0, particles);

    double e_kin = thermostat->calculate_kinetic_energy();

    // apply our formula again
    double expected_temperature = (e_kin * 2) / (3*1);

    thermostat->calculate_current_temperature();
    ASSERT_EQ(thermostat->get_current_temperature(), expected_temperature);
}

// tests the case if dimensions is set to 0
TEST_F(ThermostatTest, CurrentTemperatureZeroDimensionsTest) {
    ParticleGenerator::insertCuboid(
        std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{1, 1, 1}, 1.0, 1.0,
        std::array<double, 3>{1, 1, 1}, 0, particles);

  //  thermostat->calculate_current_temperature(0);
    // must be the initial temperature value
    ASSERT_EQ(thermostat->get_current_temperature(), 100);
}

// tests the case if dimensions is set to 4 -> invalid value for dimensions
TEST_F(ThermostatTest, CurrentTemperatureFourDimensionsTest) {
    ParticleGenerator::insertCuboid(
        std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{1, 1, 1}, 1.0, 1.0,
        std::array<double, 3>{1, 1, 1}, 0, particles);

   // thermostat->calculate_current_temperature(4);
    // must be the initial temperature value
    ASSERT_EQ(thermostat->get_current_temperature(), 100);
}

// case: particles have velocity of 0
TEST_F(ThermostatTest, CurrentTemperatureZeroVelocityParticlesTest) {
    // All particles have zero velocity, no kinetic energy should be present
    ParticleGenerator::insertCuboid(
        std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{3, 3, 3}, 1.0, 0.0,
        std::array<double, 3>{0.0, 0.0, 0.0}, 0, particles);

    thermostat->calculate_current_temperature();

    // Since the kinetic energy is 0, the temperature should also be 0
    ASSERT_EQ(thermostat->get_current_temperature(), 0.0);
}


// simple test for calculate_current_temperature()
TEST_F(ThermostatTest, CurrentTemperatureSimpleTest) {
    ParticleGenerator::insertCuboid(
    std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{3, 3, 3}, 1.0, 1.0,
    std::array<double, 3>{1.0, 1.0, 1.0}, 0, particles);

    // we know from SimpleKineticEnergyTestCase that E_kin = 40.5
    // we can now apply our formula T = (E_kin * 2) / (dimensions * number of particles)
    //
    constexpr auto expected_result = 81 / 79;
    thermostat->calculate_current_temperature();
    ASSERT_EQ(thermostat->get_current_temperature(), expected_result);
}

// in this test we use a more complex cuboid. we use the one from KineticEnergyMoreComplexParticleTest
 TEST_F(ThermostatTest, CurrentTemperatureMediumTest) {
    ParticleGenerator::insertCuboid(
    std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{12, 5, 7}, 1.0, 3.7,
    std::array<double, 3>{-4.9, 9.8, 0.2}, 0, particles);

    auto expected_result_e_kin = thermostat->calculate_kinetic_energy();

    // we can now apply our formula T = (E_kin * 2) / (dimensions * number of particles)
    auto expected_result = (expected_result_e_kin * 2) / (3 * 12*5*7);

    thermostat->calculate_current_temperature();
    ASSERT_EQ(thermostat->get_current_temperature(), expected_result);
}

// use case for high velocity values
TEST_F(ThermostatTest, CurrentTemperatureHighVelocityTest) {
    ParticleGenerator::insertCuboid(
        std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{5, 5, 5}, 1.0, 1.0,
        std::array<double, 3>{100, 100, 100}, 0, particles);

    // Calculate expected kinetic energy with high velocities
    double expected_kinetic_energy = thermostat->calculate_kinetic_energy();

    double expected_temperature = (expected_kinetic_energy * 2) / (3 * particles.size());

    thermostat->calculate_current_temperature();
    ASSERT_NEAR(thermostat->get_current_temperature(), expected_temperature, 1e-6);
}

// this tests for dimensions = 2
TEST_F(ThermostatTest, TwoDimensionalCurrentTemperatureTest) {


    ParticleGenerator::insertCuboid(
        std::array<double, 3>{0, 0}, std::array<size_t, 3>{3, 3, 3}, 1.0, 1.0,
        std::array<double, 3>{1.0, 1.0}, 0, particles);

    double expected_kinetic_energy = thermostat->calculate_kinetic_energy();
    // we change the dim param to 2  in the temperature formula
    double expected_temperature = (expected_kinetic_energy * 2) / (2 * particles.size());

 //   thermostat->calculate_current_temperature(2);
    ASSERT_NEAR(thermostat->get_current_temperature(), expected_temperature, 1e-6);
}

//-------------------------------------- tests for apply() -------------------------------------------------------------


// check this
TEST_F(ThermostatTest, ApplyZeroVelocityTest) {
    ParticleGenerator::insertCuboid(
        std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{3, 3, 3}, 1.0, 1.0,
        std::array<double, 3>{0, 0, 0}, 0, particles);

    Thermostat unit_thermostat(particles, 100, 150, 0.5,false, 3);

    unit_thermostat.apply(true);

    ASSERT_EQ(unit_thermostat.get_current_temperature(), 0);
}

TEST_F(ThermostatTest, ApplyCheckGradualTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150, 0.5,true,3);

    unit_thermostat.apply(true);

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100.5, 1e-8);
}

TEST_F(ThermostatTest, ApplyCheckGradualTwoTimesTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150, 0.5,true,3);

    unit_thermostat.apply(true);
    unit_thermostat.apply(true);


    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 101, 1e-8);
}

TEST_F(ThermostatTest, ApplyDirectlyTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150, 0.5,false,3);

    unit_thermostat.apply(false);


    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

TEST_F(ThermostatTest, ApplyDirectlyTwoTimesTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150, 0.5,false,3);

    unit_thermostat.apply(false);
    unit_thermostat.apply(false);


    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}