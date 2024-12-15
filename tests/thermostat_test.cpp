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
    LinkedCellContainer particles{};

    void SetUp() override {
        thermostat = std::make_unique<Thermostat>(particles,300, 350, 3, 0.5, false, false);
    }


};

// ------------------------- constructor and brownian initialization tests ---------------------------------------------

TEST_F(ThermostatTest, CheckValidConstructorTest) {
    Thermostat unit_thermostat(particles, 300, 350, 3, 0.5, false, true);

    EXPECT_EQ(unit_thermostat.get_current_temperature(), 0); // Initially 0 until particles' velocity is calculated
    EXPECT_EQ(unit_thermostat.get_dimensions(), 3);
    EXPECT_EQ(unit_thermostat.get_gradual(), false);
    EXPECT_EQ(unit_thermostat.get_target_temperature(), 350.0);
}

// tests if invalid dim argument gets recognized
TEST_F(ThermostatTest, InvalidDimensionFourThrowsExceptionTest) {
    EXPECT_THROW(
        Thermostat unit_thermostat(particles, 300, 350, 4), // Invalid dimension: 4
        std::invalid_argument
    );
}

// tests if invalid dim argument gets recognized
TEST_F(ThermostatTest, InvalidDimensionZeroThrowsExceptionTest) {
    EXPECT_THROW(
        Thermostat unit_thermostat(particles, 300, 350, 0), // Invalid dimension: 0
        std::invalid_argument
    );
}



// this tests checks if the current temperature gets set near to the initial temperature
// means it checks if the initialization with brownian motion is correct
TEST_F(ThermostatTest, BoltzmannFirstTest) {
    ParticleGenerator::insertCuboid(
   std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{9, 13, 9}, 1.0, 10.0,
   std::array<double, 3>{0, 0, 0}, 0, particles);

    Thermostat unit_thermostat(particles, 100, 110, 3, 0.5, true, true);

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100, 5);
}

TEST_F(ThermostatTest, BoltzmannSecondTest) {
    ParticleGenerator::insertCuboid(
   std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{9, 13, 9}, 1.0, 10.0,
   std::array<double, 3>{0, 0, 0}, 0, particles);

    Thermostat unit_thermostat(particles, 500, 600, 3, 0.5, true, true);

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 500, 15);
}

// checks for a smaller value
TEST_F(ThermostatTest, BoltzmannThirdTest) {
    ParticleGenerator::insertCuboid(
   std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{5, 5, 5}, 1.0, 5,
   std::array<double, 3>{0, 0, 0}, 0, particles);

    Thermostat unit_thermostat(particles, 5, 100, 3, 0.5, true, true);

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 5, 0.5);
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
    Particle particle{{0,0,0},{2.0, 0, 0},5.0, 0};
    particles.insert(particle);
    // E_kin = 1 * 5 * (2^2 +0+0) / 2 = 10
    ASSERT_EQ(thermostat->calculate_kinetic_energy(), 10.0);
}

// checks if E_kin is 0 when mass 0
TEST_F(ThermostatTest, KineticEnergyZeroMassTest) {
    // insert particle with
    Particle particle{{10,10,10},{2.0, 3.0, -1.0},0.0, 0};
    particles.insert(particle);
    // must be 0
    ASSERT_EQ(thermostat->calculate_kinetic_energy(), 0.0);
}

// tests for particle with 0 velocity -> E_kin should be 0
TEST_F(ThermostatTest, KineticEnergyZeroVelocityTest) {
    // insert particle with
    Particle particle{{-10,5,30},{0.0, 0.0, 0.0},2.7, 0};
    particles.insert(particle);
    // must be 0
    ASSERT_EQ(thermostat->calculate_kinetic_energy(), 0.0);
}

// tests for dimension of 2
TEST_F(ThermostatTest, KineticEnergyTwoDimensionsTest) {
    Particle particle1{{0.0, 0.0}, {1.5, 2.0}, 2.0, 0};
    particles.insert(particle1);

    // we take our formula: E_kin =
    double e_kin = 1 * 0.5 * 2 * ((1.5 * 1.5) + (2 * 2));

    ASSERT_EQ(thermostat->calculate_kinetic_energy(), e_kin);
}

// tests for two different particles with different values
TEST_F(ThermostatTest, KineticEnergyMixedParticlesTest) {
    Particle particle1{{0.0, 0.0, 0.0}, {1.5, 2.0, 0.5}, 2.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {3.0, 1.0, -2.0}, 1.0, 0};
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

// now we test here with 5 different particles
TEST_F(ThermostatTest, KineticEnergyFiveParticlesTest) {
    // we insert 5 particles with all of them having different values
    Particle particle1{{0.5, 0.5, 0.5}, {2.0, -1.5, 0.5}, 1.5, 0};
    particles.insert(particle1);

    Particle particle2{{-1.0, 0.0, 1.0}, {0.5, 3.0, 1.0}, 2.0, 1};
    particles.insert(particle2);

    Particle particle3{{1.0, 2.0, 3.0}, {1.0, 1.0, -2.0}, 2.5, 2};
    particles.insert(particle3);

    Particle particle4{{-1.5, -1.5, -1.5}, {3.5, 2.0, -1.0}, 1.0, 3};
    particles.insert(particle4);

    Particle particle5{{2.5, -2.5, 0.0}, {-2.0, 0.0, 2.0}, 1.8, 4};
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
    thermostat->calculate_current_temperature();

    ASSERT_EQ(thermostat->get_current_temperature(), 0);  // must be 0 since 0 particles in system
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

    Thermostat unit_thermostat(particles, 100, 150, 2, 0.5,true,false);

    double expected_kinetic_energy = thermostat->calculate_kinetic_energy();
    // we change the dim param to 2  in the temperature formula
    double expected_temperature = (expected_kinetic_energy * 2) / (2 * particles.size());

 //   thermostat->calculate_current_temperature(2);
    ASSERT_NEAR(unit_thermostat.get_current_temperature(), expected_temperature, 1e-6);
}


// ------------------------------- tests for initialize() --------------------------------------------------------------

// simple test to check if particle's velocity gets scaled according to init_temp
TEST_F(ThermostatTest, InitializeSimpleTest) {

    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 110, 3, 0.5, true, false);

    // ensures that the particles have the velocity that applies to the initial temperature
    unit_thermostat.initialize();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100, 1e-6);
}

// with random more complex values
TEST_F(ThermostatTest, InitializeMediumTest) {

    Particle particle1{{0,0,0}, {-37, 43.5, 12.4}, 54.5, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 537, 110, 3, 10, true, false);

    // ensures that the particles have the velocity that applies to the initial temperature
    unit_thermostat.initialize();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 537, 1e-6);
}

// with different velocities of two particles
TEST_F(ThermostatTest, InitializeDifferentVelocitiesTest) {

    Particle particle1{{0,0,0}, {-37, 43.5, 12.4}, 54.5, 0};
    particles.insert(particle1);

    Particle particle2{{0,0,0}, {-43, -20, 64.4}, 70.3, 1};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 129, 110, 3, 10, true, false);

    // ensures that the particles have the velocity that applies to the initial temperature
    unit_thermostat.initialize();


    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 129, 1e-6);
}

// tests for multiple particles with random values
TEST_F(ThermostatTest, InitializeManyParticlesTest) {
    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 7, 13}, 1.0, 53.9,
         std::array<double, 3>{13.3, 412.41, -4123.2}, 0, particles);

    Thermostat unit_thermostat(particles, 1231, 10000, 3, 0.000002, true, false);

    // ensures that the particles have the velocity that applies to the initial temperature
    unit_thermostat.initialize();


    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 1231, 1e-6);
}



//-------------------------------------- tests for apply() -------------------------------------------------------------


// we check here if the gradual application works, firstly with a large delta_t
TEST_F(ThermostatTest, ApplyCheckGradualTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150,3, 0.5,true,false);

    // this method is only used to test the apply() method (at least at the moment)
    // this method ensures to scale the velocities that lead to a kinetic energy
    // that leads to a current temperature of 100
    unit_thermostat.initialize();

    // now the thermostat gets applied. since we use gradual, the temperature must be 100.5 now. (init_t + delta_t)
    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100.5, 1e-8);
}

// checks for dimension 2
TEST_F(ThermostatTest, ApplyCheckGradualDimensionTwoTest) {
    Particle particle1{{0,0}, {1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150,2, 0.5,true,false);

    // this method is only used to test the apply() method (at least at the moment)
    // this method ensures to scale the velocities that lead to a kinetic energy
    // that leads to a current temperature of 100
    unit_thermostat.initialize();

    // now the thermostat gets applied. since we use gradual, the temperature must be 100.5 now. (init_t + delta_t)
    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100.5, 1e-8);
}

// we check here if the gradual application works with initialization with brown
TEST_F(ThermostatTest, ApplyCheckGradualBrownTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150,3, 0.5, true, true);

    // this is the current temp before application
    double current_temp = unit_thermostat.get_current_temperature();

    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), current_temp + 0.5, 1e-8);
}

// here we apply 2 times with gradual
TEST_F(ThermostatTest, ApplyCheckGradualTwoTimesTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.000001,true,false);

    unit_thermostat.initialize();

    unit_thermostat.apply();
    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100.000002, 1e-12);
}

// now we turned gradual off -> after one application we must be directly at temperature of 150
TEST_F(ThermostatTest, ApplyDirectlyTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5,false,false);

    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// the same as above with dimension 2
TEST_F(ThermostatTest, ApplyDirectlyDimensionTwoTest) {
    Particle particle1{{0,0}, {1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150, 2, 0.5, false, false);

    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// now the same as above but with brown initialization
TEST_F(ThermostatTest, ApplyDirectlyBrownTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5, false, true);

    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// tests if after applying 2 times thermostat if temperature stays the same
TEST_F(ThermostatTest, ApplyDirectlyTwoTimesTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5,false,false);

    unit_thermostat.initialize();

    unit_thermostat.apply();
    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// tests with many particles
TEST_F(ThermostatTest, ApplyDirectlyManyParticlesTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 7, 13}, 1.0, 53.9,
         std::array<double, 3>{13.3, 412.41, -4123.2}, 0, particles);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5, false, false);

    unit_thermostat.initialize();

    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// checks for a smaller delta_t
TEST_F(ThermostatTest, ApplySmallerGradualTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 7, 13}, 1.0, 53.9,
         std::array<double, 3>{13.3, 412.41, -4123.2}, 0, particles);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.0005, true, false);

    unit_thermostat.initialize();

    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100 + 0.0005, 1e-11);
}

// checks if delta_t is not set -> delta_t = infinity
// in theory gradual application would act as direct application
TEST_F(ThermostatTest, ApplyDeltaTInfTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 7, 13}, 1.0, 53.9,
         std::array<double, 3>{13.3, 412.41, -4123.2}, 0, particles);

    // here gradual is set to true
    Thermostat unit_thermostat(particles, 100, 150, 3);

    unit_thermostat.initialize();

    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-11);
}


//---------------------- tests for cooling, heating and holding the temperature ----------------------------------------

// for the following tests we use multiple particles small delta T and dimension 3

// we heat the temperature directly to 150
TEST_F(ThermostatTest, HeatingDirectlyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, 0, particles);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.0001, false, true);

    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// tests if it cools directly to 50
TEST_F(ThermostatTest, CoolingDirectlyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, 0, particles);

    Thermostat unit_thermostat(particles, 100, 50, 3, 0.0001, false, true);

    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 50, 1e-8);
}

// tests if temperature holds
TEST_F(ThermostatTest, HoldingDirectlyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, 0, particles);

    Thermostat unit_thermostat(particles, 100, 100, 3, 0.0001, false, true);

    unit_thermostat.apply();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100, 1e-8);
}

// tests if temperature gets heated gradually
TEST_F(ThermostatTest, HeatingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, 0, particles);

    Thermostat unit_thermostat(particles, 100, 110, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize();

    // here the thermostat gets 100_000 times applied. with our delta t the temperature must rise from 100 to 110
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply();
    }

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 110, 1e-8);
}

// same as above but with cooling down
TEST_F(ThermostatTest, CoolingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, 0, particles);

    Thermostat unit_thermostat(particles, 100, 90, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize();

    // here the thermostat gets 100_000 times applied. with our delta t the temperature must drop from 100 to 90
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply();
    }

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 90, 1e-8);
}

// now with holding the temperature
TEST_F(ThermostatTest, HoldingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, 0, particles);

    Thermostat unit_thermostat(particles, 100, 100, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize();

    // here the thermostat gets 100_000 times applied
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply();
    }

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100, 1e-8);
}

// now it’s interesting if the temperature gets hold at the target temperature,
// since target gets reached after 50_000 applications
TEST_F(ThermostatTest, HeatingAndHoldingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, 0, particles);

    Thermostat unit_thermostat(particles, 100, 105, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize();

    // here the thermostat gets 100_000 times applied
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply();
    }

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 105, 1e-8);
}

// now same with cooling
TEST_F(ThermostatTest, CoolingAndHoldingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, 0, particles);

    Thermostat unit_thermostat(particles, 100, 95, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize();

    // here the thermostat gets 100_000 times applied
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply();
    }

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 95, 1e-8);
}
