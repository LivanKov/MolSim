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
#include "io/input/cli/SimParams.h"

class ThermostatTest : public testing::Test {
protected:

    std::unique_ptr<Thermostat> thermostat;
    LinkedCellContainer particles{};

    void SetUp() override {
        SimParams::enable_thermo = true;
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
   std::array<double, 3>{0, 0, 0}, particles);

    Thermostat unit_thermostat(particles, 100, 110, 3, 0.5, true, true);

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100, 5);
}

TEST_F(ThermostatTest, BoltzmannSecondTest) {
    ParticleGenerator::insertCuboid(
   std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{9, 13, 9}, 1.0, 10.0,
   std::array<double, 3>{0, 0, 0}, particles);

    Thermostat unit_thermostat(particles, 500, 600, 3, 0.5, true, true);

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 500, 15);
}

// checks for a smaller value
TEST_F(ThermostatTest, BoltzmannThirdTest) {
    ParticleGenerator::insertCuboid(
   std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{5, 5, 5}, 1.0, 5,
   std::array<double, 3>{0, 0, 0},  particles);

    Thermostat unit_thermostat(particles, 5, 100, 3, 0.5, true, true);

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 5, 0.5);
}



// ---------------------------- tests for calculate_kinetic_energy() ---------------------------------------------------


// tests the kinetic energy for a simple cuboid with simple parameters for mass = 1.0 and velocity = 1.0
TEST_F(ThermostatTest, KineticEnergySimpleTest) {
    ParticleGenerator::insertCuboid(
    std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{3, 3, 3}, 1.0, 1.0,
    std::array<double, 3>{1.0, 1.0, 1.0}, particles);
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
    std::array<double, 3>{1.0, 1.0, 1.0}, particles);
    // E_kin must be 1_000 * 3 / 2
    ASSERT_EQ(thermostat->calculate_kinetic_energy(), 1000 * 3 / 2);
}

// here we test with an arbitrary cuboid with arbitrary more complex values
TEST_F(ThermostatTest, KineticEnergyMoreComplexParticleTest) {
    ParticleGenerator::insertCuboid(
    std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{12, 5, 7}, 1.0, 3.7,
    std::array<double, 3>{-4.9, 9.8, 0.2}, particles);
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
        std::array<double, 3>{1, 1, 1}, particles);

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
    std::array<double, 3>{1.0, 1.0, 1.0}, particles);

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
    std::array<double, 3>{-4.9, 9.8, 0.2}, particles);

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
        std::array<double, 3>{100, 100, 100}, particles);

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
        std::array<double, 3>{1.0, 1.0}, particles);

    Thermostat unit_thermostat(particles, 100, 150, 2, 0.5,true,false);

    double expected_kinetic_energy = thermostat->calculate_kinetic_energy();
    // we change the dim param to 2  in the temperature formula
    double expected_temperature = (expected_kinetic_energy * 2) / (2 * particles.size());


    unit_thermostat.calculate_current_temperature();
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

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100, 1e-6);
}

// with random more complex values
TEST_F(ThermostatTest, InitializeMediumTest) {

    Particle particle1{{0,0,0}, {-37, 43.5, 12.4}, 54.5, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 537, 110, 3, 10, true, false);

    // ensures that the particles have the velocity that applies to the initial temperature
    unit_thermostat.initialize();

    unit_thermostat.calculate_current_temperature();

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

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 129, 1e-6);
}

// tests for multiple particles with random values
TEST_F(ThermostatTest, InitializeManyParticlesTest) {
    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 7, 13}, 1.0, 53.9,
         std::array<double, 3>{13.3, 412.41, -4123.2}, particles);

    Thermostat unit_thermostat(particles, 1231, 10000, 3, 0.000002, true, false);

    // ensures that the particles have the velocity that applies to the initial temperature
    unit_thermostat.initialize();

    unit_thermostat.calculate_current_temperature();

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

    unit_thermostat.calculate_current_temperature();

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

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100.5, 1e-8);
}

// we check here if the gradual application works with initialization with brown
TEST_F(ThermostatTest, ApplyCheckGradualBrownTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150,3, 0.5, true, true);

    unit_thermostat.calculate_current_temperature();
    // this is the current temp before application
    double current_temp = unit_thermostat.get_current_temperature();

    unit_thermostat.apply();

    unit_thermostat.calculate_current_temperature();

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

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100.000002, 1e-12);
}

// now we turned gradual off -> after one application we must be directly at temperature of 150
TEST_F(ThermostatTest, ApplyDirectlyTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5,false,false);

    unit_thermostat.apply();

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// the same as above with dimension 2
TEST_F(ThermostatTest, ApplyDirectlyDimensionTwoTest) {
    Particle particle1{{0,0}, {1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150, 2, 0.5, false, false);

    unit_thermostat.apply();

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// now the same as above but with brown initialization
TEST_F(ThermostatTest, ApplyDirectlyBrownTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5, false, true);

    unit_thermostat.apply();

    unit_thermostat.calculate_current_temperature();

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

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// tests with many particles
TEST_F(ThermostatTest, ApplyDirectlyManyParticlesTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 7, 13}, 1.0, 53.9,
         std::array<double, 3>{13.3, 412.41, -4123.2}, particles);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5, false, false);

    unit_thermostat.initialize();

    unit_thermostat.apply();

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// checks for a smaller delta_t
TEST_F(ThermostatTest, ApplySmallerGradualTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 7, 13}, 1.0, 53.9,
         std::array<double, 3>{13.3, 412.41, -4123.2}, particles);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.0005, true, false);

    unit_thermostat.initialize();

    unit_thermostat.apply();

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100 + 0.0005, 1e-11);
}

// checks if delta_t is not set -> delta_t = infinity
// in theory gradual application would act as direct application
TEST_F(ThermostatTest, ApplyDeltaTInfTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 7, 13}, 1.0, 53.9,
         std::array<double, 3>{13.3, 412.41, -4123.2}, particles);

    // here gradual is set to true
    Thermostat unit_thermostat(particles, 100, 150, 3);

    unit_thermostat.initialize();

    unit_thermostat.apply();

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-11);
}


//---------------------- tests for cooling, heating and holding the temperature ----------------------------------------

// for the following tests we use multiple particles small delta T and dimension 3

// we heat the temperature directly to 150
TEST_F(ThermostatTest, HeatingDirectlyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.0001, false, true);

    unit_thermostat.apply();

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// tests if it cools directly to 50
TEST_F(ThermostatTest, CoolingDirectlyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 50, 3, 0.0001, false, true);

    unit_thermostat.apply();

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 50, 1e-8);
}

// tests if temperature holds
TEST_F(ThermostatTest, HoldingDirectlyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 100, 3, 0.0001, false, true);

    unit_thermostat.apply();

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100, 1e-8);
}

// tests if temperature gets heated gradually
TEST_F(ThermostatTest, HeatingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 110, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize();

    // here the thermostat gets 100_000 times applied. with our delta t the temperature must rise from 100 to 110
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply();
    }

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 110, 1e-8);
}

// same as above but with cooling down
TEST_F(ThermostatTest, CoolingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 90, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize();

    // here the thermostat gets 100_000 times applied. with our delta t the temperature must drop from 100 to 90
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply();
    }

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 90, 1e-8);
}

// now with holding the temperature
TEST_F(ThermostatTest, HoldingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 100, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize();

    // here the thermostat gets 100_000 times applied
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply();
    }

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100, 1e-8);
}

// now it’s interesting if the temperature gets hold at the target temperature,
// since target gets reached after 50_000 applications
TEST_F(ThermostatTest, HeatingAndHoldingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 105, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize();

    // here the thermostat gets 100_000 times applied
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply();
    }

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 105, 1e-8);
}

// now same with cooling
TEST_F(ThermostatTest, CoolingAndHoldingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 95, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize();

    // here the thermostat gets 100_000 times applied
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply();
    }

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 95, 1e-8);
}


















































// ---------------------- tests for thermostat modification--------------------------------


// ------------------------- calculate_average_velocity() ---------------------------------

// tests with basic values
TEST_F(ThermostatTest, CalculateAverageVelocitySimpleTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 1, 1}, 1.0, 1.0,
         std::array<double, 3>{1.0, 1.0, 1.0}, particles);

    auto actual_average_velocity = thermostat->calculate_average_velocity();

    for (size_t i = 0; i < 3; ++i) {
        ASSERT_EQ(actual_average_velocity[i], 1.0);
    }
}

// scenario: no particles in the system (edge case)
TEST_F(ThermostatTest, CalculateAverageVelocityNoParticleTest) {

    auto actual_average_velocity = thermostat->calculate_average_velocity();

    for (size_t i = 0; i < 3; ++i) {
        ASSERT_EQ(actual_average_velocity[i], 0);
    }
}

// this tests with different velocities values
TEST_F(ThermostatTest, CalculateAverageVelocityDifferentVelocitiesSimpleTest) {

    Particle particle1{{0.0, 0.0, 0.0}, {1.5, 2.0, 0.5}, 2.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {3.0, 1.0, -2.0}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0, 2.0}, {3.0, 6.0, 0.0}, 1.0, 0};
    particles.insert(particle3);

    auto actual_average_velocity = thermostat->calculate_average_velocity();

    ASSERT_EQ(actual_average_velocity[0], 2.5);
    ASSERT_EQ(actual_average_velocity[1], 3.0);
    ASSERT_EQ(actual_average_velocity[2], -0.5);
}

// tests with particles of dimension two
TEST_F(ThermostatTest, CalculateAverageVelocityTwoDimensionalTest) {

    Particle particle1{{0.0, 0.0}, {1.5, 2.0}, 2.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0}, {3.0, 1.0}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0}, {3.0, 6.0}, 1.0, 0};
    particles.insert(particle3);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5,false,false);

    auto actual_average_velocity = unit_thermostat.calculate_average_velocity();

    ASSERT_EQ(actual_average_velocity[0], 2.5);
    ASSERT_EQ(actual_average_velocity[1], 3.0);
    // is 2D, therefore must be 0
    ASSERT_EQ(actual_average_velocity[2], 0);
}

// edge case: no particles in system
TEST_F(ThermostatTest, CalculateAverageVelocityNoParticlesTest) {

    auto actual_average_velocity = thermostat->calculate_average_velocity();

    for (size_t i = 0; i < 3; ++i) {
        ASSERT_EQ(actual_average_velocity[i], 0.0);
    }
}

// tests for particles with zero velocity
TEST_F(ThermostatTest, CalculateAverageVelocityZeroVelocityTest) {

    Particle particle1{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0, 2.0}, {0.0, 0.0, 0.0}, 1.0, 0};
    particles.insert(particle3);

    auto actual_average_velocity = thermostat->calculate_average_velocity();

    for (size_t i = 0; i < 3; ++i) {
        ASSERT_EQ(actual_average_velocity[i], 0.0);
    }
}

// all velocities are negative
TEST_F(ThermostatTest, CalculateAverageVelocityNegativeVelocitiesTest) {

    Particle particle1{{0.0, 0.0, 0.0}, {-1.0, -1.0, -1.0}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {-2.0, -2.0, -2.0}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0, 2.0}, {-3.0, -3.0, -3.0}, 1.0, 0};
    particles.insert(particle3);

    auto actual_average_velocity = thermostat->calculate_average_velocity();

    ASSERT_EQ(actual_average_velocity[0], -2.0);
    ASSERT_EQ(actual_average_velocity[1], -2.0);
    ASSERT_EQ(actual_average_velocity[2], -2.0);
}

// more complex test with 10 particles with different values
TEST_F(ThermostatTest, CalculateAverageVelocityComplexTest) {
    // Insert 10 particles with different velocities in each dimension
    Particle particle1{{0.0, 0.0, 0.0}, {1.0, 2.0, 3.0}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {2.0, -1.0, 1.0}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0, 2.0}, {-1.0, 3.0, -1.0}, 1.0, 0};
    particles.insert(particle3);

    Particle particle4{{3.0, 3.0, 3.0}, {4.0, -2.0, 2.0}, 1.0, 0};
    particles.insert(particle4);

    Particle particle5{{4.0, 4.0, 4.0}, {1.0, 1.0, 1.0}, 1.0, 0};
    particles.insert(particle5);

    Particle particle6{{5.0, 5.0, 5.0}, {3.0, 0.0, -2.0}, 1.0, 0};
    particles.insert(particle6);

    Particle particle7{{6.0, 6.0, 6.0}, {0.0, 2.0, 4.0}, 1.0, 0};
    particles.insert(particle7);

    Particle particle8{{7.0, 7.0, 7.0}, {-2.0, -1.0, -3.0}, 1.0, 0};
    particles.insert(particle8);

    Particle particle9{{8.0, 8.0, 8.0}, {1.5, 1.5, 1.5}, 1.0, 0};
    particles.insert(particle9);

    Particle particle10{{9.0, 9.0, 9.0}, {0.5, 0.5, -2.5}, 1.0, 0};
    particles.insert(particle10);


    auto actual_average_velocity = thermostat->calculate_average_velocity();

    // (1.0 + 2.0 + - 1.0 + 4.0 + 1.0 + 3.0 + 0 - 2.0 + 1.5 + 0.5) / 10 = 1.0
    ASSERT_EQ(actual_average_velocity[0], 1.0);
    // (2.0 - 1.0 + 3.0 + -2.0 + 1.0 + 0 + 2.0 - 1.0 + 1.5 + 0.5) / 10 = 0.6
    ASSERT_EQ(actual_average_velocity[1], 0.6);
    // (3.0 + 1.0 - 1.0 + 2.0 + 1.0 - 2.0 + 4.0 - 3.0 + 1.5 - 2.5) / 10 = 0.4
    ASSERT_EQ(actual_average_velocity[2], 0.4);
}

// ------------ determine_thermal_motion() -------------------------

// simple thermal motion test with simple values
TEST_F(ThermostatTest, DetermineThermalMotionSimpleTest) {
    Particle particle1{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0, 2.0}, {3.0, 3.0, 3.0}, 1.0, 0};
    particles.insert(particle3);

    thermostat -> determine_thermal_motion();

    ASSERT_EQ(particles[0].getThermalMotion()[0], -1.0);
    ASSERT_EQ(particles[1].getThermalMotion()[0], 0.0);
    ASSERT_EQ(particles[2].getThermalMotion()[0], 1.0);

    ASSERT_EQ(particles[0].getThermalMotion()[1], -1.0);
    ASSERT_EQ(particles[1].getThermalMotion()[1], 0.0);
    ASSERT_EQ(particles[2].getThermalMotion()[1], 1.0);

    ASSERT_EQ(particles[0].getThermalMotion()[2], -1.0);
    ASSERT_EQ(particles[1].getThermalMotion()[2], 0.0);
    ASSERT_EQ(particles[2].getThermalMotion()[2], 1.0);
}

TEST_F(ThermostatTest, DetermineThermalMotionDifferentValuesTest) {
    // Initialize particles with different velocities
    Particle particle1{{0.0, 0.0, 0.0}, {4, 5, 6}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {7, 8, 9}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0, 2.0}, {10, 11, 12}, 1.0, 0};
    particles.insert(particle3);

    thermostat->determine_thermal_motion();


    // average_velocity[0] must be (4 + 7 + 10) / 3 = 21 / 3 = 7
    ASSERT_EQ(particles[0].getThermalMotion()[0], 4 - 7);
    ASSERT_EQ(particles[1].getThermalMotion()[0], 7 - 7);
    ASSERT_EQ(particles[2].getThermalMotion()[0], 10 - 7);

    // average_velocity[1] must be (5 + 8 + 11) / 3 = 24 / 3 = 8
    ASSERT_EQ(particles[0].getThermalMotion()[1], 5 - 8);
    ASSERT_EQ(particles[1].getThermalMotion()[1], 8 - 8);
    ASSERT_EQ(particles[2].getThermalMotion()[1], 11 - 8);

    // average_velocity[2] must be (6 + 9 + 12) / 3 = 27 / 3 = 9
    ASSERT_EQ(particles[0].getThermalMotion()[2], 6 - 9);
    ASSERT_EQ(particles[1].getThermalMotion()[2], 9 - 9);
    ASSERT_EQ(particles[2].getThermalMotion()[2], 12 - 9);
}

// this test is more complex. We use random values for every particle
TEST_F(ThermostatTest, DetermineThermalMotionComplexTest) {
    Particle particle1{{0.0, 0.0, 0.0}, {13.2, -20.6, 10.1}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {-5.4, 7.8, -12.3}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0, 2.0}, {3.3, -8.7, 6.6}, 1.0, 0};
    particles.insert(particle3);

    Particle particle4{{3.0, 3.0, 3.0}, {9.1, 1.2, -4.4}, 1.0, 0};
    particles.insert(particle4);

    Particle particle5{{4.0, 4.0, 4.0}, {-11.0, 0.5, 2.7}, 1.0, 0};
    particles.insert(particle5);

    thermostat->determine_thermal_motion();


    // average_velocity[0] = (13.2 + (-5.4) + 3.3 + 9.1 - 11) / 5 = 9.2 / 5 = 1.84
    ASSERT_NEAR(particles[0].getThermalMotion()[0], 13.2 - 1.84, 1e-6);
    ASSERT_NEAR(particles[1].getThermalMotion()[0], -5.4 - 1.84, 1e-6);
    ASSERT_NEAR(particles[2].getThermalMotion()[0], 3.3 - 1.84, 1e-6);
    ASSERT_NEAR(particles[3].getThermalMotion()[0], 9.1 - 1.84, 1e-6);
    ASSERT_NEAR(particles[4].getThermalMotion()[0], -11.0 - 1.84, 1e-6);

    // average_velocity[1] = (-20.6 + 7.8 - 8.7 + 1.2 + 0.5) / 5 = -19.8 / 5 = -3.96
    ASSERT_NEAR(particles[0].getThermalMotion()[1], -20.6 - (-3.96), 1e-6);
    ASSERT_NEAR(particles[1].getThermalMotion()[1], 7.8 - (-3.96), 1e-6);
    ASSERT_NEAR(particles[2].getThermalMotion()[1], -8.7 - (-3.96), 1e-6);
    ASSERT_NEAR(particles[3].getThermalMotion()[1], 1.2 - (-3.96), 1e-6);
    ASSERT_NEAR(particles[4].getThermalMotion()[1], 0.5 - (-3.96), 1e-6);

    // average_velocity[2] = (10.1 - 12.3 + 6.6 - 4.4 + 2.7) / 5 = 2.7 / 5 = 0.54
    ASSERT_NEAR(particles[0].getThermalMotion()[2], 10.1 - 0.54, 1e-6);
    ASSERT_NEAR(particles[1].getThermalMotion()[2], -12.3 - 0.54, 1e-6);
    ASSERT_NEAR(particles[2].getThermalMotion()[2], 6.6 - 0.54, 1e-6);
    ASSERT_NEAR(particles[3].getThermalMotion()[2], -4.4 - 0.54, 1e-6);
    ASSERT_NEAR(particles[4].getThermalMotion()[2], 2.7 - 0.54, 1e-6);
}

// again but in 2D
TEST_F(ThermostatTest, DetermineThermalMotionComplexTwoDimensionalTest) {
    Particle particle1{{0.0, 0.0}, {13.2, -20.6}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0}, {-5.4, 7.8}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0}, {3.3, -8.7}, 1.0, 0};
    particles.insert(particle3);

    Particle particle4{{3.0, 3.0}, {9.1, 1.2}, 1.0, 0};
    particles.insert(particle4);

    Particle particle5{{4.0, 4.0}, {-11.0, 0.5}, 1.0, 0};
    particles.insert(particle5);

    thermostat->determine_thermal_motion();


    // average_velocity[0] = (13.2 + (-5.4) + 3.3 + 9.1 - 11) / 5 = 9.2 / 5 = 1.84
    ASSERT_NEAR(particles[0].getThermalMotion()[0], 13.2 - 1.84, 1e-6);
    ASSERT_NEAR(particles[1].getThermalMotion()[0], -5.4 - 1.84, 1e-6);
    ASSERT_NEAR(particles[2].getThermalMotion()[0], 3.3 - 1.84, 1e-6);
    ASSERT_NEAR(particles[3].getThermalMotion()[0], 9.1 - 1.84, 1e-6);
    ASSERT_NEAR(particles[4].getThermalMotion()[0], -11.0 - 1.84, 1e-6);

    // average_velocity[1] = (-20.6 + 7.8 - 8.7 + 1.2 + 0.5) / 5 = -19.8 / 5 = -3.96
    ASSERT_NEAR(particles[0].getThermalMotion()[1], -20.6 - (-3.96), 1e-6);
    ASSERT_NEAR(particles[1].getThermalMotion()[1], 7.8 - (-3.96), 1e-6);
    ASSERT_NEAR(particles[2].getThermalMotion()[1], -8.7 - (-3.96), 1e-6);
    ASSERT_NEAR(particles[3].getThermalMotion()[1], 1.2 - (-3.96), 1e-6);
    ASSERT_NEAR(particles[4].getThermalMotion()[1], 0.5 - (-3.96), 1e-6);

    // must be 0 for every z-thermal_motion
    ASSERT_EQ(particles[0].getThermalMotion()[2], 0);
    ASSERT_EQ(particles[1].getThermalMotion()[2], 0);
    ASSERT_EQ(particles[2].getThermalMotion()[2], 0);
    ASSERT_EQ(particles[3].getThermalMotion()[2], 0);
    ASSERT_EQ(particles[4].getThermalMotion()[2], 0);
}

// this tests with many particles (100). Since every particle moves in the same direction,
// the thermal motion must be 0 for every particle
TEST_F(ThermostatTest, DetermineThermalMotionManyParticlesTest) {

    ParticleGenerator::insertCuboid(
            std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 1}, 1.0, 1.0,
            std::array<double, 3>{1.0, 1.0, 1.0}, particles);

    thermostat->determine_thermal_motion();

    for (size_t i = 0; i < particles.size(); ++i) {
        ASSERT_EQ(particles[i].getThermalMotion()[0], 0);
        ASSERT_EQ(particles[i].getThermalMotion()[1], 0);
        ASSERT_EQ(particles[i].getThermalMotion()[2], 0);
    }
}


// ------------ calculate_kinetic_energy_new() -------------------------

// simple test with one moving particle -> must be 0
TEST_F(ThermostatTest, CalculateKineticEnergyNewSimpleTest) {

    Particle particle1{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1.0, 0};
    particles.insert(particle1);

    auto kinetic_energy = thermostat->calculate_kinetic_energy_new();

    ASSERT_EQ(kinetic_energy, 0);
}

// all particles move in same direction -> no kinetic energy
TEST_F(ThermostatTest, CalculateKineticEnergyNewMultipleParticlesSimpleTest) {

    Particle particle1{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{10.0, 10.0, 10.0}, {1.0, 1.0, 1.0}, 1.0, 0};
    particles.insert(particle1);

    Particle particle3{{20.0, 20.0, 20.0}, {1.0, 1.0, 1.0}, 1.0, 0};
    particles.insert(particle1);

    auto kinetic_energy = thermostat->calculate_kinetic_energy_new();

    ASSERT_EQ(kinetic_energy, 0);
}

// simple tests with particles with other velocity values
TEST_F(ThermostatTest, CalculateKineticEnergyNewDifferentValuesSimpleTest) {

    Particle particle1{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0, 2.0}, {3.0, 3.0, 3.0}, 1.0, 0};
    particles.insert(particle3);


    auto kinetic_energy_actual = thermostat->calculate_kinetic_energy_new();

    //   particle 1     particle 2   particle 3
    // (1 * 3*(-1)^2) / 2  +   0   +   (1 * 3*1^2) / 2 = 3
    auto kinetic_energy_expected = 3;

    ASSERT_EQ(kinetic_energy_actual, kinetic_energy_expected);
}

// more complex test, with other velocity values for each particle
TEST_F(ThermostatTest, CalculateKineticEnergyNewComplexTest) {
    // Initialize particles with different velocities
    Particle particle1{{0.0, 0.0, 0.0}, {4, 5, 6}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {7, 8, 9}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0, 2.0}, {10, 11, 12}, 1.0, 0};
    particles.insert(particle3);

    auto kinetic_energy_actual = thermostat->calculate_kinetic_energy_new();

    // (1 * (3^2+3^2)) / 2 + (1 * (3^2+3^2)) / 2 + (1 * (3^2+3^2)) / 2 = 27
    auto kinetic_energy_expected = 27;

    ASSERT_EQ(kinetic_energy_actual, kinetic_energy_expected);
}

// here we calculate for more complex values
TEST_F(ThermostatTest, CalculateKineticEnergyNewComplexOtherMassTest) {
    // Initialize particles with different velocities
    Particle particle1{{0.0, 0.0, 0.0}, {4, 5, 6}, 4.7, 0};
    particles.insert(particle1);

    Particle particle2{{10.0, 10.0, 10.0}, {7, 8, 9}, 2.4, 0};
    particles.insert(particle2);

    Particle particle3{{20.0, 20.0, 20.0}, {10, 11, 12}, 0.6, 0};
    particles.insert(particle3);

    auto kinetic_energy_actual = thermostat->calculate_kinetic_energy_new();


    // particle1_part: 0.5 * 4.7 * (3^2 + 3^2 + 3^2)
    // particle2_part: 0.5 * 2.4 * 0
    // particle3_part: 0.5 * 0.6 * (3^2 + 3^2 + 3^2)
    // particle1_part + particle2_part + particle3_part = kin_energy
    auto kinetic_energy_expected = 63.45 + 8.1;

    ASSERT_EQ(kinetic_energy_actual, kinetic_energy_expected);
}

// here we have many particles. but we must have 0 kin_energy because all particles have the same velocity
TEST_F(ThermostatTest, CalculateKineticEnergyManyParticlesSameVelocityTest) {

    ParticleGenerator::insertCuboid(
            std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 1}, 5.0, 1.0,
            std::array<double, 3>{1.0, 1.0, 1.0}, particles);

    auto kinetic_energy_actual = thermostat->calculate_kinetic_energy_new();

    auto kinetic_energy_expected = 0;

    ASSERT_EQ(kinetic_energy_actual, kinetic_energy_expected);
}


// --------------- calculate_current_temperature_new() --------------------

// simple test with one moving particle -> must be 0
TEST_F(ThermostatTest, CalculateCurrentTemperatureNewSimpleTest) {

    Particle particle1{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1.0, 0};
    particles.insert(particle1);

    thermostat->calculate_current_temperature_new();

    auto temperature_expected = 0;

    ASSERT_EQ(thermostat->get_current_temperature(), temperature_expected);
}

// all particles move in same direction -> no temperature
TEST_F(ThermostatTest, CalculateCurrentTemperatureNewMultipleParticlesSimpleTest) {

    Particle particle1{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{10.0, 10.0, 10.0}, {1.0, 1.0, 1.0}, 1.0, 0};
    particles.insert(particle1);

    Particle particle3{{20.0, 20.0, 20.0}, {1.0, 1.0, 1.0}, 1.0, 0};
    particles.insert(particle1);

    thermostat->calculate_current_temperature_new();

    auto temperature_expected = 0;

    ASSERT_EQ(thermostat->get_current_temperature(), temperature_expected);
}

// simple tests with particles with other velocity values
TEST_F(ThermostatTest, CalculateCurrentTemperatureNewDifferentValuesSimpleTest) {

    Particle particle1{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0, 2.0}, {3.0, 3.0, 3.0}, 1.0, 0};
    particles.insert(particle3);


    thermostat->calculate_current_temperature_new();

    // (3*2)/(3*3)
    auto temperature_expected = 6.0 / 9.0;

    ASSERT_EQ(thermostat->get_current_temperature(), temperature_expected);
}

// more complex test, with other velocity values for each particle
TEST_F(ThermostatTest,  CalculateCurrentTemperatureNewComplexTest) {
    // Initialize particles with different velocities
    Particle particle1{{0.0, 0.0, 0.0}, {4, 5, 6}, 1.0, 0};
    particles.insert(particle1);

    Particle particle2{{1.0, 1.0, 1.0}, {7, 8, 9}, 1.0, 0};
    particles.insert(particle2);

    Particle particle3{{2.0, 2.0, 2.0}, {10, 11, 12}, 1.0, 0};
    particles.insert(particle3);

    thermostat->calculate_current_temperature_new();

    auto temperature_expected = (27.0 * 2.0) / (3.0 * 3.0);

    ASSERT_EQ(thermostat->get_current_temperature(), temperature_expected);
}

// here we test the method with many particles (200) that have different velocities
TEST_F(ThermostatTest, CalculateCurrentTemperatureNewManyParticlesDifferentVelocityTest) {

    ParticleGenerator::insertCuboid(
            std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 1}, 2.0, 1.0,
            std::array<double, 3>{1.0, 1.0, 1.0}, particles);

    ParticleGenerator::insertCuboid(
            std::array<double, 3>{100, 0, 0}, std::array<size_t, 3>{10, 10, 1}, 2.0, 1.0,
            std::array<double, 3>{3.0, 1.0, 1.0}, particles);

    thermostat->calculate_current_temperature_new();

    // kin energy is 100, because each particle has 99 that have same velocity
    // and 100 that have other (x-variable: 3-1 = 2)
    // means: 100 * (2 * 1^2) / 2 = 100
    auto temperature_expected = (100.0 * 2.0) / (3.0 * 200.0);

    ASSERT_EQ(thermostat->get_current_temperature(), temperature_expected);
}

// --------------------- initialize_new() ------------------------

// we test our new initialize method. this method only scales the thermal motion
// according to the initial temperature
// this is just a simple test with two particles with slightly different velocities
TEST_F(ThermostatTest, InitializeNewSimpleTest) {

    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);
    Particle particle2{{0,0,0}, {3.0, 3.0, 3.0}, 1, 0};
    particles.insert(particle2);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5,true,false);

    // we have a current temp of 1 in the moment

    // this method ensures to scale the velocities that lead to a kinetic energy
    // that leads to a current temperature of 100
    unit_thermostat.initialize_new();

    // we calculate the current temperature
    unit_thermostat.calculate_current_temperature_new();

    // now we check the temperature, it should be the initial temp
    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100, 1e-6);
}

// here we check if the velocities get scaled to the supposed ones
TEST_F(ThermostatTest, InitializeNewSimpleVelocityTest) {

    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);
    Particle particle2{{0,0,0}, {3.0, 3.0, 3.0}, 1, 0};
    particles.insert(particle2);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5,true,false);

    // we have a current temp of 1 in the moment

    // this method ensures to scale the velocities that lead to a kinetic energy
    // that leads to a current temperature of 100
    unit_thermostat.initialize_new();

    // we calculate the current temperature
    unit_thermostat.calculate_current_temperature_new();

    // scaling factor must be sqrt(100/1) = 10
    // thermal motion of particle1 is -1 and of particle2 it is 1
    // means velocity for each coordinate must be updated for
    // particle1: 2 - 1*10 = -8
    // particle2: 2 + 1*10 = 12

    ASSERT_EQ(particles[0].getV()[0], -8.0);
    ASSERT_EQ(particles[0].getV()[1], -8.0);
    ASSERT_EQ(particles[0].getV()[2], -8.0);

    ASSERT_EQ(particles[1].getV()[0], 12);
    ASSERT_EQ(particles[1].getV()[1], 12);
    ASSERT_EQ(particles[1].getV()[2], 12);
}

// particles have no thermal motion in the beginning but should have one
// when initializing with brownian motion
TEST_F(ThermostatTest, InitializeNewZeroThermalMotionButBrownianTest) {

    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);
    Particle particle2{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle2);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5, false, true);

    // we calculate the current temperature
    unit_thermostat.calculate_current_temperature_new();

    // now we check the temperature, it should be greater than 0
    ASSERT_GT(unit_thermostat.get_current_temperature(), 0);
}




// ---------------------- apply_new() ---------------------------------


// we check here if the gradual application works, firstly with a large delta_t
TEST_F(ThermostatTest, ApplyNewCheckGradualTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);

    Particle particle2{{0,0,0}, {3.0, 3.0, 3.0}, 1, 0};
    particles.insert(particle2);

    Thermostat unit_thermostat(particles, 100, 150,3, 0.5,true,false);

    // this method is only used to test the apply() method (at least at the moment)
    // this method ensures to scale the velocities that lead to a kinetic energy
    // that leads to a current temperature of 100
    unit_thermostat.initialize_new();

    // now the thermostat gets applied. since we use gradual, the temperature must be 100.5 now. (init_t + delta_t)
    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100.5, 1e-8);
}

// checks for dimension 2
TEST_F(ThermostatTest, ApplyNewCheckGradualDimensionTwoTest) {
    Particle particle1{{0,0}, {1.0, 1.0}, 1, 0};
    particles.insert(particle1);
    Particle particle2{{0,0,0}, {3.0, 3.0, 3.0}, 1, 0};
    particles.insert(particle2);


    Thermostat unit_thermostat(particles, 100, 150,2, 0.5,true,false);

    // this method is only used to test the apply() method (at least at the moment)
    // this method ensures to scale the velocities that lead to a kinetic energy
    // that leads to a current temperature of 100
    unit_thermostat.initialize_new();

    // now the thermostat gets applied. since we use gradual, the temperature must be 100.5 now. (init_t + delta_t)
    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100.5, 1e-8);
}


// we check here if the gradual application works with initialization with brown
TEST_F(ThermostatTest, ApplyNewCheckGradualBrownTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);
    Particle particle2{{0,0,0}, {3.0, 3.0, 3.0}, 1, 0};
    particles.insert(particle2);

    Thermostat unit_thermostat(particles, 100, 150,3, 0.5, true, true);

    unit_thermostat.calculate_current_temperature_new();
    // this is the current temp before application
    double current_temp = unit_thermostat.get_current_temperature();

    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), current_temp + 0.5, 1e-8);
}

// here we apply 2 times with gradual
TEST_F(ThermostatTest, ApplyNewCheckGradualTwoTimesTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);
    Particle particle2{{0,0,0}, {3.0, 3.0, 3.0}, 1, 0};
    particles.insert(particle2);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.000001,true,false);

    unit_thermostat.initialize_new();

    unit_thermostat.apply_new();
    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100.000002, 1e-12);
}

// here we test the new apply method with direct application of thermostat
TEST_F(ThermostatTest, ApplyNewDirectlyTest) {

    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);
    Particle particle2{{0,0,0}, {3.0, 3.0, 3.0}, 1, 0};
    particles.insert(particle2);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5,false,false);


    // this method is only used to test the apply() method (at least at the moment)
    // this method ensures to scale the velocities that lead to a kinetic energy
    // that leads to a current temperature of 100
    unit_thermostat.initialize_new();


    // now the thermostat gets applied.
    unit_thermostat.apply_new();

    // we calculate the new temperature
    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// the same as above with dimension 2
TEST_F(ThermostatTest, ApplyNewDirectlyDimensionTwoTest) {
    Particle particle1{{0,0}, {1.0, 1.0}, 1, 0};
    particles.insert(particle1);
    Particle particle2{{0,0,0}, {3.0, 3.0, 3.0}, 1, 0};
    particles.insert(particle2);

    Thermostat unit_thermostat(particles, 100, 150, 2, 0.5, false, false);

    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// now the same as above but with brown initialization
TEST_F(ThermostatTest, ApplyNewDirectlyBrownTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);
    Particle particle2{{0,0,0}, {3.0, 3.0, 3.0}, 1, 0};
    particles.insert(particle2);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5, false, true);

    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// tests if after applying 2 times thermostat if temperature stays the same
TEST_F(ThermostatTest, ApplyNewDirectlyTwoTimesTest) {
    Particle particle1{{0,0,0}, {1.0, 1.0, 1.0}, 1, 0};
    particles.insert(particle1);
    Particle particle2{{0,0,0}, {3.0, 3.0, 3.0}, 1, 0};
    particles.insert(particle2);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5,false,false);

    unit_thermostat.initialize_new();

    unit_thermostat.apply_new();
    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// tests with many particles with same velocity, but one with other values
TEST_F(ThermostatTest, ApplyNewDirectlyManyParticlesTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 7, 13}, 1.0, 53.9,
         std::array<double, 3>{13.3, 412.41, -4123.2}, particles);
    Particle particle2{{0,0,0}, {3.0, 3.0, 3.0}, 1, 0};
    particles.insert(particle2);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.5, false, false);

    unit_thermostat.initialize_new();

    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-8);
}

// checks for a smaller delta_t
TEST_F(ThermostatTest, ApplyNewSmallerGradualTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 7, 13}, 1.0, 53.9,
         std::array<double, 3>{13.3, 412.41, -4123.2}, particles);
    Particle particle2{{0,0,0}, {3.0, 3.0, 3.0}, 1, 0};
    particles.insert(particle2);

    Thermostat unit_thermostat(particles, 100, 150, 3, 0.0005, true, false);

    unit_thermostat.initialize_new();

    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100 + 0.0005, 1e-8);
}

// checks if delta_t is not set -> delta_t = infinity
// in theory gradual application would act as direct application
TEST_F(ThermostatTest, ApplyNewDeltaTInfTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 7, 13}, 1.0, 53.9,
         std::array<double, 3>{13.3, 412.41, -4123.2}, particles);

    // here gradual is set to true
    Thermostat unit_thermostat(particles, 100, 150, 3);

    unit_thermostat.initialize_new();

    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 150, 1e-11);
}













// ---------------- heating, cooling and holding temperature tests of thermostat modification --------------------------

// we heat the temperature directly with one large cuboid init with brownian motion
TEST_F(ThermostatTest, NewHeatingDirectlyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 40, 60, 3, 0.0001, false, true);

    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 60, 1e-8);
}

// noe the same but with multiple cuboids with different velocities
TEST_F(ThermostatTest, NewHeatingDirectlyMultipleCuboidsTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{1.7, -35.5, 3.3}, particles);
    ParticleGenerator::insertCuboid(
         std::array<double, 3>{100, 100, 100}, std::array<size_t, 3>{1, 10, 25}, 1.0, 5.7,
         std::array<double, 3>{99.7, 0.7, -23.1}, particles);
    ParticleGenerator::insertCuboid(
         std::array<double, 3>{200, 200, 200}, std::array<size_t, 3>{6, 12, 8}, 1.0, 3.1,
         std::array<double, 3>{10.0, -5.5, -5.8}, particles);

    Thermostat unit_thermostat(particles, 40, 60, 3, 0.0001, false, true);

    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 60, 1e-8);
}

// tests if it cools directly to 50
TEST_F(ThermostatTest, NewCoolingDirectlyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 50, 3, 0.0001, false, true);

    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 50, 1e-8);
}

// again, the same with more different particles
TEST_F(ThermostatTest, NewCoolingDirectlyMultipleCuboidsTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{1.7, -35.5, 3.3}, particles);
    ParticleGenerator::insertCuboid(
         std::array<double, 3>{100, 100, 100}, std::array<size_t, 3>{1, 10, 25}, 1.0, 5.7,
         std::array<double, 3>{99.7, 0.7, -23.1}, particles);
    ParticleGenerator::insertCuboid(
         std::array<double, 3>{200, 200, 200}, std::array<size_t, 3>{6, 12, 8}, 1.0, 3.1,
         std::array<double, 3>{10.0, -5.5, -5.8}, particles);

    Thermostat unit_thermostat(particles, 100, 50, 3, 0.0001, false, true);

    unit_thermostat.apply();

    unit_thermostat.calculate_current_temperature();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 50, 1e-8);
}

// tests if temperature holds
TEST_F(ThermostatTest, NewHoldingDirectlyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 100, 3, 0.0001, false, true);

    unit_thermostat.apply_new();

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100, 1e-8);
}

// tests if temperature gets heated gradually
TEST_F(ThermostatTest, NewHeatingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 110, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize_new();

    // here the thermostat gets 100_000 times applied. with our delta t the temperature must rise from 100 to 110
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply_new();
    }

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 110, 1e-8);
}

// same as above but with cooling down
TEST_F(ThermostatTest, NewCoolingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 90, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize_new();

    // here the thermostat gets 100_000 times applied. with our delta t the temperature must drop from 100 to 90
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply_new();
    }

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 90, 1e-8);
}

// now with holding the temperature
TEST_F(ThermostatTest, NewHoldingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 100, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize_new();

    // here the thermostat gets 100_000 times applied
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply_new();
    }

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 100, 1e-8);
}

// now it’s interesting if the temperature gets hold at the target temperature,
// since target gets reached after 50_000 applications
TEST_F(ThermostatTest, NewHeatingAndHoldingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 105, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize_new();

    // here the thermostat gets 100_000 times applied
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply_new();
    }

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 105, 1e-8);
}

// now same with cooling
TEST_F(ThermostatTest, NewCoolingAndHoldingGraduallyTest) {

    ParticleGenerator::insertCuboid(
         std::array<double, 3>{0, 0, 0}, std::array<size_t, 3>{10, 10, 5}, 1.0, 5.0,
         std::array<double, 3>{10.0, -5.5, 3.3}, particles);

    Thermostat unit_thermostat(particles, 100, 95, 3, 0.0001, true, true);

    // with this call we ensure that the real init temp is very near to 100
    // better for testing purposes
    unit_thermostat.initialize_new();

    // here the thermostat gets 100_000 times applied
    for (int i = 0; i < 100000; ++i) {
        unit_thermostat.apply_new();
    }

    unit_thermostat.calculate_current_temperature_new();

    ASSERT_NEAR(unit_thermostat.get_current_temperature(), 95, 1e-8);
}