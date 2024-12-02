#include <gtest/gtest.h>
#include "simulator/particle/LinkedCellContainer.h"
#include "simulator/particle/ParticleGenerator.h"
#include "utils/logger/Logger.h"
#include <array>


class LinkedCellTest : public ::testing::Test {
protected:
    LinkedCellTest() : container{{9.0,9.0}, 3.0, {0.0, 0.0, 0.0}}, container_3d{{9.0, 9.0, 9.0}, 3.0, {0.0, 0.0, 0.0}} {}
    Logger &logger = Logger::getInstance("debug");
    LinkedCellContainer container;
    LinkedCellContainer container_3d;
};

TEST_F(LinkedCellTest, LocationTest) {
    for(size_t i = 0; i < container.cells.size(); ++i){
        EXPECT_EQ(container.cells[i].size(), 0);
    }
    
    EXPECT_TRUE(container.cells.size() == 9);    

    // single particles insertion test

    Particle p1(std::array<double, 3>{5.0, 4.0, 0.0}, std::array<double, 3>{0.0, 0.0, 0.0}, 1.0, 0);

    container.insert(p1);

    EXPECT_EQ(container.cells[4].size(), 1);

    for(size_t i = 0; i < container.cells.size(); ++i){
        if(i != 4){
            EXPECT_EQ(container.cells[i].size(), 0);
        }
    }

    //update particle location, make sure it is removed from old cell and inserted into new cell

    std::array<double,3>old_position = p1.getX();

    container.cells[4][0]->updateX(7.0, 7.0, 0.0);

    container.update_particle_location(container.cells[4][0], old_position);


    //EXPECT_EQ(container.cells[8].size(), 1);

    for(size_t i = 0; i < container.cells.size(); ++i){
        if(i != 8){
            EXPECT_EQ(container.cells[i].size(), 0);
        }
    }

    // leave the domain

    
    std::array<double, 3>another_old_position = container.cells[8][0]->getX();
    container.cells[8][0]->updateX(10.0, 10.0, 0.0);
    container.update_particle_location(container.cells[8][0], another_old_position);

    //Ensure that the particle is still within the container

    EXPECT_TRUE(container.size() == 1);
    
    // Ensure that the corresponding flag has been set

    EXPECT_TRUE(container[0].left_domain);

    // Ensure that the particle is no longer in the cells
    
    for(size_t i = 0; i < container.cells.size(); ++i){
        EXPECT_EQ(container.cells[i].size(), 0);
    }
}


TEST_F(LinkedCellTest, CuboidTest){

    // Insert a 3x3x1 cuboid of particles with a side length of 3.0 and a mass of 1.0
    ParticleGenerator::insertCuboid(std::array<double, 3>{1.5, 1.5, 0.0}, std::array<size_t, 3>{3, 3, 1}, 3.0, 1.0, std::array<double,3>{0.0, 0.0, 0.0}, 0.0, container);

    EXPECT_TRUE(container.size() == 9);

    for(size_t i = 0; i < container.cells.size(); ++i){
        EXPECT_EQ(container.cells[i].size(), 1);
    }

    container.clear();


    // Insert a 3x3x3 cuboid of particles with a side length of 3.0 and a mass of 1.0  

    ParticleGenerator::insertCuboid(std::array<double, 3>{1.5, 1.5, 0.0}, std::array<size_t, 3>{3, 3, 3}, 3.0, 1.0, std::array<double,3>{0.0, 0.0, 0.0}, 0.0, container);

    EXPECT_TRUE(container.size() == 27);

    for(size_t i = 0; i < container.cells.size(); ++i){
        EXPECT_EQ(container.cells[i].size(), 1);
    }
}

TEST_F(LinkedCellTest, NeighbourTest){

    ParticleGenerator::insertCuboid(std::array<double, 3>{1.5, 1.5, 0.0}, std::array<size_t, 3>{3, 3, 1}, 3.0, 1.0, std::array<double,3>{0.0, 0.0, 0.0}, 0.0, container);

    EXPECT_TRUE(container.size() == 9);
    
    Particle p_1 = container[0];

    EXPECT_TRUE(p_1.getX()[0] == 1.5 && p_1.getX()[1] == 1.5 && p_1.getX()[2] == 0.0);

    EXPECT_TRUE(container.get_neighbours(p_1).size() == 3);

    Particle p_2 = container[4];

    EXPECT_TRUE(p_2.getX()[0] == 4.5 && p_2.getX()[1] == 4.5 && p_2.getX()[2] == 0.0);

    EXPECT_TRUE(container.get_neighbours(p_2).size() == 8);

    Particle p_3 = container[7];

    EXPECT_TRUE(p_3.getX()[0] == 4.5 && p_3.getX()[1] == 7.5 && p_3.getX()[2] == 0.0);

    EXPECT_TRUE(container.get_neighbours(p_3).size() == 5);


    //verify every single neigbour for posterity's sake

    auto neighbours = container.get_neighbours(p_3);

    // check all surrounding coordinates
    auto it = std::find_if(neighbours.begin(), neighbours.end(), [](ParticlePointer p) {
        return p->getX()[0] == 7.5 && p->getX()[1] == 7.5 && p->getX()[2] == 0.0;
    });

    EXPECT_TRUE(it != neighbours.end());
    it = std::find_if(neighbours.begin(), neighbours.end(), [](ParticlePointer p) {
        return p->getX()[0] == 7.5 && p->getX()[1] == 4.5 && p->getX()[2] == 0.0;
    });
    EXPECT_TRUE(it != neighbours.end());

    it = std::find_if(neighbours.begin(), neighbours.end(), [](ParticlePointer p) {
        return p->getX()[0] == 4.5 && p->getX()[1] == 4.5 && p->getX()[2] == 0.0;
    });
    EXPECT_TRUE(it != neighbours.end());

    it = std::find_if(neighbours.begin(), neighbours.end(), [](ParticlePointer p) {
        return p->getX()[0] == 1.5 && p->getX()[1] == 4.5 && p->getX()[2] == 0.0;
    });
    EXPECT_TRUE(it != neighbours.end());

    it = std::find_if(neighbours.begin(), neighbours.end(), [&](ParticlePointer p) {
        return p->getX()[0] == 1.5 && p->getX()[1] == 7.5 && p->getX()[2] == 0.0;
    });
    EXPECT_TRUE(it != neighbours.end());
}

    // Insert a disc of particles with a radius of 1.5 and a mass of 1.0