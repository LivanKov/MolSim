
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"
#include "utils/GetoptWrapper.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <unordered_map>


/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 * calculate the position for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

constexpr double start_time = 0;
double end_time = 1;
double delta_t = 0.014;

//e : time-end
//d : delta
//i : input path
//o : output  path
//s : sparse output : 1(default) -> enabled, 0 -> disabled

//help: need to implement later


// If the input file is not found -> exception
// If the output file is not found: Create it or throw exception?

constexpr std::string_view optstring {"e:d:i:o:s:"};

// this variable should be set via input, it can be helpful for testing
// controls how often the simulator state is plotted
// sparse = every 10th iteration
// !sparse = every iteration
bool sparse_output = true;

// TODO: what data structure to pick?
std::list<Particle> particles;

int main(int argc, char *argsv[]) {

  /*std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc < 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename" << std::endl;
  }*/

  GetoptWrapper wrapper{optstring};
  std::unordered_map<char,char*> opt_map = wrapper.parse(argc,argsv);

  for(auto const& [a, b] : opt_map){
    std::cout << a << ", " << b << std::endl;
  }
  
  end_time = atof(opt_map['e']);
  delta_t = atof(opt_map['d']);

  FileReader fileReader;
  fileReader.readFile(particles, opt_map['i']);

  double current_time = start_time;

  int iteration = 0;

  for(auto &p : particles) {
    std::cout << p.toString() << std::endl;
  }


  // Calculate position based on velocity: Störmer Verlet
  // Calculate velocity in the next step: Störmer Verlet
  // Calculate force: Simple force calculation

  // Force is initially set to 0
  // Need a way to handle collisions

  // Consider gravitational pull of bigger bodies

  // for this loop, we assume: current x, current f and current v are known
  std::cout << "Starting a simulation with:\n"
            << "\tStart time: " << start_time << "\n"
            << "\tEnd time: " << end_time << "\n"
            << "\tDelta: " << delta_t << "\n";
  while (current_time < end_time) {
    // calculate new x
    calculateX();
    // calculate new f
    calculateF();
    // calculate new v
    calculateV();

    iteration++;
    if(sparse_output && iteration % 10 == 0)
      plotParticles(iteration);
    else
      plotParticles(iteration);
    std::cout << "Iteration " << iteration << " finished." << std::endl;

    current_time += delta_t;
  }

  std::cout << "output written. Terminating..." << std::endl;
  return 0;
}

void calculateF() {
  std::list<Particle>::iterator iterator;
  iterator = particles.begin();

  for (auto &p1 : particles) {
    std::array<double, 3> force_copy = p1.getF();
    for (auto &p2 : particles) {
      double f_x, f_y, f_z = 0;
      double distance = std::sqrt(std::pow(p1.getX()[0] - p2.getX()[0], 2) +
                                  std::pow(p1.getX()[1] - p2.getX()[1], 2) +
                                  std::pow(p1.getX()[2] - p2.getX()[2], 2));
      if (!(p1 == p2)) {
        f_x = (p2.getX()[0] - p1.getX()[0]) * (p1.getM() * p2.getM()) /
              pow(distance, 3);
        f_y = (p2.getX()[1] - p1.getX()[1]) * (p1.getM() * p2.getM()) /
              pow(distance, 3);
        f_z = (p2.getX()[2] - p1.getX()[2]) * (p1.getM() * p2.getM()) /
              pow(distance, 3);
        p1.updateF(p1.getF()[0] + f_x, p1.getF()[1] + f_y, p1.getF()[2] + f_z);
      }
    }
    p1.updateOldF(force_copy[0], force_copy[1], force_copy[2]);
  }
}

void calculateX() {
  for (auto &p : particles) {
    std::array<double, 3> location_copy(p.getX());
    if (p.getF()[0] != 0 || p.getF()[1] != 0 || p.getF()[2] != 0) {
      location_copy[0] = pow(delta_t, 2) * p.getF()[0] / (2 * p.getM());
      location_copy[1] = pow(delta_t, 2) * p.getF()[1] / (2 * p.getM());
      location_copy[2] = pow(delta_t, 2) * p.getF()[2] / (2 * p.getM());
    }
    location_copy[0] = location_copy[0] + delta_t * p.getV()[0];
    location_copy[1] = location_copy[1] + delta_t * p.getV()[1];
    location_copy[2] = location_copy[2] + delta_t * p.getV()[2];
    p.updateX(location_copy[0], location_copy[1], location_copy[2]);
  }
}

void calculateV() {
  for (auto &p : particles) {
    std::array<double, 3> velocity_copy(p.getV());
    velocity_copy[0] = velocity_copy[0] +
                       delta_t * (p.getOldF()[0] + p.getF()[0]) / 2 * p.getM();
    velocity_copy[1] = velocity_copy[1] +
                       delta_t * (p.getOldF()[1] + p.getF()[1]) / 2 * p.getM();
    velocity_copy[2] = velocity_copy[2] +
                       delta_t * (p.getOldF()[2] + p.getF()[2]) / 2 * p.getM();
    p.updateV(velocity_copy[0], velocity_copy[1], velocity_copy[2]);
  }
}

void plotParticles(int iteration) {

  std::string out_name("MD_vtk");
  std::string output_path("../output");
  outputWriter::XYZWriter writer;
  writer.plotParticles(particles, out_name, output_path, iteration);
}
