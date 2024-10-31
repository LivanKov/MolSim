MolSim (WS24/25, Group H)
===
This repository contains the code written for the Molecular Dynamics Lab (IN0012) during WS24/25 at TUM.

![CI Build](https://github.com/LivanKov/MolSim/actions/workflows/ci.yml/badge.svg?branch=dev-StatusBadge)


### Build and Run on Linux 

We recommend using the following compilers:

- **GCC** (GNU Compiler Collection) version 11.1 or later
- **Clang** version 14.0.0 or later

Building this project will require the following tools:

- **CMake**: Version 3.10 or later.
  We recommend installing CMake from the [official website](https://cmake.org/download/), or using your package manager.

- **Xerces-C++ XML Parser**:
  We recommend installing and building the most recent version from the [official website](https://xerces.apache.org/xerces-c/), or using your package manager. 
   
Build the project

```
$ mkdir build
$ cd build
$ cmake ..
``` 

Compile and run the executable

```
$ make
$ ./MolSim
``` 

### Options 


| Option            | Description                                                         |
|-------------------|---------------------------------------------------------------------|
| `-h`              | Show this help message and exit                                     |
| `-o <file_path>`  | Specify the output file path                                        |
| `-i <file_path>`  | Specify the input file path                                         |
| `-e <end_time>`   | Specify the end time for the simulation to run                      |
| `-d <time_delta>` | Specify the time increments for each simulation step                |
| `-t`              | Enable testing mode (writes a file for each iteration of the run)   |
| `-x`              | Output files in `.xyz` format instead of the default `.vtu` format  |


### Examples

- `./MolSim -i data/input.txt -o results/output.txt`: Run with specified input and output paths
- `./MolSim -e 100 -d 0.01`: Run for a specific time duration with specific time steps: 
- `./MolSim -t -x`: Run while outputting a .xyz file for each iteration: 


### Branch naming conventions

- `dev-<sheet number>`: For specific work sheets handed out during the course; should be merged into master branch by deadline.
- `dev-<sheet number/task number>`: For specific tasks on a work sheet; should be merged into specific work sheet branch throughout the week.
- `dev-<feature-name>`: For specific features under development.
- `dev-fix-<bug name>`: For fixing bugs that occur during the project.
- `dev-misc`: For miscellaneous updates like documentation, configuration, or small fixes.

### Participants

- [Bangrui Wu](https://github.com/BangruiW)
- [Sebastian John](https://github.com/sebastian-j-john)
- [Ivan Lomakov](https://github.com/LivanKov)

### Supervisors

- [Manish Kumar Mishra](https://github.com/manishmishra6016)
- [Samuel Newcome](https://github.com/SamNewcome)
