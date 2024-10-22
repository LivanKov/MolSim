MolSim (WS24/25, Group H)
===

This repository contains the code written for the Molecular Dynamics Lab (IN0012) during WS24/25 at TUM


### Build and Run on Linux 

We recommend using the following compilers:

- **GCC** (GNU Compiler Collection) version 11.1 or later
- **Clang** version 14.0.0 or later

Building this project will require the following tools:

- **CMake**: Version 3.10 or later.
 We recommend installing CMake from the [official website](https://cmake.org/download/), or using your package manager.

Build the project

```
$ mkdir build
$ cd build
$ cmake ..
``` 

Compile and run the executable

```
$ make
$ ./<PROJECT_NAME>
``` 

### Branch naming conventions

- `dev-misc`: For miscellaneous updates like documentation, configuration, or small fixes.
- `dev-<sheet number, task number>`: For specific tasks.
- `dev-<feature-name>`: For specific features under development.
- `dev-fix-<bug name>`: For fixing bugs that occur during the project.

### Participants

- [BangruiW](https://github.com/BangruiW)
- [sebastian3007](https://github.com/sebastian3007)
- [LivanKov](https://github.com/LivanKov)