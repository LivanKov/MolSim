MolSim (WS24/25, Group H)
===
This repository contains the code written for the Molecular Dynamics Lab (IN0012) during WS24/25 at TUM.


![CI Build](https://github.com/LivanKov/MolSim/actions/workflows/ci.yml/badge.svg?branch=dev-sheet_2)

![GIF Description](/assets/example1.gif)

![GIF Description](/assets/example2.gif)




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

## XML Input Configuration

The MolSim application now supports using an XML input file for configuring simulation parameters. This XML file allows users to define the simulation setup in a structured way, making it easier to manage complex configuration.

### XML Parameters
The XML input file provides the following parameters:
- **end_time**: The total time for which the simulation should run.
- **time_delta**: The time increment between simulation steps.
- **output_path**: The base name for output files (including filepath).
- **write_frequency**: Specifies how often to write output files (e.g., every N iterations).
- **r_cutoff_radius**: Specifies the distance, from which the force set to 0.
- **domain_size**: size of the inner area including the boundary of the domain.

These parameters are defined in the XML file and are loaded when running the simulation.

### Command Line Arguments
Users can override the parameters defined in the XML file using command line arguments. The command line arguments available are as follows.

| Option            | Description                                                         |
|-------------------|---------------------------------------------------------------------|
| `-h`              | Show this help message and exit                                     |
| `-i <file_path>`  | Specify the input file path                                         |
| `-o <output_name>`| The base name for output files including filepath (override XML)    |
| `-e <end_time>`   | Specify the end time for the simulation to run (override XML)       |
| `-d <time_delta>` | Specify the time increments for each simulation step (override XML) |
| `-t <write_freq>` | write a file for every -t iteration of the run (override XML)       |
| `-x`              | Output files in `.xyz` format instead of the default `.vtu` format  |
| `-l <log_level>`  | Option to choose the logging level                                  |
| `-f`              | Calculate Gravitational Force instead of Lennard-Jones Force        |
| `-n`              | Disable all file outputs                                            | 

### Examples
To run the MolSim program with an XML input file and additional command line arguments:

```sh
./MolSim -i ../input/input_file.xml -e 100.0 -d 0.01 -o ../output/simulation_output -t 5
```

In this example:
- **-i ../input_file.xml** is the XML file that provides the base configuration for the simulation.
- **-e 100.0** overrides the `end_time` specified in the XML file to run the simulation for 100.0 units of time.
- **-d 0.01** overrides the `delta_time` specified in the XML file to set the time increment to 0.01.
- **-o ../output/simulation_output** overrides the `output_path` specified in the XML file to set the output file path.
- **-t 5** overrides the `write_frequence` specified in the XML file to set the write frequency to write a file every 5 iterations.

### XML Example
Here is an example of an XML input file:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<MolSim>
    <simulation_parameters>
        <end_time>5</end_time>
        <delta_time>0.005</delta_time>
        <output_basename>../output/MD_vtk</output_basename>
        <write_frequency>10</write_frequency>
        <r_cutoff_radius>3.0</r_cutoff_radius>
        <gravity>-12.44</gravity>
        <enable_brownian>true</enable_brownian>
        <domain_size>
            <x>180.0</x>
            <y>90.0</y>
            <z>0.0</z>
        </domain_size>
    </simulation_parameters>
    <thermostats>
        <initial_temp>40</initial_temp>
        <target_temp>40</target_temp>
        <delta_temp>0.0005</delta_temp>
        <is_gradual>true</is_gradual>
        <n_thermostats>1000</n_thermostats>
    <thermostats>
    <boundary_conditions>
        <left>Outflow</left>
        <right>Outflow</right>
        <top>Outflow</top>
        <bottom>Outflow</bottom>
        <front>Reflecting</front>
        <back>Reflecting</back>
    </boundary_conditions>
    <discs>
        <disc>
            <center>
                <x>60.0</x>
                <y>25.0</y>
                <z>0.0</z>
            </center>
            <initial_velocity>
                <x>0.0</x>
                <y>-10.0</y>
                <z>0.0</z>
            </initial_velocity>
            <radius>15</radius>
            <mesh_width>1.1225</mesh_width>
            <mass>1.0</mass>
            <epsilon>5.0</epsilon>
            <sigma>1.0</sigma>
        </disc>
    </discs>
    <cuboids>
        <cuboid>
            <coordinate>
                <x>0.0</x>
                <y>0.0</y>
                <z>0.0</z>
            </coordinate>
            <dimensions>
                <x>3</x>
                <y>3</y>
                <z>3</z>
            </dimensions>
            <mesh_width>1.0</mesh_width>
            <mass>1.0</mass>
            <epsilon>5.0</epsilon>
            <sigma>1.0</sigma>
            <initial_velocity>
                <x>0.0</x>
                <y>0.0</y>
                <z>0.0</z>
            </initial_velocity>
        </cuboid>
    </cuboids>
    <particles>
        <particle>
            <position>
                <x>0.0</x>
                <y>0.0</y>
                <z>0.0</z>
            </position>
            <velocity>
                <x>0.0</x>
                <y>0.0</y>
                <z>0.0</z>
            </velocity>
            <mass>1.0</mass>
        </particle>
    </particles>
</MolSim>
```
### Notes
- The XML file is parsed first, and all defined parameters are loaded. If command line arguments are provided, they take precedence over the values in the XML file.
- Using an XML input file makes it easier to manage and reuse complex configurations, especially when working with multiple simulations or running parameter sweeps.

### Scripts Overview

```bash
./scripts/cleanup.sh
```
Deletes all files in the `output` directory. Use this to clear generated files or results from previous runs.

```bash
./scripts/format.sh
```
 Formats all .cpp and .h files in src using clang-format.


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
