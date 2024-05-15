# Parallel Global VLSI Router

This repository contains the implementation of a Parallel Global Router for Very Large Scale Integration (VLSI) design. The focus is on optimizing the routing of connections between pins on a chip during the physical design stage. This involves converting an abstract netlist and logic level design into a chip layout that can be fabricated.

## Background

In VLSI design, routing metal interconnects between component pins is crucial. These interconnects are made of straight metal wires and vertical connections called vias, arranged in multiple layers. Wires travel in one direction: horizontally on odd-numbered layers and vertically on even-numbered layers. Given the NP-Hard nature of this routing problem, we utilize heuristics and OpenMP for parallelization to simplify and expedite the process.

## Authors

- **Deepak S**
  - Email: deepaksridhar13@gmail.com

- **Anurag Manoj**
  - Email: armj2357@gmail.com

## Repository Structure

The repository includes C, C++, and Python files necessary for creating nets, routing, evaluating, and visualizing routes.

## Getting Started

### Step 1: Clone the Repository

Clone the repository to your local machine to access all necessary files.

```bash
git clone https://github.com/mm22b011-deepaks/Parallel-Global-VLSI-Router.git
```

### Step 2: Compile and run the code

You will need OpenMP installed in your computer to compile this code and run.

To compile the Random Test Case Generator, we will use the following command to compile and run

```bash
gcc -fopenmp rand_test_generator.c -o rand_test_generator.o -lm && ./rand_test_generator.o <Number of threads>
```

Replace the number of threads with the number of threads to run the generator code, if no argument is passed, it defaults to 8 threads

To compile the VLSI Router, use the below commands

```bash
g++ -g -fopenmp main.cpp netlist.cpp batch.cpp Grid_Graph.cpp -o VLSI_Router.o
```

VLSI_Router.o takes the following arguments as inputs

```bash
Correct format is ./VLSI_Router.o <input filename> <output filename> <Bounding Box Dimensions> <Number of Iterations for Maze Route> <Number of Threads for Parallelising>
```

The input filename is the output of rand_test_generator.o

Bouding box dimensions are usually in the order of 5,7 or 9

### Step 3: Visualising the routes

Using the visualise.py script given in the repository, we can visualise the routes. To make routes more legible, there is a certain offset that has been given to overlapping routes, you can have the value of the offset in the visualise.py file at your wish

```bash
python3 visualise.py <output filename>
```

### Step 4: Evaluate the routes

Run the input and the output files through the Evaluator.py script 

```bash
python3 Evaluator.py <input filename> <output filename>
```

## Using Makefile
There is a Makefile present in the repository, that can also be used by running the following command

```bash
make
```
