# TISCC: Trapped-Ion Surface Code Compiler
[![arXiv](https://img.shields.io/badge/arXiv-####.#####-b31b1b.svg)](https://arxiv.org)

Tool to compile a universal set of surface code operations into low-level operations native to a trapped-ion quantum computer and estimate resources for them.

TISCC compiles patch operations using primitives such as transversal operations over data qubits, rounds of error correction over stabilizer plaquettes, lattice surgery operations between pairs of patches, and corner movements. Circuits generated using TISCC have been verified to yield the appropriate logical transformations through simulation and tomography. See our paper (link above) for more details..

TISCC can be used either as an executable or as a library.

## TISCC Executable

### Build

Standard CMake build:
```shell
$ mkdir build
$ cd build
$ cmake ..
```

The `TISCC` executable will be at the top level of the `build` directory. 

### Usage

Example usage:

```shell
# Generate a circuit and estimate resources for a single round of error correction on a d = 5 surface code patch.
./TISCC -x 5 -z 5 -t 1 -o idle -p -r

# Generate a circuit for a merge operation between two horizontally-adjacent d = 5 surface code patches.
./TISCC -x 5 -z 5 -t 5 -o merge -s double-horiz -p

# Visualize the hardware grid allocated to two horizontally-adjacent d = 5 surface code patches
./TISCC -x 5 -z 5 -t 1 -s double-horiz -i grid
```

Full usage: 

```
Usage: ./TISCC [options...]
Options:
    -x, --dx               Single-tile code distance for X errors (Pauli weight of the minimum-weight logical X operator) (Required)
    -z, --dz               Single-tile code distance for Z errors (Pauli weight of the minimum-weight logical Z operator) (Required)
    -t, --dt               Code distance along the time dimension (number of surface code cycles for a round of error correction) (Required)
    -i, --info             Information to be queried. Options: {instructions, plaquettes, grid, parity}
    -o, --operation        Surface code operation to be compiled. Options: {idle, prepz, prepx, measz, measx, pauli_x, pauli_y, pauli_z, hadamard, inject_y, inject_t, flip_patch, move_right, swap_left, extension, contraction, move, merge, split, jointmeas, mergecontract, extendsplit, bellprep, bellmeas}
    -s, --tile_spec        Size of grid. Options: {single (default), double-vert, double-horiz}
    -d, --debug            Provide extra output for the purpose of debugging
    -p, --printgates       Print the time-resolved gate sequence for the operation given
    -r, --resources        Provide resource analysis for the operation given
    -h, --help             Shows this page 
```

### Example Output

The first layer of CNOT gates from a round of error-correction over a d = 2 patch:

```     
        Time (us)   Instruction  qsite1,qsite2
            0.0      Prepare_Z             36
            0.0      Prepare_Z             29
            0.0      Prepare_Z             43
           10.0        Y_-pi/4             29
           10.0        Y_-pi/4             43
           20.0         Z_pi/2             29
           20.0         Z_pi/2             43
           23.0        Y_-pi/4             36
           23.0        Y_-pi/4             40
           33.0           Move          36,37
           33.0           Move          43,44
           38.2           Move          37,34
           38.2           Move          44,41
          148.5             ZZ          33,34
          148.5             ZZ          40,41
         2148.5           Move          34,37
         2148.5           Move          41,44
         2258.8           Move          37,36
         2258.8           Move          44,43
         2264.0        Z_-pi/4             33
         2264.0        X_-pi/4             36
         2264.0        Z_-pi/4             43
         2264.0        X_-pi/4             40
         2274.0        Z_-pi/4             36
         2274.0        Z_-pi/4             40
```

Resource analysis and gate counts for one round of error correction over a d = 5 patch:

```
Grid area: 1.016e-04 m^2.
Computation time: 9.613e-03 s.
Space-time volume: 9.767e-07 s*m^2.
Trapping zones: 252
Trapping zone-seconds: 2.422e+00 zone*s
Trapping zone-seconds (active): 2.176e-01 zone*s
Active zone-seconds (%): 8.982
Junction: 160
Measure_Z: 24
Move: 160
Prepare_Z: 24
X_-pi/4: 80
Y_-pi/4: 104
ZZ: 80
Z_-pi/4: 160
Z_pi/2: 24
```

ASCII representation of the hardware grid allocated to two horizontally-adjacent dx = 5, dz = 3 surface code patches. The grid is made of three types of trapping zones (qsites): `M' (memory), `O' (operation), and `J' (junction). 

```
  0       1       2       3       4       5       6       7       8       9       10       11       
0 J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M 
  M       M       M       M       M       M       M       M       M       M       M       M       
  O       O       O       O       O       O       O       O       O       O       O       O       
  M       M       M       M       M       M       M       M       M       M       M       M       
1 J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M 
  M       M       M       M       M       M       M       M       M       M       M       M       
  O       O       O       O       O       O       O       O       O       O       O       O       
  M       M       M       M       M       M       M       M       M       M       M       M       
2 J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M 
  M       M       M       M       M       M       M       M       M       M       M       M       
  O       O       O       O       O       O       O       O       O       O       O       O       
  M       M       M       M       M       M       M       M       M       M       M       M       
3 J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M J M O M 
  M       M       M       M       M       M       M       M       M       M       M       M       
  O       O       O       O       O       O       O       O       O       O       O       O       
  M       M       M       M       M       M       M       M       M       M       M       M 
```

Another representation of the grid where all sites that are occupied by default are labeled 'O' except for the syndrome measurement qsites of a merged patch, which are instead labeled by their operator type.

```
  0       1       2       3       4       5       6       7       8       9       10       11       
0 - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - 
  |       |       |       |       |       |       |       |       |       |       |       |       
  O       O       Z       O       Z       O       Z       O       Z       O       Z       O       
  |       |       |       |       |       |       |       |       |       |       |       |       
1 - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - 
  |       |       |       |       |       |       |       |       |       |       |       |       
  X       Z       X       Z       X       Z       X       Z       X       Z       X       O       
  |       |       |       |       |       |       |       |       |       |       |       |       
2 - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - 
  |       |       |       |       |       |       |       |       |       |       |       |       
  O       X       Z       X       Z       X       Z       X       Z       X       Z       X       
  |       |       |       |       |       |       |       |       |       |       |       |       
3 - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - - - O - 
  |       |       |       |       |       |       |       |       |       |       |       |       
  O       Z       O       Z       O       Z       O       Z       O       Z       O       O       
  |       |       |       |       |       |       |       |       |       |       |       | 
```

## TISCC API

## Example Usage

Compile a Bell state preparation circuit through a joint XX measurement of two (vertically-adjacent) patches prepared in the |0> state:

```
// Designate the number rows and columns for a single tile
// Note: We add one or two buffer strips of qubits depending on whether the code distance is even or odd
unsigned int nrows = dz+1+!(dz%2); 
unsigned int ncols = dx+1+!(dx%2);

// Construct hardware grid with appropriate dimensions to hold two tiles top and bottom 
GridManager* grid = new TISCC::GridManager(2*nrows, ncols);

// Initialize logical qubit object using the grid
LogicalQubit* lq1 = new TISCC::LogicalQubit(dx, dz, 0, 0, *grid);

// Initialize second logical qubit object to the bottom of the first
LogicalQubit* lq2 = new TISCC::LogicalQubit(dx, dz, nrows, 0, *grid);

// Create a merged qubit
LogicalQubit* lq = lq1->get_merged_lq(*lq2, *grid);

// Initialize time tracker
double time = 0;

// Initialize hardware circuit
std::vector<TISCC::HW_Instruction> hw_master;

// Prepare both patches in the |0> state
lq1->transversal_op("prepz", *grid, hw_master, time);
lq2->transversal_op("prepz", *grid, hw_master, time);

// Perform merge operation
time = lq->merge(cycles, *grid, hw_master, time);

// Perform split operation
time = lq->split(*grid, hw_master, time);

// Resolve any junction conflicts
grid->enforce_hw_master_validity(hw_master);

// Count resources
grid->resource_counter(hw_master);
```

If this Bell preparation circuit were fed into a simulator (or real hardware), one might wish to perform a Pauli correction depending on the measured value of the joint XX operator. To aid in these types of situations, TISCC contains functions that help collect the qsite labels for which resulting measurement outcomes need to be multiplied.

```
// Collect measurement qsites for all X stabilizers lying between the two default-edge logical X operators
std::vector<unsigned int> deformation_qsites = lq->get_logical_deformation_operator_movement('X', nrows, *grid);
```

### Documentation

Documentation has been generated using doxygen. To view it, use your web browser to open the index.html within the docs/html directory.

## Regression Tests

TISCC uses the same regression test suite as [Lattice Surgery Compiler](https://github.com/latticesurgery-com/liblsqecc/tree/main). See the README.md file within the regression_tests directory for usage.

## Acknowledgements

TISCC was developed at Oak Ridge National Laboratory as part of the DARPA Quantum Benchmarking program. TISCC is one member of an end-to-end compilation and resource estimation pipeline for fault-tolerant quantum computing on trapped-ion processors. 