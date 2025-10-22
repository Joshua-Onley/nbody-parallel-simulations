# N-body Parallel Simulation

A high-performance N-body gravitational simulation implementing serial, OpenMP, and MPI parallelisation strategies for performance comparison and analysis.

## Overview

This project simulates the gravitational interactions between N bodies in 3D space using Newton's law of universal gravitation. The simulation implements four different approaches:

- **Serial implementation**: Baseline single-threaded execution for performance benchmarking
- **OpenMP implementation**: Shared-memory parallelisation using multi-threading
- **MPI implementation**: Distributed-memory parallelisation across multiple processes
- **Visualisation implementation**: Real-time graphical display using SFML library

The project was developed as part of a High-Performance Computing (HPC) module to explore and compare different parallelisation paradigms, with an additional visualisation component for educational demonstration.

## Physics Model

The simulation uses a direct N-body algorithm where every body interacts with every other body through gravitational forces:

```
F = G * (m₁ * m₂) / r²
```

Where:
- `G` is the gravitational constant (normalised to 1.0)
- `m₁` and `m₂` are the masses of two bodies
- `r` is the distance between them

The system evolves using a simple Euler integration scheme with time step `dt = 0.1`.

## Implementation Details

### Serial Version (`serial.cpp`)

The baseline implementation uses a nested loop structure with O(N²) complexity:
- Outer loop: iterates through all bodies
- Inner loop: calculates pairwise forces with remaining bodies
- Newton's third law is exploited to reduce calculations by half

### OpenMP Version (`OpenMP.cpp`)

Parallelises the force calculation loop using:
- Thread-local force arrays to avoid race conditions
- `#pragma omp parallel` for thread creation
- `#pragma omp for schedule(static, 2)` for work distribution
- `#pragma omp critical` section for force accumulation
- Configurable thread count (default: 2 threads)

### MPI Version (`MPI.cpp`)

**Note**: The current MPI implementation is identical to the OpenMP version and requires proper MPI parallelisation to be added.

### Visualisation Version (`visualisation.cpp`)

Provides real-time graphical rendering of the N-body simulation using SFML (Simple and Fast Multimedia Library):
- Displays bodies as coloured circles in a 800×600 window
- Uses random red-tinted colours for visual distinction
- Maps 3D positions to 2D screen coordinates (projects x-y plane)
- Smaller time step (`dt = 0.001`) for smoother animation
- Adjustable animation speed via sleep delay (default: 10ms per frame)
- Implements forward Euler integration with symplectic properties
- Increased body count to 900 for more dramatic visual effects

## Simulation Parameters

### Performance Testing Versions (Serial, OpenMP, MPI)

| Parameter | Value | Description |
|-----------|-------|-------------|
| N | 840 | Number of bodies |
| T | 2000.0 | Total simulation time |
| dt | 0.1 | Time step size |
| G | 1.0 | Gravitational constant (normalised) |

### Visualisation Version

| Parameter | Value | Description |
|-----------|-------|-------------|
| N | 900 | Number of bodies |
| T | 2000.0 | Total simulation time |
| dt | 0.001 | Time step size (smaller for smoother animation) |
| G | 1.0 | Gravitational constant (normalised) |
| Window Size | 800×600 | Display resolution |
| Particle Radius | 5 pixels | Visual size of each body |

### Initial Conditions

All quantities are randomly generated from normal distributions:
- **Masses**: μ = 1.0, σ = 0.0 (uniform mass)
- **Positions**: μ = 0.0, σ = 1.0 (centred distribution)
- **Velocities**: μ = 0.0, σ = 1.0 (random initial velocities)

## Compilation

### Serial Version
```bash
g++ -O3 -o serial serial.cpp
```

### OpenMP Version
```bash
g++ -O3 -fopenmp -o openmp OpenMP.cpp
```

### MPI Version
```bash
mpic++ -O3 -fopenmp -o mpi MPI.cpp
```

### Visualisation Version
```bash
g++ -O3 -o visualisation visualisation.cpp -lsfml-graphics -lsfml-window -lsfml-system
```

**Note**: SFML library must be installed for the visualisation version. See [SFML installation guide](https://www.sfml-dev.org/tutorials/)

## Execution

### Serial
```bash
./serial
```

### OpenMP
```bash
./openmp
```

To change the number of threads, modify `num_threads` in the source code or set the environment variable:
```bash
export OMP_NUM_THREADS=4
./openmp
```

### MPI
```bash
mpirun -np 4 ./mpi
```

### Visualisation
```bash
./visualisation
```

The visualisation window will open automatically. Close the window to terminate the simulation. To adjust animation speed, modify the `sleep_for` milliseconds value in the source code (line ~145).

## Performance Analysis

The simulation outputs execution time in seconds. Key metrics for comparison:

- **Speedup**: S(p) = T₁ / Tₚ
- **Efficiency**: E(p) = S(p) / p
- **Parallel Overhead**: Communication and synchronisation costs

### Expected Performance Characteristics

- **Serial**: O(N²) per time step, baseline reference
- **OpenMP**: Scales with thread count, limited by memory bandwidth and critical section overhead
- **MPI**: Requires substantial communication for position/force data distribution

## Known Issues and Limitations

1. **MPI Implementation**: Currently identical to OpenMP version - proper MPI parallelisation needed
2. **Force Calculation Bug**: In OpenMP/MPI versions, `local_force[j]` subtracts force (Newton's 3rd law) but should add the opposite force
3. **Critical Section Bottleneck**: The force accumulation critical section can limit OpenMP scalability
4. **Visualisation Performance**: The visualisation version is not optimised for performance (single-threaded with rendering overhead)
5. **2D Projection**: Visualisation only displays x-y coordinates; z-axis depth not represented
6. **Numerical Integration**: Simple Euler method; more accurate schemes (e.g., Verlet, RK4) would improve accuracy
7. **Collision Handling**: No softening parameter or collision detection implemented

## Future Improvements

- Implement proper MPI domain decomposition
- Add Barnes-Hut or Fast Multipole Method for O(N log N) complexity
- Enhance visualisation with 3D rendering (OpenGL/Vulkan)
- Add interactive camera controls and zoom functionality
- Use reduction clauses instead of critical sections in OpenMP
- Add energy conservation checks for validation
- Implement more sophisticated integration schemes (Verlet, leapfrog, RK4)
- Add gravitational softening to prevent singularities
- Implement trajectory trails in visualisation
- Add ability to save/load simulation states

## Requirements

### All Versions
- C++11 compatible compiler
- Standard C++ libraries

### Parallel Versions
- OpenMP support (gcc 4.9+, clang 3.8+, or equivalent)
- MPI library (OpenMPI, MPICH, etc.)

### Visualisation Version
- SFML 2.5+ (Simple and Fast Multimedia Library)
  - Graphics module
  - Window module
  - System module

## Licence

This project was developed for educational purposes as part of an HPC module.

## Author

Submitted as coursework for High-Performance Computing module.

---

**Computational Complexity**: O(N² × T/dt) per simulation  
**Memory Usage**: O(N) for positions, velocities, forces, and masses