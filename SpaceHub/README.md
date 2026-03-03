# SpaceHub

[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fmnras%2Fstab1189-blue)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.1053W/abstract)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/std/the-standard)

SpaceHub is a high-performance, header-only C++ template library for extremely high precision few-body and large scale N-body simulations in astrophysics. It provides a comprehensive toolkit for gravitational dynamics simulations with support for various integration methods, particle interactions, and advanced features like tidal forces and scattering experiments.

## Features

- **Header-only library**: Easy integration into existing projects
- **High precision**: Optimized for extremely accurate few-body simulations
- **Multiple integrators**: Support for various numerical methods (IAS15, Bulirsch-Stoer, Gauss-Radau, etc.)
- **Advanced physics**: Newtonian gravity, post-Newtonian effects, tidal forces, drag forces
- **Flexible particle types**: Point particles, finite-size particles, tidal particles
- **Scattering experiments**: Built-in tools for orbital scattering simulations
- **Multi-threading support**: Parallel execution for large-scale simulations
- **Comprehensive callbacks**: Extensive output and monitoring capabilities
- **Convenient constants**: Built-in astronomical unit system

## Quick Start

### Requirements

- C++17 compatible compiler
- Optional: MPFR library for arbitrary precision arithmetic

### Basic Usage

```cpp
#include "src/spaceHub.hpp"
using namespace hub;
using namespace unit;

int main() {
    // Create particles
    using Solver = methods::DefaultMethod<>;
    using Particle = Solver::Particle;
    
    Particle star{1_Ms};  // 1 solar mass
    Particle planet{1_Me}; // 1 earth mass
    
    // Create a Keplerian orbit
    auto orbit = orbit::Elliptic(star.mass, planet.mass, 
                                1_AU, 0.2,        // semi-major axis, eccentricity
                                0_deg, 0_deg,     // inclination, longitude of ascending node
                                0_deg, 0_deg);    // argument of periapsis, true anomaly
    
    // Apply orbital motion to planet
    orbit::move_particles(orbit, planet);
    orbit::move_to_COM_frame(star, planet);
    
    // Create solver and run simulation
    Solver solver{0, star, planet};
    Solver::RunArgs args;
    args.add_stop_condition(1000_year);
    args.add_operation(callback::DefaultWriter("output.txt"));
    
    solver.run(args);
    return 0;
}
```

## Tutorial Examples

The `tutorial/` directory contains comprehensive examples organized by complexity:

### 1. Basic System Creation
- **1_1**: Direct initial conditions - Manual specification of particle states
- **1_2**: Keplerian orbits - Creating elliptical orbital configurations
- **1_3**: Random orbits - Generating random orbital parameters
- **1_4**: Planetary systems - Multi-planet configurations
- **1_5**: Hierarchical triples - Three-body hierarchical systems
- **1_6**: Hierarchical multipoles - Complex multi-body systems
- **1_7**: Many-body systems - Large particle ensembles

### 2. Simulation Control & Output
- **2_1**: Output callbacks - Basic data writing
- **2_2**: Time-sliced callbacks - Periodic output at specific intervals
- **2_3**: Lambda callbacks - Custom operations using lambda functions
- **2_4**: Start/end point callbacks - Operations at simulation boundaries
- **2_5**: Stop conditions - Custom termination criteria
- **2_6**: Variable capture - Advanced callback data handling
- **2_7**: Tick-tock timing - Performance monitoring

### 3. Advanced Physics & Solvers
- **3_1**: Alternative solvers - Using Bulirsch-Stoer and other integrators
- **3_2**: Additional forces - Post-Newtonian corrections
- **3_3**: Multiple forces - Combining different physical effects
- **3_4**: Finite-size particles - Extended body simulations
- **3_5**: Multiple solvers - Using different integrators in one simulation
- **3_6**: Tidal forces - Tidal interactions and tidal particles
- **3_7**: Drag forces - Atmospheric and gas drag effects

### 4. Scattering & Advanced Applications
- **4_1**: Hyperbolic orbits - Unbound orbital configurations
- **4_2**: Scattering tools - Basic scattering experiment setup
- **4_3**: Single-binary scattering - Two-body vs. single-body encounters
- **4_4**: Binary-binary scattering - Four-body scattering experiments
- **4_5**: Planetary system scattering - Complex multi-body encounters
- **4_6**: Multi-threading - Parallel execution for large experiments

## Building and Running

### Using Tutorial Makefile

```bash
cd tutorial
make all  # Compile all examples
./1_1direct-initial-conditions  # Run specific example
```

### Manual Compilation

```bash
g++ -std=c++17 -O3 -pthread your_program.cpp -o your_program
```

## Integration Methods

SpaceHub provides multiple high-precision integrators including IAS15, Bulirsch-Stoer, Gauss-Radau, and symplectic methods. 

Example solver selection:
```cpp
using Solver = methods::BS<>;          // Bulirsch-Stoer
using Solver = methods::GaussRadau<>;  // Gauss-Radau  
using Solver = methods::DefaultMethod<>; // Default method
```

## Physical Interactions

### Force Models
- **Newtonian gravity**: Standard gravitational interactions
- **Post-Newtonian**: Relativistic corrections (1PN, 2.5PN)
- **Tidal forces**: Static and dynamical tidal effects
- **Drag forces**: Gas drag and atmospheric effects

### Particle Types
- **Point particles**: Standard gravitational point masses
- **Finite-size particles**: Extended bodies with physical radii
- **Tidal particles**: Bodies with tidal deformation properties
- **Drag particles**: Particles experiencing drag forces

Example with tidal forces:
```cpp
using namespace force;
using f = Interactions<NewtonianGrav, Tidal>;
using Solver = methods::DefaultMethod<f, particles::TideParticles>;

// Create tidal particle: (mass, radius, apsidal_constant, lag_time)
Particle planet{1_Mj, 1_Rj, 0.25, 10_sec};
```

## Callbacks and Output

### Built-in Callbacks
- `DefaultWriter`: Standard particle state output
- `EnergyErrWriter`: Energy conservation monitoring
- `TimeStepWriter`: Integration step size logging

### Custom Callbacks
```cpp
// Lambda callback example
args.add_operation([&](auto& particles, auto step_size) {
    // Custom operation on particles
    print(std::cout, "Time: ", particles.time(), "\n");
});

// Conditional callbacks
args.add_step_slice_operation(100_year, [&](auto& ptc, auto h) {
    // Execute every 100 years
});
```

## Performance Optimization

### Multi-threading
```cpp
#include "src/taskflow/taskflow.hpp"

tf::Executor executor;
for (size_t i = 0; i < job_count; ++i) {
    executor.silent_async(simulation_job, i, parameters);
}
executor.wait_for_all();
```

### High-Precision Arithmetic
Enable MPFR support by compiling with MPFR flags:
```bash
g++ -std=c++17 -O3 -pthread -lmpfr -lgmp your_program.cpp -o your_program
```

## Documentation

- **API Documentation**: Generated with Doxygen
- **Website**: https://yihanwangastro.github.io/SpaceHub/
- **Paper**: [MNRAS 505, 1053 (2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.1053W/abstract)

## Citation

If you use SpaceHub in your research, please cite:

```bibtex
@article{SpaceHub2021,
    title={SpaceHub: A high-performance gravity integration toolkit for few-body problems in astrophysics},
    author={Wang, Yihan and Leigh, Nathan and Liu, Bin and Perna, Rosalba},
    journal={Monthly Notices of the Royal Astronomical Society},
    volume={505},
    number={1},
    pages={1053--1070},
    year={2021},
    doi={10.1093/mnras/stab1189}
}
```

## License

SpaceHub is licensed under the GNU General Public License v3.0. See [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## Acknowledgments

SpaceHub integrates several external libraries:
- **Taskflow**: Modern C++ parallel computing framework
- **Catch2**: Unit testing framework
- **MPFR**: Multiple precision floating-point library (optional)