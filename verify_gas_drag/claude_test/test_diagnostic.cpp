// Phase A: Diagnostic test to reproduce and characterize the DiskModel hang
// Mirrors stae321-12-pagn.cpp setup but with per-step stderr logging
// and short stop_time (10 yr) + 2-minute wall-clock timeout

#include "../../SpaceHub/src/spaceHub.hpp"
#include <sstream>
#include <chrono>
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;

using f = Interactions<NewtonianGrav, DiskModel>;
using Solver = methods::DefaultMethod<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {

    // Load disk table
    DiskModel::init_from_file("../../SpaceHub/src/interaction/disk_tab/disk_stae321_pagn.csv");
    DiskModel::enable_dynamical_friction = true;
    DiskModel::enable_aerodynamic_drag = false;
    DiskModel::enable_bondi_hoyle = false;

    // Orbit parameters (identical to stae321-12-pagn.cpp)
    Scalar m1 = 1e8_Ms;
    Scalar m2 = 30_Ms;
    Scalar r1 = 2.0 * consts::G * m1 / (consts::C * consts::C);
    Scalar sma = 0.01_PC;
    auto ecc = 0.67;
    auto inclination = 5_deg;
    Scalar Rs_30Msun = 2.0 * consts::G * m2 / (consts::C * consts::C);

    Particle p1{m1, r1};
    Particle p2{m2, Rs_30Msun};

    auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, ecc, inclination, 0_deg, 0_deg, 0_deg);
    orbit::move_particles(orb, p2);
    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};
    Solver::RunArgs args;
    args.rtol = 1e-10;

    auto stop_time = 10_year;
    args.add_stop_condition(stop_time);

    // Per-step diagnostic: print step count, sim time, and step size to stderr
    size_t step_count = 0;
    auto wall_start = std::chrono::steady_clock::now();

    auto diagnostic = [&step_count, &wall_start](auto& ptc, auto step_size) {
        step_count++;
        // Print every 100 steps to avoid I/O bottleneck
        if (step_count % 100 == 0 || step_count <= 10) {
            auto wall_now = std::chrono::steady_clock::now();
            double wall_sec = std::chrono::duration<double>(wall_now - wall_start).count();

            // Get simulation time in years (1 year = 2*pi in SpaceHub units)
            double t_sim = ptc.time();
            double t_years = t_sim / (2.0 * M_PI);

            std::cerr << "step=" << step_count
                      << " t=" << t_sim
                      << " t_yr=" << t_years
                      << " dt=" << step_size
                      << " wall=" << wall_sec << "s"
                      << "\n";
        }
    };

    args.add_operation(diagnostic);

    // Also write output file for comparison
    args.add_operation(TimeSlice(DefaultWriter("output/diagnostic.dat"), 0_year, stop_time, 500));

    print(std::cout, "Starting Phase A diagnostic: i=5deg, stop=10yr\n");
    print(std::cout, "Watch stderr (diagnostic.log) for per-step output\n");

    tools::Timer timer;
    timer.start();
    solver.run(args);
    double elapsed = timer.get_time();

    print(std::cout << "Completed in " << elapsed << "s, " << step_count << " steps\n");

    return 0;
}

// Compile: g++ -std=c++17 -O3 -pthread test_diagnostic.cpp -o diag
// Run:     timeout 120 ./diag 2>output/diagnostic.log
