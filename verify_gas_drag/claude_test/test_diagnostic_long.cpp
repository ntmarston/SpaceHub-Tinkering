// Phase A (extended): Reproduce actual hang with 1.5e7 yr stop time
// Per-step diagnostic to stderr, 2-minute wall-clock timeout

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

    DiskModel::init_from_file("../../SpaceHub/src/interaction/disk_tab/disk_stae321_pagn.csv");
    DiskModel::enable_dynamical_friction = true;
    DiskModel::enable_aerodynamic_drag = false;
    DiskModel::enable_bondi_hoyle = false;

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

    // SAME stop_time as the real simulation
    auto stop_time = 1.5e7_year;
    args.add_stop_condition(stop_time);

    // Per-step diagnostic
    size_t step_count = 0;
    auto wall_start = std::chrono::steady_clock::now();

    auto diagnostic = [&step_count, &wall_start](auto& ptc, auto step_size) {
        step_count++;
        // Print every 1000 steps, and the first 20, and whenever dt drops below some threshold
        bool print_step = (step_count <= 20) || (step_count % 1000 == 0);

        // Also print if step_size is very small (potential stiffness indicator)
        double t_sim = ptc.time();
        if (step_size < 1e-5 && step_count > 20) print_step = true;

        if (print_step) {
            auto wall_now = std::chrono::steady_clock::now();
            double wall_sec = std::chrono::duration<double>(wall_now - wall_start).count();
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
    args.add_operation(TimeSlice(DefaultWriter("output/diagnostic_long.dat"), 0_year, stop_time, 500));

    print(std::cout, "Diagnostic (long): i=5deg, stop=1.5e7yr — same as real sim\n");

    tools::Timer timer;
    timer.start();
    solver.run(args);
    print(std::cout << "Completed in " << timer.get_time() << "s, " << step_count << " steps\n");

    return 0;
}
// g++ -std=c++17 -O3 -pthread test_diagnostic_long.cpp -o diag_long
// timeout 120 ./diag_long 2>output/diagnostic_long.log
