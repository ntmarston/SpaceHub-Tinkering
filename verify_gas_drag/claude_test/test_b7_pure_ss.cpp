// Phase B7: Pure Shakura-Sunyaev disk (no Q=1 zone boundary)
//
// Hypothesis test: the 43% grad_P discontinuity at the Standard/Self-Reg
// boundary in disk_stae321_pagn.csv causes the BS integrator timestep
// collapse at ~172,000 yr. A pure S-S disk has no such boundary, so grad_P
// should be smooth throughout the orbit, and the simulation should complete.
//
// If PASS (<120s): hypothesis confirmed — zone boundary is the root cause.
// If TIMEOUT: something else also causes stiffness; check b7_pure_ss.log.

#include "../../SpaceHub/src/spaceHub.hpp"
#include <chrono>
#include <sstream>
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;

using f = Interactions<NewtonianGrav, DiskModel>;
using Solver = methods::DefaultMethod<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {

    DiskModel::init_from_file("../../SpaceHub/src/interaction/disk_tab/disk_pure_ss_1e8.csv");
    DiskModel::enable_dynamical_friction = true;
    DiskModel::enable_aerodynamic_drag   = false;
    DiskModel::enable_bondi_hoyle        = false;

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

    auto stop_time = 250000_year;   // must exceed 172,000 yr crash point
    args.add_stop_condition(stop_time);

    // Per-step diagnostic (mirrors test_diagnostic_long.cpp)
    size_t step_count = 0;
    auto wall_start = std::chrono::steady_clock::now();

    auto diagnostic = [&step_count, &wall_start](auto& ptc, auto step_size) {
        step_count++;
        bool print_step = (step_count <= 20) || (step_count % 1000 == 0);

        double t_sim = ptc.time();
        if (step_size < 1e-5 && step_count > 20) print_step = true;

        if (print_step) {
            auto wall_now  = std::chrono::steady_clock::now();
            double wall_sec = std::chrono::duration<double>(wall_now - wall_start).count();
            double t_years  = t_sim / (2.0 * M_PI);

            // Also log cylindrical R so we can see when orbit crosses old boundary (1791 AU)
            auto const& pos = ptc.pos();
            double dx = pos[1].x - pos[0].x;
            double dy = pos[1].y - pos[0].y;
            double R_cyl = std::sqrt(dx*dx + dy*dy);

            std::cerr << "step=" << step_count
                      << " t_yr=" << t_years
                      << " dt="   << step_size
                      << " R_cyl="<< R_cyl
                      << " wall=" << wall_sec << "s\n";
        }
    };

    args.add_operation(diagnostic);
    args.add_operation(TimeSlice(DefaultWriter("output/b7_pure_ss.dat"), 0_year, stop_time, 500));

    print(std::cout, "B7: Pure S-S disk (no Q=1 zone), dynfric only, 250kyr\n");
    print(std::cout, "    Loaded: disk_pure_ss_1e8.csv\n");
    print(std::cout, "    Orbit: M=1e8 Msun, m=30 Msun, sma=0.01 PC, e=0.67, i=5 deg\n");

    tools::Timer timer;
    timer.start();
    solver.run(args);
    print(std::cout << "B7: Complete in " << timer.get_time() << "s, " << step_count << " steps\n");

    return 0;
}
// Compile: g++ -std=c++17 -O3 -pthread test_b7_pure_ss.cpp -o test_b7
// Run:     timeout 120 ./test_b7 2>output/b7_pure_ss.log
