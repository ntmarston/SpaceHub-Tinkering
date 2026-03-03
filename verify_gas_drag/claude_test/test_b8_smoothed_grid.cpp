// Phase B8: Smoothed ZengAndPan disk — exhaustive parameter grid
//
// Tests whether gradient smoothing at the Q=1 zone boundary resolves the
// Bulirsch-Stoer timestep collapse for a range of (sma, ecc, inc) values.
//
// Usage: ./test_b8 <sma_PC> <ecc> <inc_deg>
//   Example: ./test_b8 0.01 0.67 5
//
// Disk: disk_ZengAndPan_pagn.csv (smoothed, Q=1 boundary at R/Rg ~ 3350)
// Forces: dynamical friction only
// Duration: 250 kyr
//
// PASS: completes within 120s wall-clock with no dt collapse
// FAIL: timeout or dt collapse at boundary

#include "../../SpaceHub/src/spaceHub.hpp"
#include <chrono>
#include <sstream>
#include <cstdlib>
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;

using f = Interactions<NewtonianGrav, DiskModel>;
using Solver = methods::DefaultMethod<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <sma_PC> <ecc> <inc_deg>\n";
        return 1;
    }

    double sma_pc  = std::atof(argv[1]);
    double ecc_val = std::atof(argv[2]);
    double inc_deg = std::atof(argv[3]);

    DiskModel::init_from_file("../../SpaceHub/src/interaction/disk_tab/disk_ZengAndPan_pagn.csv");
    DiskModel::enable_dynamical_friction = true;
    DiskModel::enable_aerodynamic_drag   = false;
    DiskModel::enable_bondi_hoyle        = false;

    Scalar m1 = 1e8_Ms;
    Scalar m2 = 30_Ms;
    Scalar r1 = 2.0 * consts::G * m1 / (consts::C * consts::C);
    Scalar sma = sma_pc * 1_PC;
    Scalar inclination = inc_deg * 1_deg;
    Scalar Rs_30Msun = 2.0 * consts::G * m2 / (consts::C * consts::C);

    Particle p1{m1, r1};
    Particle p2{m2, Rs_30Msun};

    auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, ecc_val, inclination, 0_deg, 0_deg, 0_deg);
    orbit::move_particles(orb, p2);
    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};
    Solver::RunArgs args;
    args.rtol = 1e-10;

    auto stop_time = 250000_year;
    args.add_stop_condition(stop_time);

    // Build output filename from parameters
    std::ostringstream dat_name;
    dat_name << "output/b8_a" << sma_pc << "_e" << ecc_val << "_i" << inc_deg << ".dat";

    // Per-step diagnostic logging to stderr
    size_t step_count = 0;
    double min_dt = 1e30;
    auto wall_start = std::chrono::steady_clock::now();

    auto diagnostic = [&](auto& ptc, auto step_size) {
        step_count++;
        if (step_size < min_dt) min_dt = step_size;

        bool print_step = (step_count <= 20) || (step_count % 1000 == 0);
        if (step_size < 1e-5 && step_count > 20) print_step = true;

        if (print_step) {
            auto wall_now = std::chrono::steady_clock::now();
            double wall_sec = std::chrono::duration<double>(wall_now - wall_start).count();
            double t_years = ptc.time() / (2.0 * M_PI);

            auto const& pos = ptc.pos();
            double dx = pos[1].x - pos[0].x;
            double dy = pos[1].y - pos[0].y;
            double R_cyl = std::sqrt(dx * dx + dy * dy);

            std::cerr << "step=" << step_count
                      << " t_yr=" << t_years
                      << " dt=" << step_size
                      << " R_cyl=" << R_cyl
                      << " wall=" << wall_sec << "s\n";
        }
    };

    args.add_operation(diagnostic);
    args.add_operation(TimeSlice(DefaultWriter(dat_name.str()), 0_year, stop_time, 500));

    std::cout << "B8: sma=" << sma_pc << "PC e=" << ecc_val << " i=" << inc_deg
              << "deg | ZengAndPan_pagn (smoothed)\n";

    tools::Timer timer;
    timer.start();
    solver.run(args);

    double elapsed = timer.get_time();
    std::cout << "B8: Complete in " << elapsed << "s, " << step_count
              << " steps, min_dt=" << min_dt << "\n";

    return 0;
}
// Compile: g++ -std=c++17 -O3 -pthread test_b8_smoothed_grid.cpp -o test_b8
// Run:     timeout 120 ./test_b8 0.01 0.67 5 2>output/b8_a0.01_e0.67_i5.log
