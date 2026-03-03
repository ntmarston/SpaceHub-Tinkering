// Diagnostic with position output to track orbital state at timestep crash
// Short time slice (1000 years) to capture orbital elements near transition

#include "../../SpaceHub/src/spaceHub.hpp"
#include <sstream>
#include <chrono>
#include <cmath>
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

    auto stop_time = 200000_year;  // Past the 172k yr transition
    args.add_stop_condition(stop_time);

    // Per-step: print position of particle 2 relative to particle 1
    size_t step_count = 0;
    auto wall_start = std::chrono::steady_clock::now();

    // Open position log
    std::ofstream pos_log("output/position_log.txt");
    pos_log << "step,t_yr,dt,R_cyl,z,R3d\n";
    pos_log << std::scientific;

    auto diagnostic = [&](auto& ptc, auto step_size) {
        step_count++;
        double t_yr = ptc.time() / (2.0 * M_PI);

        // Always log near transition region, sample sparsely elsewhere
        bool do_log = false;
        if (t_yr > 169000 && t_yr < 175000) do_log = true;
        if (step_count % 10000 == 0) do_log = true;

        if (do_log) {
            auto const& p = ptc.pos();
            auto dr = p[1] - p[0];
            double R_cyl = std::sqrt(dr.x*dr.x + dr.y*dr.y);
            double z = dr.z;
            double R3d = std::sqrt(R_cyl*R_cyl + z*z);

            pos_log << step_count << "," << t_yr << "," << step_size << ","
                    << R_cyl << "," << z << "," << R3d << "\n";

            if (step_count <= 20 || step_count % 50000 == 0 || (t_yr > 169000 && step_count % 1000 == 0)) {
                auto wall_now = std::chrono::steady_clock::now();
                double wall_s = std::chrono::duration<double>(wall_now - wall_start).count();
                std::cerr << "step=" << step_count << " t_yr=" << t_yr
                          << " dt=" << step_size << " R_cyl=" << R_cyl
                          << " R3d=" << R3d << " wall=" << wall_s << "s\n";
            }
        }
    };

    args.add_operation(diagnostic);

    // Fast output: every 10,000 years (not 30,000) to capture orbital evolution
    args.add_operation(TimeSlice(DefaultWriter("output/pos_diag.dat"), 0_year, stop_time, 20000));

    print(std::cout, "Position diagnostic: i=5deg, stop=200kyr\n");

    tools::Timer timer;
    timer.start();
    solver.run(args);
    print(std::cout << "Complete in " << timer.get_time() << "s, steps=" << step_count << "\n");

    return 0;
}
