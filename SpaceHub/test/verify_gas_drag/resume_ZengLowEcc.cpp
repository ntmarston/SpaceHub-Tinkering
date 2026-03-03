
#include "../../src/spaceHub.hpp"
#include "../../src/nick-tools.hpp"
#include <sstream>
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
/*--------------------------------------------------New-----------------------------------------------------------*/
// Resume timed-out ZengLowEcc-pagn simulations from their last snapshot.
// Physics identical to ZengLowEcc-pagn.cpp: e_ini=0.3, p_ini=300 Rg, dynamical friction only.
using f = Interactions<NewtonianGrav, DiskModel>;

using Solver = methods::DefaultMethod<f, particles::SizeParticles>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {

    // Parse command-line arguments for inclination
    double inclination_deg = 0;  // Default value
    if (argc > 1) {
        inclination_deg = std::atof(argv[1]);
        print(std::cout << "Using inclination: " << inclination_deg << " degrees\n");
    } else {
        print(std::cout << "No inclination specified, using default: " << inclination_deg << " degrees\n");
        print(std::cout << "Usage: " << argv[0] << " <inclination_in_degrees>\n");
    }

    // Load pagn-generated disk file (Sirko-Goodman, le=0.5, alpha=0.1)
    DiskModel::init_from_file("../../src/interaction/disk_tab/disk_ZengAndPan_pagn.csv");

    // Configure force toggles: ONLY dynamical friction enabled
    DiskModel::enable_dynamical_friction = true;
    DiskModel::enable_aerodynamic_drag = false;
    DiskModel::enable_bondi_hoyle = false;

    // Read last snapshot from the original (timed-out) output file
    std::ostringstream input_path;
    input_path << "ZengLowEcc-pagn/incl-" << static_cast<int>(inclination_deg) << ".dat";

    print(std::cout << "Reading last snapshot from: " << input_path.str() << "\n");

    auto state = hub::nick::read_last_snapshot_com<Particle>(input_path.str());

    print(std::cout << "Resume time: " << state.resume_time << "\n");
    print(std::cout << "Number of particles: " << state.particles.size() << "\n");

    // Initialize solver from resumed state
    Solver solver{state.resume_time, state.particles};

    Solver::RunArgs args;
    args.rtol = 1e-10;

    // Same collision detection as original
    auto collision_detect = [](auto &ptc, auto h)
            {
                size_t particle_num = ptc.number();
                for (size_t i = 0; i < particle_num; ++i)
                {
                    for (size_t j = i + 1; j < particle_num; ++j)
                    {
                        if (distance(ptc.pos(i), ptc.pos(j)) < ptc.radius(i) + ptc.radius(j))
                        {
                            return true;
                        }
                    }
                }
                return false;
            };

    args.add_stop_condition(collision_detect);

    // Same sma shrink stop: a < 0.5 * 325 AU (original initial sma)
    Scalar sma = 325_AU;
    Scalar m1 = 1e8_Ms;
    Scalar m2 = 30_Ms;

    auto sma_shrink_stop = [sma, m1, m2](auto &ptc, auto h) {
        auto dr = ptc.pos(1) - ptc.pos(0);
        auto dv = ptc.vel(1) - ptc.vel(0);
        auto r = norm(dr);
        auto v2 = dot(dv, dv);
        auto eps = 0.5 * v2 - consts::G * (m1 + m2) / r;
        auto a_cur = -consts::G * (m1 + m2) / (2.0 * eps);
        return a_cur < 0.5 * sma;
    };
    args.add_stop_condition(sma_shrink_stop);

    // Same stop time as original
    auto stop_time = 1.561e8_year;
    args.add_stop_condition(stop_time);

    // Output to NEW file with -resumed suffix
    std::ostringstream output_path;
    output_path << "ZengLowEcc-pagn/incl-" << static_cast<int>(inclination_deg) << "-resumed.dat";

    args.add_operation(TimeSlice(DefaultWriter(output_path.str()), state.resume_time, stop_time, 10000));

    tools::Timer timer;
    timer.start();

    solver.run(args);

    double elapsed_time = timer.get_time();

    print(std::cout << "Resumed simulation complete in " << elapsed_time << "s with no errors!\n");

    return 0;
}


// Commands to compile and run this simulation (run from verify_gas_drag/):
// g++ -std=c++17 -O3 -pthread resume_ZengLowEcc.cpp -o sim-resume-ZengLowEcc
// ./sim-resume-ZengLowEcc <inclination_in_degrees>
// Example: ./sim-resume-ZengLowEcc 45
