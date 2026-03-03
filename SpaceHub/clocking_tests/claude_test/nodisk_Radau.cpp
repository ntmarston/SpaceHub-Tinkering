
#include "../../src/spaceHub.hpp"
#include <sstream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
/*--------------------------------------------------New-----------------------------------------------------------*/
// Newtonian gravity + tabulated disk model — but ALL disk forces DISABLED
// Baseline: non-regularized Radau for comparison with AR_Radau in conservative case
using f = Interactions<NewtonianGrav, DiskModel>;

using Solver = methods::Radau<f, particles::SizeParticles>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {

    // Parse command-line arguments for inclination
    double inclination_deg = 5.0;  // Default value
    if (argc > 1) {
        inclination_deg = std::atof(argv[1]);
        print(std::cout << "Using inclination: " << inclination_deg << " degrees\n");
    } else {
        print(std::cout << "No inclination specified, using default: " << inclination_deg << " degrees\n");
        print(std::cout << "Usage: " << argv[0] << " <inclination_in_degrees>\n");
    }

    // Load pagn-generated disk file (Sirko-Goodman model, alpha=0.01, le=1.0)
    DiskModel::init_from_file("../../src/interaction/disk_tab/disk_default_pagn.csv");

    // DISABLE all disk forces — conservative gravity only
    DiskModel::enable_dynamical_friction = false;
    DiskModel::enable_aerodynamic_drag = false;
    DiskModel::enable_bondi_hoyle = false;

    // Central AGN black hole and orbiting 30 Msun star
    Scalar m1 = 1e8_Ms;  // DiskModel assumes particle 0 is the central mass
    Scalar m2 = 30_Ms;

    // Schwarzschild radius for central black hole: R_s = 2GM/c^2
    Scalar r1 = 2.0 * consts::G * m1 / (consts::C * consts::C);

    // Estimate using main-sequence mass-radius relation: R = R_sun * (M/M_sun)^(3/4)
    Scalar r2 = pow(m2/Ms, 0.75) * Rs;

    // Orbital parameters: 5000 AU (self-reg zone, Q=1 boundary at ~3306 AU), e=0.67
    Scalar sma = 5000_AU;
    auto ecc = 0.67;

    auto inclination = inclination_deg * 1_deg;
    auto longitude_of_ascending_node = 0_deg;
    auto argument_of_periapsis = 0_deg;  // cos(0) = +1, satisfies cos(omega_0) = ±1
    auto true_anomaly = 0_deg;

    // Physical radii for particles
    Particle p1{m1, r1};  // Central BH with Schwarzschild radius
    Particle p2{m2, r2};  // 30 Msun main-sequence star with M-R relation

    auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, ecc, inclination, longitude_of_ascending_node, argument_of_periapsis, true_anomaly);

    orbit::move_particles(orb, p2);
    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};

    Solver::RunArgs args;
    args.rtol = 1e-9;

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

    auto stop_time = 3e4_year;
    args.add_stop_condition(stop_time);

    // Build output and log filenames with inclination
    std::ostringstream output_filename;
    std::ostringstream log_filename;
    output_filename << "nodisk_Radau" << "-clocktest-" << inclination_deg << ".dat";
    log_filename << "nodisk_Radau" << "-clocktest-" << inclination_deg << ".log";
    args.add_operation(TimeSlice(DefaultWriter(output_filename.str()), 0_year, stop_time, 10000));

    // --- Logging: start time, 9 checkpoints at 10%-90%, final elapsed ---
    auto log_file = std::make_shared<std::ofstream>(log_filename.str());

    // Log START with system time
    auto start_sys_time = std::chrono::system_clock::now();
    auto start_tt = std::chrono::system_clock::to_time_t(start_sys_time);
    (*log_file) << "START " << std::put_time(std::localtime(&start_tt), "%Y-%m-%d %H:%M:%S") << std::endl;

    tools::Timer wall_timer;
    wall_timer.start();

    // 9 checkpoints at 10%, 20%, ..., 90% of stop_time
    auto checkpoint_logger = [&wall_timer, log_file](auto& ptc, auto step_size) {
        double wall_s = wall_timer.get_time();
        auto now = std::chrono::system_clock::now();
        auto now_tt = std::chrono::system_clock::to_time_t(now);
        (*log_file) << "CHECK sim_t=" << ptc.time() << " wall_t=" << wall_s << "s "
                    << std::put_time(std::localtime(&now_tt), "%H:%M:%S") << std::endl;
    };
    args.add_operation(TimeSlice(checkpoint_logger, stop_time * 0.1, stop_time, 9));

    solver.run(args);

    // Log STOP with total elapsed time
    double elapsed_time = wall_timer.get_time();
    (*log_file) << "STOP elapsed=" << elapsed_time << "s" << std::endl;

    print(std::cout << "Simulation complete in " << elapsed_time << "s with no errors!\n");

    return 0;
}

// Conservative (no-disk-forces) baseline Radau test
// Compile: g++ -std=c++17 -O3 -pthread nodisk_Radau.cpp -o ct_nodisk_Radau
// Run: ./ct_nodisk_Radau 5.0
