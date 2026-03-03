
#include "../../src/spaceHub.hpp"
#include <sstream>
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
/*--------------------------------------------------New-----------------------------------------------------------*/
// Newtonian gravity + tabulated disk model with dynamical friction ONLY
// Uses pagn (Sirko-Goodman) disk model instead of disktab
using f = Interactions<NewtonianGrav, DiskModel>;

using Solver = methods::DefaultMethod<f, particles::SizeParticles>;
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
    DiskModel::init_from_file("../../src/interaction/disk_tab/disk_stae321_pagn.csv");

    // Configure force toggles: ONLY dynamical friction enabled
    DiskModel::enable_dynamical_friction = true;
    DiskModel::enable_aerodynamic_drag = false;
    DiskModel::enable_bondi_hoyle = false;

    // Central AGN black hole and orbiting 30 Msun black hole
    Scalar m1 = 1e8_Ms;  // DiskModel assumes particle 0 is the central mass
    Scalar m2 = 30_Ms;   // Stellar-mass black hole

    Scalar r1 = 2.0 * consts::G * m1 / (consts::C * consts::C);
    // Orbital parameters: 0.01 pc, e=0.67, omega=0 (cos(omega)=+1)
    Scalar sma = 0.01_PC;
    auto ecc = 0.67;

    auto inclination = inclination_deg * 1_deg;
    auto longitude_of_ascending_node = 0_deg;
    auto argument_of_periapsis = 0_deg;  // cos(0) = +1, satisfies cos(omega_0) = ±1
    auto true_anomaly = 0_deg;

    // Physical radii for particles
    // Schwarzschild radius: Rs = 2GM/c^2
    // For 30 Msun: Rs = 2 * G * 30Ms / c^2
    Scalar Rs_30Msun = 2.0 * consts::G * m2 / (consts::C * consts::C);

    Particle p1{m1, r1};       // Central BH (effectively point mass)
    Particle p2{m2, Rs_30Msun};      // 30 Msun BH with Schwarzschild radius

    auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, ecc, inclination, longitude_of_ascending_node, argument_of_periapsis, true_anomaly);

    orbit::move_particles(orb, p2);
    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};

    Solver::RunArgs args;
    args.rtol = 1e-10;

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

    auto stop_time = 1.5e7_year;
    args.add_stop_condition(stop_time);

    // Build output filename with inclination
    std::ostringstream output_filename;
    output_filename << "fig12-pagn/DynFriction-i" << inclination_deg << ".dat";

    args.add_operation(TimeSlice(DefaultWriter(output_filename.str()), 0_year, stop_time, 500));

    tools::Timer timer;
    timer.start();

    solver.run(args);

    double elapsed_time = timer.get_time();

    print(std::cout << "Simulation complete in " << elapsed_time << "s with no errors!\n");
    print(std::cout << "30 Msun Schwarzschild radius: " << Rs_30Msun << " AU\n");

    return 0;
}



// Commands to compile and run this simulation (run from verify_gas_drag/):
// g++ -std=c++17 -O3 -pthread stae321-12-pagn.cpp -o stae321-12-pagn
// ./stae321-12-pagn [inclination_in_degrees]
// Example: ./stae321-12-pagn 10.0
