
#include "../../src/spaceHub.hpp"
#include <sstream>
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
/*--------------------------------------------------New-----------------------------------------------------------*/
// Newtonian gravity + tabulated disk model with dynamical friction ONLY
// Uses pagn (Sirko-Goodman) disk model: le=0.5, alpha=0.1
// Reproduces Zeng & Pan fig 10: e_ini=0.7, p_ini=300 Rg
using f = Interactions<NewtonianGrav, DiskMigration>;

using Solver = methods::Sym4<f, particles::SizeParticles>;
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

    // Load pagn-generated disk file
    DiskMigration::init_from_file("../../src/interaction/disk_tab/SG_01Edd.csv");

    // Configure force toggles: ONLY dynamical friction enabled
    DiskMigration::DISABLE_MIGRATION          = true;
    DiskMigration::DISABLE_E_DAMPING          = true;
    DiskMigration::DISABLE_I_DAMPING          = true;

    DiskMigration::DISABLE_DYNAMICAL_FRICTION = false;
    
    DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
    DiskMigration::DISABLE_BONDI_HOYLE        = true;
    DiskMigration::ignore_dynamical_friction  = false;

    // Central AGN black hole and orbiting 30 Msun black hole
    Scalar m1 = 1e8_Ms;  // DiskModel assumes particle 0 is the central mass
    Scalar m2 = 30_Ms;   // Stellar-mass black hole

    Scalar r1 = 2.0 * consts::G * m1 / (consts::C * consts::C);

    // Orbital parameters: p_ini = 300 Rg, e=0.7
    // a = p/(1-e^2) = 295.6/0.51 ≈ 580 AU
    Scalar sma = 580_AU;
    auto ecc = 0.7;

    auto inclination = inclination_deg * 1_deg;
    auto longitude_of_ascending_node = 0_deg;
    auto argument_of_periapsis = 0_deg;
    auto true_anomaly = 0_deg;

    // Schwarzschild radius for 30 Msun BH
    Scalar Rs_30Msun = 2.0 * consts::G * m2 / (consts::C * consts::C);

    Particle p1{m1, r1};            // Central BH
    Particle p2{m2, Rs_30Msun};     // 30 Msun BH

    auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, ecc, inclination, longitude_of_ascending_node, argument_of_periapsis, true_anomaly);

    orbit::move_particles(orb, p2);
    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};

    Solver::RunArgs args;
    args.rtol = 1e-9;

    std::string exit_reason = "time_limit";

    auto collision_detect = [&exit_reason](auto &ptc, auto h)
            {
                size_t particle_num = ptc.number();
                for (size_t i = 0; i < particle_num; ++i)
                {
                    for (size_t j = i + 1; j < particle_num; ++j)
                    {
                        if (distance(ptc.pos(i), ptc.pos(j)) < ptc.radius(i) + ptc.radius(j))
                        {
                            exit_reason = "collision";
                            return true;
                        }
                    }
                }
                return false;
            };

    args.add_stop_condition(collision_detect);

    // Stop when semi-major axis has shrunk by 50% from its initial value (a < 0.5 * sma_ini)
    auto sma_shrink_stop = [sma, m1, m2, &exit_reason](auto &ptc, auto h) {
        auto dr = ptc.pos(1) - ptc.pos(0);
        auto dv = ptc.vel(1) - ptc.vel(0);
        auto r = norm(dr);
        auto v2 = dot(dv, dv);
        auto eps = 0.5 * v2 - consts::G * (m1 + m2) / r;
        auto a_cur = -consts::G * (m1 + m2) / (2.0 * eps);
        if (a_cur < 0.5 * sma) { exit_reason = "sma_shrink"; return true; }
        return false;
    };
    args.add_stop_condition(sma_shrink_stop);

    std::ofstream logfile;
    tools::Timer wall_timer;

    auto write_log = [&logfile, &wall_timer](auto &ptc, auto h) {
        logfile << ptc.time() << "," << h << "," << wall_timer.get_time() << "\n";
    };

    // t_stop = 1e15 M_bullet in geometrized units (G=c=1), converted to years:
    // T_M = G*M/c^3 = 4.9255e-6 s * 1e8 = 492.55 s = 1.561e-5 yr
    // => 1e15 * 1.561e-5 yr = 1.561e10 yr
    auto stop_time = 1e9_year;
    args.add_stop_condition(stop_time);

    // Build output filename with inclination
    std::ostringstream output_filename;
    std::ostringstream log_filename;
    output_filename << "out/ZengHighEcc/incl-" << inclination_deg << ".dat";
    log_filename << "out/ZengHighEcc/incl-" << inclination_deg << ".log";

    logfile.open(log_filename.str());
    logfile << "time,step_size,system_time\n";
    wall_timer.start();

    args.add_operation(TimeSlice(DefaultWriter(output_filename.str()), 0_year, stop_time, 10000));
    args.add_operation(TimeSlice(write_log, 0_year, stop_time, 10000));

    tools::Timer timer;
    timer.start();

    solver.run(args);

    double elapsed_time = timer.get_time();

    print(std::cout << "Exit condition: " << exit_reason << "\n");
    print(std::cout << "Simulation complete in " << elapsed_time << "s with no errors!\n");

    return 0;
}



// Commands to compile and run this simulation (starting from project root):
// cd SpaceHub/test/verify_gas_drag
// mkdir -p simulations/bin out/ZengHighEcc
// g++ -std=c++17 -O3 -pthread ZengHighEcc.cpp -o simulations/bin/ZengHighEcc
// ./simulations/bin/ZengHighEcc <inclination_in_degrees>
// Example: ./simulations/bin/ZengHighEcc 20.0
