
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
// Pure Newtonian gravity — NO DiskModel in force template at all
// Tests whether AR_Radau works correctly without vel_dependent forces
using f = Interactions<NewtonianGrav>;

using Solver = methods::AR_Radau<f, particles::SizeParticles>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {

    double inclination_deg = 5.0;
    if (argc > 1) {
        inclination_deg = std::atof(argv[1]);
        print(std::cout << "Using inclination: " << inclination_deg << " degrees\n");
    } else {
        print(std::cout << "No inclination specified, using default: " << inclination_deg << " degrees\n");
    }

    Scalar m1 = 1e8_Ms;
    Scalar m2 = 30_Ms;
    Scalar r1 = 2.0 * consts::G * m1 / (consts::C * consts::C);
    Scalar r2 = pow(m2/Ms, 0.75) * Rs;
    Scalar sma = 5000_AU;
    auto ecc = 0.67;
    auto inclination = inclination_deg * 1_deg;

    Particle p1{m1, r1};
    Particle p2{m2, r2};

    auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, ecc, inclination, 0_deg, 0_deg, 0_deg);
    orbit::move_particles(orb, p2);
    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};
    Solver::RunArgs args;
    args.rtol = 1e-9;

    auto collision_detect = [](auto &ptc, auto h) {
        size_t n = ptc.number();
        for (size_t i = 0; i < n; ++i)
            for (size_t j = i + 1; j < n; ++j)
                if (distance(ptc.pos(i), ptc.pos(j)) < ptc.radius(i) + ptc.radius(j))
                    return true;
        return false;
    };
    args.add_stop_condition(collision_detect);

    auto stop_time = 3e4_year;
    args.add_stop_condition(stop_time);

    std::ostringstream output_filename, log_filename;
    output_filename << "puregrav_AR_Radau-clocktest-" << inclination_deg << ".dat";
    log_filename << "puregrav_AR_Radau-clocktest-" << inclination_deg << ".log";
    args.add_operation(TimeSlice(DefaultWriter(output_filename.str()), 0_year, stop_time, 10000));

    auto log_file = std::make_shared<std::ofstream>(log_filename.str());
    auto start_sys_time = std::chrono::system_clock::now();
    auto start_tt = std::chrono::system_clock::to_time_t(start_sys_time);
    (*log_file) << "START " << std::put_time(std::localtime(&start_tt), "%Y-%m-%d %H:%M:%S") << std::endl;

    tools::Timer wall_timer;
    wall_timer.start();

    auto checkpoint_logger = [&wall_timer, log_file](auto& ptc, auto step_size) {
        double wall_s = wall_timer.get_time();
        auto now = std::chrono::system_clock::now();
        auto now_tt = std::chrono::system_clock::to_time_t(now);
        (*log_file) << "CHECK sim_t=" << ptc.time() << " wall_t=" << wall_s << "s "
                    << std::put_time(std::localtime(&now_tt), "%H:%M:%S") << std::endl;
    };
    args.add_operation(TimeSlice(checkpoint_logger, stop_time * 0.1, stop_time, 9));

    solver.run(args);

    double elapsed_time = wall_timer.get_time();
    (*log_file) << "STOP elapsed=" << elapsed_time << "s" << std::endl;
    print(std::cout << "Simulation complete in " << elapsed_time << "s with no errors!\n");

    return 0;
}
// Pure gravity AR_Radau test (no DiskModel template parameter)
// Compile: g++ -std=c++17 -O3 -pthread puregrav_AR_Radau.cpp -o ct_puregrav_AR_Radau
