
#include "../../../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
/*--------------------------------------------------New-----------------------------------------------------------*/
// using the Newtonian gravity + first order Post-Newtonian correction
using f = Interactions<NewtonianGrav, StaticGasField>;

using Solver = methods::DefaultMethod<f, particles::SizeParticles>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {

    // Default values
    double mratio = 10.0;      // m2/m1 = q
    double m1_val = 1.0;       // solar masses
    double m2_val = -1.0;      // flagged to use q instead
    double sma_val = 0.1;      // AU
    double ecc = 0.1;
    double inc_val = 0.0;      // inclination deg
    double lan_val = 0.0;      // longitude of ascending node deg
    double aop_val = 0.0;      // argument of periapsis deg
    double ta_val = 0.0;       // true anomaly in deg
    //std::string outfile = "SpaceHub/test/nick_test/basic_tests/KeplerSimple.dat";

    // cmd arg processing
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if ((arg == "-m1" || arg == "--m1") && i + 1 < argc) {
            m1_val = std::stod(argv[++i]);
        } else if ((arg == "-m2" || arg == "--m2") && i + 1 < argc) {
            m2_val = std::stod(argv[++i]);
        } else if ((arg == "-q" || arg == "--q") && i + 1 < argc) {
            mratio = std::stod(argv[++i]);
        } else if ((arg == "-sma" || arg == "--sma") && i + 1 < argc) {
            sma_val = std::stod(argv[++i]);
        } else if ((arg == "-ecc" || arg == "--ecc") && i + 1 < argc) {
            ecc = std::stod(argv[++i]);
        } else if ((arg == "-inc" || arg == "--inc") && i + 1 < argc) {
            inc_val = std::stod(argv[++i]);
        } else if ((arg == "-lan" || arg == "--lan") && i + 1 < argc) {
            lan_val = std::stod(argv[++i]);
        } else if ((arg == "-aop" || arg == "--aop") && i + 1 < argc) {
            aop_val = std::stod(argv[++i]);
        } else if ((arg == "-ta" || arg == "--ta") && i + 1 < argc) {
            ta_val = std::stod(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]\n"
                      << "Options:\n"
                      << "  -m1 <value>   Primary mass in solar masses (default: 1.0)\n"
                      << "  -m2 <value>   Secondary mass in solar masses (overrides -q)\n"
                      << "  -q <value>    Mass ratio m2/m1 (default: 10.0)\n"
                      << "  -sma <value>  Semi-major axis in AU (default: 0.1)\n"
                      << "  -ecc <value>  Eccentricity (default: 0.1)\n"
                      << "  -inc <value>  Inclination in degrees (default: 0.0)\n"
                      << "  -lan <value>  Longitude of ascending node in degrees (default: 0.0)\n"
                      << "  -aop <value>  Argument of periapsis in degrees (default: 0.0)\n"
                      << "  -ta <value>   True anomaly in degrees (default: 0.0)\n";
            return 0;
        }
    }

    
    if (m2_val < 0) {
        m2_val = mratio * m1_val;
    }

    // debug
    std::cout << "Running Simulation with:\n"
              << "  m1  = " << m1_val << " Ms\n"
              << "  m2  = " << m2_val << " Ms\n"
              << "  q   = " << (m2_val / m1_val) << "\n"
              << "  sma = " << sma_val << " AU\n"
              << "  ecc = " << ecc << "\n"
              << "  inc = " << inc_val << " deg\n"
              << "  lan = " << lan_val << " deg\n"
              << "  aop = " << aop_val << " deg\n"
              << "  ta  = " << ta_val << " deg\n";

    Scalar m1 = m1_val * 1_Ms;
    Scalar m2 = m2_val * 1_Ms;
    Scalar sma = sma_val * 1_AU;
    auto inclination = inc_val * 1_deg;
    auto longitude_of_ascending_node = lan_val * 1_deg;
    auto argument_of_periapsis = aop_val * 1_deg;
    auto true_anomaly = ta_val * 1_deg;

    Particle p1{m1, 1_Rs};
    Particle p2{m2, 1_Rs};

    auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, ecc, inclination, longitude_of_ascending_node, argument_of_periapsis, true_anomaly);

    orbit::move_particles(orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};

    Solver::RunArgs args;

    auto collision_detect = [](auto &ptc, auto h)
            {
                //print(std::cout << "time:" << ptc.time() << "step:" << h << "\n");
                size_t particle_num = ptc.number();
                for (size_t i = 0; i < particle_num; ++i)
                {
                    for (size_t j = i + 1; j < particle_num; ++j)
                    {

                        if (distance(ptc.pos(i), ptc.pos(j)) < ptc.radius(i) + ptc.radius(i))
                        {
                            return true;
                        }
                    }
                }
                return false;
            };

    args.add_stop_condition(collision_detect);


    auto Torb = period(orb);
    Scalar orbital_endtime = 250 * Torb;
    args.add_stop_condition(orbital_endtime);
    args.add_stop_condition(1000_year);

    //args.add_operation(StepSlice(DefaultWriter(outfile), 10));

    tools::Timer timer;
    timer.start();

    solver.run(args);

    double elapsed_time = timer.get_time();

    print(std::cout << "Simulation complete in " << elapsed_time << "s with no errors!\n");

    return 0;
}



// Commands to compile and run this simulation for copy+paste purposes:
// g++ -std=c++17 -O3 -pthread test/nick_test/KeplerSimple.cpp -o test/nick_test/KeplerSimple
// test/nick_test/KeplerSimple