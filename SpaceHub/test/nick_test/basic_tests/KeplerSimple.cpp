
#include "../../../src/spaceHub.hpp"
#include <algorithm>
#include <cctype>
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
    double mratio = 1.0;      // m2/m1 = q
    double m1_val = 1.0;       // in solar masses
    double m2_val = -1.0;      // in solar masses (-1 means use q instead)
    double sma_val = 1.0;      // in AU
    double ecc = 0.0;
    double inc_val = 0.0;      // inclination deg
    double lan_val = 0.0;      // longitude of ascending node deg
    double aop_val = 0.0;      // argument of periapsis deg
    double ta_val = 0.0;       // true anomaly deg
    double atol_val = -1.0;    // absolute tolerance (-1 means use SpaceHub default)
    std::string outfile = "SpaceHub/test/nick_test/basic_tests/KeplerSimple.dat";
    std::string rtype = "ms";  // radius type: "bh" (Schwarzschild) or "ms" (main sequence)
    double stop_time = -1.0;   // stop time in years (-1 means use 250 orbital periods)

    // deal w/ args 
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
        } else if ((arg == "-atol" || arg == "--atol") && i + 1 < argc) {
            atol_val = std::stod(argv[++i]);
        } else if ((arg == "-out" || arg == "--out") && i + 1 < argc) {
            outfile = argv[++i];
        } else if ((arg == "-rtype" || arg == "--rtype") && i + 1 < argc) {
            rtype = argv[++i];
            // Convert to lowercase for case-insensitive comparison
            std::transform(rtype.begin(), rtype.end(), rtype.begin(),
                           [](unsigned char c) { return std::tolower(c); });
        } else if ((arg == "-stop" || arg == "--stop") && i + 1 < argc) {
            stop_time = std::stod(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]\n"
                      << "Options:\n"
                      << "  -m1 <value>   Primary mass in solar masses (default: 1.0)\n"
                      << "  -m2 <value>   Secondary mass in solar masses (overrides -q)\n"
                      << "  -q <value>    Mass ratio m2/m1 (default: 1.0)\n"
                      << "  -sma <value>  Semi-major axis in AU (default: 1.0)\n"
                      << "  -ecc <value>  Eccentricity (default: 0.0)\n"
                      << "  -inc <value>  Inclination in degrees (default: 0.0)\n"
                      << "  -lan <value>  Longitude of ascending node in degrees (default: 0.0)\n"
                      << "  -aop <value>  Argument of periapsis in degrees (default: 0.0)\n"
                      << "  -ta <value>   True anomaly in degrees (default: 0.0)\n"
                      << "  -atol <value> Absolute tolerance (default: SpaceHub default)\n"
                      << "  -out <file>   Output file (default: test/nick_test/KeplerSimple.dat)\n"
                      << "  -rtype <type> Radius type: 'bh' (Schwarzschild) or 'ms' (main sequence, default)\n"
                      << "  -stop <value> Stop time in years (default: 250 orbital periods)\n"
                      << "\nDefault (no args): Circular orbit of two 1 solar mass bodies at 1 AU separation.\n";
            return 0;
        }
    }

    // if m2 is still the default value -1 (not set by args), then calculate it from q
    if (m2_val < 0) {
        m2_val = mratio * m1_val;
    }

    // Validate radius type
    if (rtype != "bh" && rtype != "ms") {
        std::cerr << "Error: Invalid radius type '" << rtype << "'. Use 'bh' or 'ms'.\n";
        return 1;
    }

    // Calculate radii based on type
    // Schwarzschild radius: Rs = 2GM/c^2 ≈ 2.95 km per solar mass ≈ 4.24e-6 R_sun per solar mass
    // Main sequence approximation: R/R_sun ≈ (M/M_sun)^0.8 for M < 1 M_sun
    //                              R/R_sun ≈ (M/M_sun)^0.57 for M > 1 M_sun
    double r1_val, r2_val;  // in solar radii
    if (rtype == "bh") {
        const double Rs_per_Msun = 4.24e-6;  // Schwarzschild radius in solar radii per solar mass
        r1_val = Rs_per_Msun * m1_val;
        r2_val = Rs_per_Msun * m2_val;
    } else {  // "ms"
        // Main sequence mass-radius relation
        auto ms_radius = [](double m) {
            if (m <= 1.0) {
                return std::pow(m, 0.8);
            } else {
                return std::pow(m, 0.57);
            }
        };
        r1_val = ms_radius(m1_val);
        r2_val = ms_radius(m2_val);
    }

    // debug
    std::cout << "Running Simulation with:\n"
              << "  m1  = " << m1_val << " Ms\n"
              << "  m2  = " << m2_val << " Ms\n"
              << "  q   = " << (m2_val / m1_val) << "\n"
              << "  r1  = " << r1_val << " Rs (" << rtype << ")\n"
              << "  r2  = " << r2_val << " Rs (" << rtype << ")\n"
              << "  sma = " << sma_val << " AU\n"
              << "  ecc = " << ecc << "\n"
              << "  inc = " << inc_val << " deg\n"
              << "  lan = " << lan_val << " deg\n"
              << "  aop = " << aop_val << " deg\n"
              << "  ta  = " << ta_val << " deg\n"
              << "  stop= " << (stop_time >= 0 ? std::to_string(stop_time) + " yr" : "250 orbital periods") << "\n"
              << "  atol= " << (atol_val >= 0 ? std::to_string(atol_val) : "default") << "\n"
              << "  out = " << outfile << "\n";

    Scalar m1 = m1_val * 1_Ms;
    Scalar m2 = m2_val * 1_Ms;
    Scalar r1 = r1_val * 1_Rs;
    Scalar r2 = r2_val * 1_Rs;
    Scalar sma = sma_val * 1_AU;
    auto inclination = inc_val * 1_deg;
    auto longitude_of_ascending_node = lan_val * 1_deg;
    auto argument_of_periapsis = aop_val * 1_deg;
    auto true_anomaly = ta_val * 1_deg;

    Particle p1{m1, r1};
    Particle p2{m2, r2};

    auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, ecc, inclination, longitude_of_ascending_node, argument_of_periapsis, true_anomaly);

    orbit::move_particles(orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};

    Solver::RunArgs args;

    if (atol_val >= 0) {
        args.atol = atol_val;
    }

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

    // Set stop time based on user input or default to 250 orbital periods
    Scalar orbital_endtime;
    if (stop_time >= 0) {
        orbital_endtime = stop_time * 1_year;
    } else {
        auto Torb = period(orb);
        orbital_endtime = 250 * Torb;
    }
    args.add_stop_condition(orbital_endtime);

    args.add_operation(TimeSlice(DefaultWriter(outfile), 0_year, stop_time, 300));

    tools::Timer timer;
    timer.start();

    solver.run(args);

    double elapsed_time = timer.get_time();
    

    print(std::cout << "Simulation complete in " << elapsed_time << "s with no errors!\n");

    return 0;
}



// Commands to compile and run this simulation for copy+paste purposes:
// g++ -std=c++17 -O3 -pthread test/nick_test/KeplerSimple.cpp -o test/nick_test/KeplerSimple
// test/nick_test/KeplerSimple -ARGUMENTS
//
// Default simulation (circular orbit, 2x 1 Ms bodies at 1 AU, black hole radii):
// test/nick_test/KeplerSimple -rtype bh
//
// Custom example:
// test/nick_test/KeplerSimple -m1 1.0 -q 10 -sma 0.1 -ecc 0.1 -incl 60 -loa 0 -aop 0 -ta 0 -out "SpaceHub/test/nick_test/basic_tests/KeplerSimpleCurrent.dat"