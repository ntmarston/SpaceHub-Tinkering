// CN08 2D sweep — eccentricity × inclination at 493042 AU
// Usage: ./sweep <eccentricity> <inclination_deg>

#include "../../src/spaceHub.hpp"
#include <sstream>
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;

using f = Interactions<NewtonianGrav, DiskMigration>;
using Solver = methods::Sym4<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {
    DiskMigration::init_from_file("../../src/interaction/disk_tab/SG_01Edd.csv");

    // Parse eccentricity
    double ecc_val = 0.1;
    if (argc > 1) {
        ecc_val = std::atof(argv[1]);
    } else {
        print(std::cout, "Usage: ", argv[0], " <eccentricity> <inclination_deg>\n");
        print(std::cout, "No eccentricity specified, using default: ", ecc_val, "\n");
    }

    // Parse inclination (degrees)
    double incl_deg = 0.0;
    if (argc > 2) {
        incl_deg = std::atof(argv[2]);
    } else {
        print(std::cout, "No inclination specified, using default: ", incl_deg, " deg\n");
    }

    print(std::cout, "ecc=", ecc_val, "  incl=", incl_deg, " deg\n");

    Scalar m1 = 1e8_Ms;
    Scalar m2 = 30_Ms;
    Scalar r1 = 2.0 * m1 / (consts::C * consts::C);
    Scalar r2 = 2.0 * m2 / (consts::C * consts::C);
    Scalar sma = 493042_AU;

    Particle p1{m1, r1};
    Particle p2{m2, r2};

    auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, ecc_val,
                                incl_deg * deg, 0_deg, 0_deg, 0_deg);
    orbit::move_particles(orb, p2);
    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};

    Solver::RunArgs args;
    args.rtol = 1e-9;

    auto stop_time = 1e7_year;
    args.add_stop_condition(stop_time);

    // Build output filename: CN08Sweep/ecc-{e}_incl-{i}.dat
    std::ostringstream output_filename;
    output_filename << "out/CN08Sweep/ecc-" << ecc_val << "_incl-" << incl_deg << ".dat";

    args.add_operation(TimeSlice(DefaultWriter(output_filename.str()), 0_year, stop_time, 1000));

    tools::Timer timer;
    timer.start();
    solver.run(args);

    print(std::cout, "CN08 sweep complete in ", timer.get_time(), "s\n");
    return 0;
}

// cd SpaceHub/test/migration_test && g++ -std=c++17 -O3 -pthread sweep.cpp -o bin/sweep
// ./bin/sweep 0.1 0.5
