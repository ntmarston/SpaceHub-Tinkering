// Jimenez Type I migration smoketest — circular, coplanar orbit at 100 Rg
// Uses DiskMigration with migration_Jimenez only (JM17 linear torques)
// Expected: semi-major axis decreases over ~25 Myr timescale

#include "../../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;

using f = Interactions<NewtonianGrav, DiskMigration>;
using Solver = methods::Sym4<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main() {
    DiskMigration::init_from_file("../../src/interaction/disk_tab/SG_01Edd.csv");

    DiskMigration::LoweccDamping_CN06 = false;
    DiskMigration::CN06_ECC_DECOUPLED = false;
    DiskMigration::CN06_MIG_DECOUPLED = false;
    DiskMigration::migration_Jimenez = true;
    DiskMigration::inclined_zhu = false;

    Scalar m1 = 1e8_Ms;
    Scalar m2 = 30_Ms;
    Scalar r1 = 2.0 * m1 / (consts::C * consts::C);
    Scalar r2 = 2.0 * m2 / (consts::C * consts::C);
    Scalar sma = 100.0 * m1 / (consts::C * consts::C);

    Particle p1{m1, r1};
    Particle p2{m2, r2};

    auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, 0.0, 0_deg, 0_deg, 0_deg, 0_deg);
    orbit::move_particles(orb, p2);
    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};

    Solver::RunArgs args;
    args.rtol = 1e-9;

    auto stop_time = 1000000_year;
    args.add_stop_condition(stop_time);
    args.add_operation(TimeSlice(DefaultWriter("out/jimenez.dat"), 0_year, stop_time, 1000));

    tools::Timer timer;
    timer.start();
    solver.run(args);

    print(std::cout, "Jimenez test complete in ", timer.get_time(), "s\n");
    return 0;
}

// g++ -std=c++17 -O3 -pthread JimenezTest.cpp -o JimenezTest
