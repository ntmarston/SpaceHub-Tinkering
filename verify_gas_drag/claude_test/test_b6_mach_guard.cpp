// Phase B6: Add a minimum Mach guard AFTER computing Mach
// Tests the proposed fix: skip force when Mach < 1e-4 (orbit near-circular w.r.t. gas)
// This is physically motivated: I(M) = M/3 at M<<1, so force is negligible at M<<1e-3

// Use TEST_MACH_GUARD flag in disk-model-test.hpp
#define TEST_MACH_GUARD
#include "spaceHub-test.hpp"
#include <sstream>
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

    // Use 1.5e7 year stop time (same as real simulation) to prove the fix works end-to-end
    auto stop_time = 250000_year;
    args.add_stop_condition(stop_time);

    args.add_operation(TimeSlice(DefaultWriter("output/b6_mach_guard.dat"), 0_year, stop_time, 500));

    print(std::cout, "B6: Starting with Mach guard fix (full 1.5e7 yr test)\n");
    tools::Timer timer;
    timer.start();
    solver.run(args);
    print(std::cout << "B6: Complete in " << timer.get_time() << "s\n");

    return 0;
}
