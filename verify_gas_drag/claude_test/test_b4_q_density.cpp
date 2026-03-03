// Phase B4: DiskModel with Q-dependent vertical density restored
// Tests whether always-Gaussian density (replacing Q-dependent logic) causes the hang

#define TEST_Q_DENSITY
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

    auto stop_time = 250000_year;  // Must exceed 172,000 yr where timestep crashes
    args.add_stop_condition(stop_time);

    args.add_operation(TimeSlice(DefaultWriter("output/b4_q_density.dat"), 0_year, stop_time, 500));

    print(std::cout, "B4: Starting (Q-dependent vertical density)\n");
    tools::Timer timer;
    timer.start();
    solver.run(args);
    print(std::cout << "B4: Complete in " << timer.get_time() << "s\n");

    return 0;
}
// g++ -std=c++17 -O3 -pthread test_b4_q_density.cpp -o test_b4
// timeout 120 ./test_b4
