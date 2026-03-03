// Fine-grained diagnostic: print position every step around the crash (steps 94k-115k)
// To find exact orbital state when timestep crashes at step ~102,000

#include "../../SpaceHub/src/spaceHub.hpp"
#include <fstream>
#include <cmath>
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;

using f = Interactions<NewtonianGrav, DiskModel>;
using Solver = methods::DefaultMethod<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main() {
    DiskModel::init_from_file("../../SpaceHub/src/interaction/disk_tab/disk_stae321_pagn.csv");
    DiskModel::enable_dynamical_friction = true;
    DiskModel::enable_aerodynamic_drag = false;
    DiskModel::enable_bondi_hoyle = false;

    Scalar m1 = 1e8_Ms, m2 = 30_Ms;
    Scalar r1 = 2.0 * consts::G * m1 / (consts::C * consts::C);
    Scalar Rs_30Msun = 2.0 * consts::G * m2 / (consts::C * consts::C);

    Particle p1{m1, r1}, p2{m2, Rs_30Msun};
    auto orb = orbit::Elliptic(m1, m2, 0.01_PC, 0.67, 5_deg, 0_deg, 0_deg, 0_deg);
    orbit::move_particles(orb, p2);
    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};
    Solver::RunArgs args;
    args.rtol = 1e-10;

    size_t step_count = 0;
    const double Rmin_val = 11.8621;

    std::ofstream fg("output/finegrained.txt");
    fg << "step,t_yr,dt_sp,R_cyl,z,R3d,in_disk\n";
    fg << std::scientific;

    // Log every step from step 94000 onward, stop at step 115000
    auto per_step_and_stop = [&](auto& ptc, auto step_size) -> bool {
        step_count++;
        if (step_count >= 94000) {
            auto const& p = ptc.pos();
            auto dr = p[1] - p[0];
            double R_cyl = std::sqrt(dr.x*dr.x + dr.y*dr.y);
            double z = dr.z;
            double R3d = std::sqrt(R_cyl*R_cyl + z*z);
            int in_disk = (R_cyl > Rmin_val) ? 1 : 0;
            double t_yr = ptc.time() / (2.0 * M_PI);
            fg << step_count << "," << t_yr << "," << step_size << ","
               << R_cyl << "," << z << "," << R3d << "," << in_disk << "\n";

            if (step_count % 1000 == 0) {
                std::cout << "step=" << step_count << " t_yr=" << t_yr
                          << " dt=" << step_size << " R_cyl=" << R_cyl << "\n";
                std::cout.flush();
            }
        }
        return step_count >= 115000;  // Stop after 115k steps
    };
    args.add_stop_condition(per_step_and_stop);

    print(std::cout, "Fine-grained: steps 94k-115k\n");
    tools::Timer timer; timer.start();
    solver.run(args);
    fg.flush();
    print(std::cout << "Done in " << timer.get_time() << "s, steps=" << step_count << "\n");
    return 0;
}
