
#include "../../../../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
/*--------------------------------------------------New-----------------------------------------------------------*/
// using the Newtonian gravity + first order Post-Newtonian correction
using f = Interactions<NewtonianGrav, DiskMigration>;

using Solver = methods::Sym6<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;


int main(int argc, char** argv) {
    DiskMigration::init_from_file("SpaceHub/src/interaction/disk_tab/SG_01Edd.csv");

    DiskMigration::DISABLE_MIGRATION          = true;
    DiskMigration::DISABLE_I_DAMPING          = true;
    DiskMigration::DISABLE_DYNAMICAL_FRICTION = false;
    DiskMigration::DISABLE_AERODYNAMIC_DRAG   = false;
    DiskMigration::DISABLE_BONDI_HOYLE        = false;
    DiskMigration::ignore_dynamical_friction  = false;
    DiskMigration::DISABLE_E_DAMPING          = true;
    
        
    
    
    Scalar m1 = 1e8_Ms;
    Scalar m2 = 10_Ms;

    Scalar Rg  = consts::G * m1 / (consts::C * consts::C);  // gravitational radius of central mass
    Scalar sma = 1e5 * Rg;
    auto ecc = 0.2;
    auto inclination = 0_deg;
    auto longitude_of_ascending_node = 0_deg;
    auto argument_of_periapsis = 0_deg; 
    auto true_anomaly = 0_deg;

    
    

    Scalar r_s1 = 2.0 * consts::G * m1 / (consts::C * consts::C);  // Schwarzschild radius of m1
    Scalar r_s2 = 2.0 * consts::G * m2 / (consts::C * consts::C);  // Schwarzschild radius of m2

    Particle p1{m1, r_s1};
    Particle p2{m2, r_s2};

    auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, ecc, inclination, longitude_of_ascending_node, argument_of_periapsis, true_anomaly);



    orbit::move_particles(orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};

    Solver::RunArgs args;
    args.rtol = 1e-9;

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


    //auto Torb = period(orb);
    //Scalar orbital_endtime = 250 * Torb;
    //args.add_stop_condition(orbital_endtime);
    auto stop_time = 10000_year;
    args.add_stop_condition(stop_time);

    args.add_operation(TimeSlice(DefaultWriter("SpaceHub/test/migration_test/AnalyticalValidation/CN08/out/test_eccdamp.dat"), 0_year, stop_time, 500));

    tools::Timer timer;
    timer.start();

    solver.run(args);

    double elapsed_time = timer.get_time();

    print(std::cout << "Simulation complete in " << elapsed_time << "s with no errors!\n");

    return 0;
}



// Commands to compile and run this simulation for copy+paste purposes (run from project root):
// g++ -std=c++17 -O3 -pthread SpaceHub/test/migration_test/AnalyticalValidation/CN08/test_eccdamp.cpp -o SpaceHub/test/migration_test/AnalyticalValidation/CN08/test_eccdamp
// SpaceHub/test/migration_test/AnalyticalValidation/CN08/test_eccdamp