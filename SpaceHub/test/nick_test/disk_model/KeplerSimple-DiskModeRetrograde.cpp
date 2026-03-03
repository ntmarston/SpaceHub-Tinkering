
#include "../../../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
/*--------------------------------------------------New-----------------------------------------------------------*/
// Newtonian gravity + tabulated disk model gas drag
using f = Interactions<NewtonianGrav, DiskModel>;

using Solver = methods::DefaultMethod<f, particles::SizeParticles>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {

    DiskModel::init_from_file("SpaceHub/src/interaction/disk_tab/disk_test.csv");  // Must initialize before solver

    Scalar m1 = 1e8_Ms;  // DiskModel assumes particle 0 is the central mass
    Scalar m2 = 10_Ms;
    Scalar sma = 1000_AU;
    auto ecc = 0.5;
    auto inclination = 180_deg;  // Retrograde orbit test
    auto longitude_of_ascending_node = 0_deg;
    auto argument_of_periapsis = 0_deg; 
    auto true_anomaly = 0_deg;

    
    

    Particle p1{m1, 1e-10_AU};
    Particle p2{m2, 2e-7_AU};

    auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, ecc, inclination, longitude_of_ascending_node, argument_of_periapsis, true_anomaly);



    orbit::move_particles(orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};

    Solver::RunArgs args;
    args.rtol = 1e-10;

    auto collision_detect = [](auto &ptc, auto h)
            {
                //print(std::cout << "time:" << ptc.time() << "step:" << h << "\n");
                size_t particle_num = ptc.number();
                for (size_t i = 0; i < particle_num; ++i)
                {
                    for (size_t j = i + 1; j < particle_num; ++j)
                    {

                        if (distance(ptc.pos(i), ptc.pos(j)) < ptc.radius(i) + ptc.radius(j))
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
    auto stop_time = 1e6 * year;
    args.add_stop_condition(stop_time);

    args.add_operation(TimeSlice(DefaultWriter("SpaceHub/test/nick_test/disk_model/AnalyticalTests/DiskModel-Retrograde.dat"), 0_year, stop_time, 500));

    tools::Timer timer;
    timer.start();

    solver.run(args);

    double elapsed_time = timer.get_time();

    print(std::cout << "Simulation complete in " << elapsed_time << "s with no errors!\n");

    return 0;
}



// Commands to compile and run this simulation for copy+paste purposes:
// g++ -std=c++17 -O3 -pthread test/nick_test/KeplerSimple-DiskModeRetrograde.cpp -o test/nick_test/KeplerSimple-DiskModeRetrograde
// test/nick_test/KeplerSimple-DiskModeRetrograde