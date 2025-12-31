
#include "../../src/spaceHub.hpp"
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
    Particle p1{1_Ms, 1_Rs};
    Particle p2{1_Ms, 1_Rs};

    Scalar sma = 1_AU;
    auto ecc = 0.7;
    auto inclination = 0_deg;
    auto longitude_of_ascending_node = 0_deg;
    auto argument_of_periapsis = 0_deg; 
    auto true_anomaly = 0_deg;

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

                            print(std::cout << "collision" << "\n");
                           
                            return true;
                        }
                    }
                }
                return false;
            };

    args.add_stop_condition(collision_detect);

    args.add_stop_condition(1000_year);

    args.add_operation(StepSlice(DefaultWriter(("test/nick_test/KeplerSimple.dat")), 25));


    

  




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