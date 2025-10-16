
#include "../../src/spaceHub.hpp"
#include "../../src/orbits/orbits.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace force;  
using f = Interactions<NewtonianGrav>;

using Solver = methods::DefaultMethod<f>;
using Scalar = Solver::Scalar;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;

int main(int argc, char** argv) {

    Particle p1{1000_Ms};
    Particle p2{1_Ms};


    //Elliptic parameters: M1, M2, initial: semimajor axis, eccentricity, inclination, longitude of ascending node, argument of periapsis, true anomaly
    auto orb = orbit::Elliptic(p1.mass, p2.mass, 1_AU, 0.0, hub::consts::pi, 0.0, 0.0_deg, 0_deg);
    //a=5AU, e = 0.1, i=45deg, \Omega=270deg, \omega = 90deg nu = 0deg
    orbit::move_particles(orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    //Solver: t_start, particle1, particle2, ...
    Solver solver{0, p1, p2};
    
    Solver::RunArgs args;

    args.add_stop_condition(1000_year);

    /*auto check_elements = [&p1, &p2](auto &ptc, auto &h){

        size_t particle_num = ptc.number();
        for (size_t i = 0; i < particle_num; ++i)
                {
                    for (size_t j = i + 1; j < particle_num; ++j)
                    {
                        auto dr = ptc.pos(i) - ptc.pos(j);
                        auto dv = ptc.vel(i) - ptc.vel(j);
                        auto m1 = p1.mass;
                        auto m2 = p2.mass;
                        
                        //auto L = orbit::calc_angular_momentum(m1, m2, dr, dv);
                        //print(std::cout << "L:" << L << '\n');
                    }
                }
        

    };*/

        auto check_elements = [&p1, &p2](auto& ptc, auto h) {
        size_t particle_num = ptc.number();
             for (size_t i = 0; i < particle_num; ++i)
                {
                    for (size_t j = i + 1; j < particle_num; ++j)
                    {
                        auto dr = ptc.pos(i) - ptc.pos(j);
                        auto dv = ptc.vel(i) - ptc.vel(j);
                        auto m1 = p1.mass;
                        auto m2 = p2.mass;
                        
                        auto L = orbit::calc_angular_momentum(m1, m2, dr, dv);
                        //print(std::cout << "L:" << L << '\n');
                        print(std::cout, "L:", L, '\n');
                    }
                }


        
        //dummy_variable += 1;
    };

    args.add_operation(StepSlice(check_elements, 1000));



    args.add_operation(DefaultWriter("simulation_results/Result.txt"));
    //args.add_operation(StepSlice(check_elements, 20));

    solver.run(args);

    print(std::cout, "Simulation complete");

    return 0;
}
