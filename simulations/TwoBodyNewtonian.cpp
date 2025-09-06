
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace force;  
using f = Interactions<NewtonianGrav>;

using Solver = methods::DefaultMethod<f>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;

int main(int argc, char** argv) {

    Particle p1{1000_Ms};
    Particle p2{1_Ms};


    //Elliptic parameters: M1, M2, initial: semimajor axis, eccentricity, inclination, longitude of ascending node, argument of periapsis, true anomaly
    auto orb = orbit::Elliptic(p1.mass, p2.mass, 5_AU, 0.3, 45_deg, 60_deg, 70_deg, 270_deg);
    //a=5AU, e = 0.1, i=45deg, \Omega=270deg, \omega = 90deg nu = 0deg
    orbit::move_particles(orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    //Solver: t_start, particle1, particle2, ...
    Solver solver{0, p1, p2};
    
    Solver::RunArgs args;

    args.add_stop_condition(1000_year);


    args.add_operation(DefaultWriter("simulation_results/Result.txt"));

    solver.run(args);

    print(std::cout, "Simulation complete");

    return 0;
}
