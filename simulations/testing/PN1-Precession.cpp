
#include "../../SpaceHub/src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace force;  // save writting force::
// add whatever external term you like, just keep the first one to be the internal Newtonian gravity. check available force at xxx
using f = Interactions<NewtonianGrav, PN1>;

using Solver = methods::DefaultMethod<f>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;

int main(int argc, char** argv) {
    Particle p1{1_Ms};
    Particle p2{1_Ms};

    auto inner_orb = orbit::Elliptic(p1.mass, p2.mass, 1_AU, 0.95, 1_deg, 1_deg, 1_deg, 1_deg);

    orbit::move_particles(inner_orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    //Solver: t_start, particle1, particle2, ...
    Solver solver{0, p1, p2};

    Solver::RunArgs args;

    args.add_stop_condition(10000_year);

    auto t_writer = TimeSlice(DefaultWriter("simulations/testing/results/Result.txt"), 0.0, 10000_year, 1000);

    //args.add_operation(DefaultWriter("simulations/testing/results/Result.txt"));

    
    args.add_operation(t_writer);
    solver.run(args);

    print(std::cout, "Simulation Complete!\n");

    return 0;
}
