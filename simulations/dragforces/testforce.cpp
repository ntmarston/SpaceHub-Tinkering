
#include "../../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace force;  // save writting force::
// add whatever external term you like, just keep the first one to be the internal Newtonian gravity. check available force at xxx
using f = Interactions<NewtonianGrav, StaticGasField>;

using Solver = methods::DefaultMethod<f, particles::SizeParticles>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;

int main(int argc, char** argv) {
    auto M1 = 100_Ms;
    auto R1 = 0.00042450051_Rs;
    auto M2 = 1_Ms;
    auto R2 = 4.2450051e-6_Rs;
    Particle p1{M1, R1};
    Particle p2{M2, R2};

    auto orb = orbit::Elliptic(p1.mass, p2.mass, 3_AU, 0.0, 0_deg, 0_deg, 0_deg, 0_deg);

    orbit::move_particles(orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    //Solver: t_start, particle1, particle2, ...
    Solver solver{0, p1, p2};
    Solver::RunArgs args;
    args.add_stop_condition(100_year);

    args.add_operation(DefaultWriter("testforce.txt"));

    solver.run(args);

    print(std::cout, "Simulation Complete!\n");

    return 0;
}
