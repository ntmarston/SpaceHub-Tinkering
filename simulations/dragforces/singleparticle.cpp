
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
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {

    std::fstream energy_file("simulations/dragforces/energylog.txt", std::ios::out);


    auto M2 = 1_Ms;
    auto R2 = 4.2450051e-6_Rs;
    Particle p1{M2, R2, 0, 0, 0, 25_kms, 0, 0};
    Particle p2{1_kg, 0.01_AU, 100,100,100,0,0,0};
     auto orb = orbit::Elliptic(p1.mass, p2.mass, 3_AU, 0.0, 0_deg, 0_deg, 0_deg, 0_deg);
    orbit::move_particles(orb, p2);

    auto kick = [](auto& ptc, auto step_size){
            auto const &m = ptc.mass();
            auto const &p = ptc.pos();
            auto const &v = ptc.vel();
            v[0] = Vec3(10_kms, 0_kms, 0_kms);
            print(std::cout << "Kicked\n");
    };
    
    //Solver: t_start, particle1, particle2, ...
    Solver solver{0, p1, p2};
    Solver::RunArgs args;
    //args.add_start_point_operation(kick);
    args.add_stop_condition(100_year);

    auto twriter = TimeSlice(DefaultWriter("sp.txt"), 0.0, 100_year, 50);
    args.add_operation(DefaultWriter("FdynSubsonicTest.txt"));
    

    solver.run(args);

    print(std::cout, "Simulation Complete!\n");

    return 0;
}
//g++ -std=c++17 -O3 -pthread simulations/dragforces/singleparticle.cpp -o simulations/dragforces/singleparticle
//simulations/dragforces/singleparticle