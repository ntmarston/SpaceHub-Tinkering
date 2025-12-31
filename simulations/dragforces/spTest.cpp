
#include "../../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace force;  

using f = Interactions<NewtonianGrav, StaticGasField>;

using Solver = methods::BS<f, particles::SizeParticles>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {

    //std::fstream energy_file("simulations/dragforces/energylog.txt", std::ios::out);

    auto M1 = 1_Ms;
    auto R1 = 4.2450051e-6_Rs;

    Particle p1{M1, R1, 0,0,0,25_kms,0,0};
    std::vector<Particle> particles;
    
    particles.push_back(p1);
    Solver solver{0, particles};
    Solver::RunArgs args;
    
    auto stop_time = 70_year;
    args.add_stop_condition(stop_time);

    //auto twriter = TimeSlice(DefaultWriter("dynamic_subsonic.txt"), 0.0, stop_time, 100);
    auto swriter = StepSlice(DefaultWriter("dynamic_subsonic.txt"), 1);
    args.add_operation(swriter);

    

    solver.run(args);

    print(std::cout, "Simulation Complete!\n");

    return 0;
}
// g++ -std=c++17 -O3 -pthread simulations/dragforces/spTest.cpp -o simulations/dragforces/spTest
// simulations/dragforces/spTest