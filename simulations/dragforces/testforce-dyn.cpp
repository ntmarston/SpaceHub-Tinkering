
#include "../../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace force;  // save writting force::
// add whatever external term you like, just keep the first one to be the internal Newtonian gravity. check available force at xxx
using f = Interactions<StaticGasField>;

using Solver = methods::Sym6<f, particles::SizeParticles>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {

    std::fstream energy_file("simulations/dragforces/energylog.txt", std::ios::out);

    auto M1 = 1_Ms;
    //auto R1 = 0.00042450051_Rs;
    auto R1 = 12_km;
    auto M2 = 0.001_Ms;
    auto R2 = 4.2450051e-6_Rs;
    Particle p1{M1, R1, 0,0,0,100_kms,0,0};
    Particle p2{M2, R2, 0, 100, 0, 0, 0, 0};


    //auto orb = orbit::Elliptic(p1.mass, p2.mass, 100_AU, 0.0, 0_deg, 0_deg, 0_deg, 0_deg);

    //orbit::move_particles(orb, p2);

    //orbit::move_to_COM_frame(p1, p2);
    /*
    auto energy_log = [&energy_file ](auto& ptc, auto step_size) {
            auto const &m = ptc.mass();
            auto const &p = ptc.pos();
            auto const &v = ptc.vel();

            size_t particle_num = ptc.number();
                for (size_t i = 0; i < particle_num; ++i)
                {
                    for (size_t j = i + 1; j < particle_num; ++j)
                    {
                        auto iid = i;
                        auto jid = j;
                        decltype(p[iid]) dr = p[iid] - p[jid];
                        Scalar potential_eng = -(m[iid] * m[jid] / norm(dr)) * consts::G;
                        Scalar kinetic_eng_i = 0.5 * dot(v[iid], v[iid]) * m[iid];
                        Scalar kinetic_eng_j = 0.5 * dot(v[jid], v[jid]) * m[jid];
                        auto kinetic_eng = kinetic_eng_i + kinetic_eng_j;
                        auto total_energy = kinetic_eng + potential_eng;
                        energy_file << i << "," << j << "," << total_energy << "\n";
                         
                    }
                }

            };
    */
    //Solver: t_start, particle1, particle2, ...
    Solver solver{0, p1, p2};
    Solver::RunArgs args;
    auto stop_time = 800_year;
    args.add_stop_condition(stop_time);

    auto twriter = TimeSlice(DefaultWriter("testdynfriction.csv"), 0.0, stop_time, 500);
    //auto eWriter = TimeSlice(energy_log, 0.0, 100.0, 50);
    args.add_operation(twriter);
    //args.add_operation(eWriter);
    

    solver.run(args);

    print(std::cout, "Simulation Complete!\n");

    return 0;
}
// g++ -std=c++17 -O3 -pthread simulations/dragforces/testforce-dyn.cpp -o simulations/dragforces/testforce-dyn
// simulations/dragforces/testforce-dyn