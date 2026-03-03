
#include "../../SpaceHub/src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace force;  // save writting force::
// add whatever external term you like, just keep the first one to be the internal Newtonian gravity. check available force at xxx
using f = Interactions<NewtonianGrav, PN2p5>;

using Solver = methods::BS<f>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {
    Particle p1{1.4_Ms};
    Particle p2{1.4_Ms};

    Scalar a = 0.016357768_AU;
    Scalar ecc = 0;
    Scalar incl = 0_deg;
    Scalar l_ascnode = 0_deg;
    Scalar arg_of_periapsis = 0_deg;
    Scalar true_anom = 0_deg;

    auto inner_orb = orbit::Elliptic(p1.mass, p2.mass, a, ecc, incl, l_ascnode, arg_of_periapsis, true_anom);

    orbit::move_particles(inner_orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    //Solver: t_start, particle1, particle2, ...
    Solver solver{0, p1, p2};
    Solver::RunArgs args;

    Scalar E0 = 0;

 

    args.add_operation( StepSlice([&](auto &ptc, auto step_size) {
        
        //------------
        int p1id = 0; int p2id = 1; 


        auto const &m = ptc.mass();
        auto const &p = ptc.pos();
        auto const &v = ptc.vel();

        
        decltype(p[p1id]) dr = p[p1id] - p[p2id];
        Scalar potential_eng = -(m[p1id] * m[p2id] / norm(dr)) * consts::G;
        Scalar kinetic_eng_i = 0.5 * dot(v[p1id], v[p1id]) * m[p1id];
        Scalar kinetic_eng_j = 0.5 * dot(v[p2id], v[p2id]) * m[p2id];
        auto kinetic_eng = kinetic_eng_i + kinetic_eng_j;
        
        //print(std::cout << "TEST T: " << kinetic_eng << "U: " << potential_eng << "\n");

        auto total_energy = kinetic_eng + potential_eng;
    

       



        //------------
        
        }, 10)   );

    args.add_stop_condition(1.5e9_year);

    auto collision_detect = [&](auto &ptc, auto h)
            {
                size_t particle_num = ptc.number();
                for (size_t i = 0; i < particle_num; ++i)
                {
                    for (size_t j = i + 1; j < particle_num; ++j)
                    {

                        if (distance(ptc.pos(i), ptc.pos(j)) < 1e-5_AU)
                        {

                            return true;
                        }
                    }
                }
                return false;
            };
    args.add_stop_condition(collision_detect);

   args.add_operation( StepSlice([&](auto &ptc, auto step_size) {
        
        //------------
        print(std::cout, "1000 Steps Completed\n");
        std::cout << ptc;

        //------------
        
        }, 1000)   );
    //args.add_operation(StepSlice(print_total_energy, 10));
    args.add_operation(StepSlice(DefaultWriter("Hulse-Taylor.dat"), 1000));

    solver.run(args);

    print(std::cout, "Simulation Complete!\n");

    return 0;
}

