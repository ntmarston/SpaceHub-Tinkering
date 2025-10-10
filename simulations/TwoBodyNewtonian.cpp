
#include "../src/spaceHub.hpp"
#include "../src/orbits/orbits.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace force;  
using f = Interactions<NewtonianGrav, PN2p5>;

using Solver = methods::DefaultMethod<f>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;

int main(int argc, char** argv) {

    Particle p1{1000_Ms};
    Particle p2{1_Ms};


    //Elliptic parameters: M1, M2, initial: semimajor axis, eccentricity, inclination, longitude of ascending node, argument of periapsis, true anomaly
    auto orb = orbit::Elliptic(p1.mass, p2.mass, 1e-2_AU, 0.8, 45_deg, 270_deg, 60_deg, 0_deg);
    //a=5AU, e = 0.1, i=45deg, \Omega=270deg, \omega = 90deg nu = 0deg
    orbit::move_particles(orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    //Solver: t_start, particle1, particle2, ...
    Solver solver{0, p1, p2};
    
    Solver::RunArgs args;

    args.add_stop_condition(1000_year);

    auto particle_ejected = [&orb](auto &ptc, auto h)
                {
                    //print(std::cout << ptc.time() << "<" << ttp <<'\n');
                    
                        
                        size_t particle_num = ptc.number();
                        for (size_t i = 0; i < particle_num; ++i)
                        {
                            for (size_t j = i + 1; j < particle_num; ++j)
                            {
                                auto u = 1 * (1001_Ms);
                                //print(std::cout << ptc.pos(i) << " - " << ptc.pos(j) << "=\n");
                                auto dr = ptc.pos(i) - ptc.pos(j);
                                //print(std::cout << dr << "\n");
                                auto dv = ptc.vel(i) - ptc.vel(j);
                                auto e = orbit::calc_eccentricity(u, dr, dv);
                                print(std::cout << e << "\n");
                                //print(std::cout << orb.nu << "\n");
                            }
                        }
                    
                    return false;
                };
    args.add_operation(DefaultWriter("simulation_results/Result.txt"));
    args.add_stop_condition(particle_ejected);

    solver.run(args);

    print(std::cout, "Simulation complete");

    return 0;
}
