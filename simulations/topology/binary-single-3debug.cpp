
#include "../../src/spaceHub.hpp"
#include "../../src/rand-generator.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace orbit; // save writing orbit::
/*----------------------------------------------------------------------------------------------------------------*/
using Solver = methods::AR_Chain_Plus<>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char **argv)
{
    Scalar v_inf = 16.314_kms;
    Scalar r_start = 100_AU;
    Scalar ab = 10_AU;
    // std::fstream incident_orb_res_file("simulation_results/2+1_incident_orbres-2d.txt", std::ios::out);
    // std::fstream inner_orb_res_file("simulation_results/2+1_inner_orbres-2d.txt", std::ios::out);
    std::fstream ptc_res_file("simulation_results/topology/3d2+1_debug.txt", std::ios::out);
    std::fstream ordered_inputs_file("simulation_results/topology/3d2+1_ICs_debug.txt", std::ios::out);
    // pre-fill column headers
    ptc_res_file << "time,id,mass,px,py,pz,vx,vy,vz" << '\n';
    // incident_orb_res_file << "m1,m2,slr,e,i,Omega,omega,nu" << '\n';
    // inner_orb_res_file << "m1,m2,slr,e,i,Omega,omega,nu" << '\n';
    ordered_inputs_file << "b,azi,incl" << '\n';

    Scalar b_min = 1_AU;
    Scalar b_max = 30_AU;

    auto incl_max = 1 * consts::pi;
    size_t n = 50; // sims per parameter will be n^3

    for (size_t i = 0; i < n; ++i)
    {
        Scalar b_i = b_min + i * ((b_max - b_min) / (n - 1));

        for (size_t j = 0; j < n; ++j)
        {

            auto azi_j = 0 + j * ((2 * consts::pi - 0) / (n - 1)); // Saving iterations since 2d sim shows symmetry

            for (size_t k = 0; k < n; ++k)
            {
                

                auto incl = 0 + k * ((incl_max - 0) / (n - 1));
                
                
                try {
                Particle p1{1_Ms}; // Create particles of equal mass at rest at origin
                Particle p2{1_Ms};
                Particle p3{1_Ms};
                /*--------------------------------------------------New-----------------------------------------------------------*/

                // Create binary elliptic orbit with M1=M2=1Ms, semimajor axis=5.0, e=0. incl, Omega, omega, nu
                /*In 2D, long. of asc. node, argument of periapsis, or true anomaly could be used to define the azimuth angle between
                the incident and binary orbit, however since I am not sure how spacehub defines longitude of asc node when incl=0,
                we will set Omega and nu to zero and vary the arg of periapsis
                --> Or this could be set in the incident orbit? Probably easier here*/
                auto binary_orb = Elliptic(p1.mass, p2.mass, ab, 0.0, 0.0, 0.0, azi_j, 0.0);

                move_particles(binary_orb, p2); // move p2 to the corresponding position/velocity of the orbit around origin(p1)

                move_to_COM_frame(p1, p2); // sets origin to the center of mass between p1/p2?

                // sig: incident_orbit(target_mass, incident_mass, v_inf, b_max, r_rel between target and incident)

                auto incident_orb = orbit::Hyperbolic(M_tot(p1, p2), p3.mass, v_inf, b_i, incl, 0.0, 0.0, r_start, orbit::Hyper::in);

                // Put p3 on the incident orbit trajectory
                move_particles(incident_orb, p3);

                // Move the origin to the com of all the particles? why?
                move_to_COM_frame(p1, p2, p3);
                /*----------------------------------------------------------------------------------------------------------------*/

                Solver solver{0, p1, p2, p3};

                Solver::RunArgs args;

                /*--------------------------------------------------New-----------------------------------------------------------*/
                // orbit::group(p1, p2) is the target object, p3 in the incident object
                
                Scalar t_end = 6 * time_to_periapsis(orbit::group(p1, p2), p3);
                
                // to do: args.rtol <-- vary this
                args.add_stop_condition(t_end);
                args.rtol = 1e-12;
                /*----------------------------------------------------------------------------------------------------------------*/
                // [&] capture all variables in lambda by reference
                args.add_stop_point_operation([&ptc_res_file, &ordered_inputs_file, &b_i, &azi_j, &incl](auto &ptc, auto h)
                                              {
                                                  // print the end state of the system into file
                                                  // time,pxyz,vxyz for each particle in the system
                                                  ordered_inputs_file << '\n'
                                                                      << b_i << "," << azi_j << "," << incl << "\n";
                                                  ptc_res_file << ptc << '\n';
                                              });
                    
                
                solver.run(args);
                }
                catch(...){
                    print(std::cout << b_i << "," << azi_j << "," << incl << '\n');
                }
            }
        }
    }


    auto output_size = n*n*n;
    std::cout << "Complete! Output size: " << n << "^3" << std::endl;
    // g++ -std=c++17 -O3 -pthread simulations/topology/binary-single-3debug.cpp -o simulations/binary-single-3debug

    // simulations/binary-single-3demo
    return 0;
}
