
#include "../../src/spaceHub.hpp"
#include "../../src/rand-generator.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace orbit; // save writing orbit::
using namespace force;
/*----------------------------------------------------------------------------------------------------------------*/
using f = Interactions<NewtonianGrav, PN1, PN2, PN2p5>;
using Solver = methods::AR_Chain_Plus<f, particles::SizeParticles>;

using Particle = Solver::Particle;

using Scalar = Solver::Scalar;

int main(int argc, char **argv)
{

    //IC range set-up
    Scalar Mijk = 10_Ms;
    Scalar Rsijk = 4.2450051e-5_Rs;
    Scalar v_inf = 16.314_kms; //Yes
    Scalar r_start = 20_AU;
    Scalar ab = 1e-2_AU;

    Scalar b_min = -4_AU;
    Scalar b_max = 4_AU;
    Scalar azi_max = 2 * consts::pi;
    Scalar azi_min = 0 * consts::pi;

    Scalar incl_range[2] = {0, 1 * consts::pi};
    int n_incl = 1;
    // std::fstream incident_orb_res_file("simulation_results/2+1_incident_orbres-2d.txt", std::ios::out);
    // std::fstream inner_orb_res_file("simulation_results/2+1_inner_orbres-2d.txt", std::ios::out);
    std::fstream ptc_res_file("simulation_results/topology/SI-fig2-ptcres.txt", std::ios::out);
    std::fstream ordered_inputs_file("simulation_results/topology/SI-fig2-ICs.txt", std::ios::out);
    // pre-fill column headers
    ptc_res_file << "time,id,mass,px,py,pz,vx,vy,vz" << '\n';
    // incident_orb_res_file << "m1,m2,slr,e,i,Omega,omega,nu" << '\n';
    // inner_orb_res_file << "m1,m2,slr,e,i,Omega,omega,nu" << '\n';
    ordered_inputs_file << "b,azi,incl,ci,cj" << '\n';


    size_t n = 10; //b-f grid size is n X n

    for (size_t i = 0; i < n; ++i)
    {
        Scalar b_i = b_min + i * ( (b_max - b_min) / (n - 1) );

        for (size_t j = 0; j < n; ++j)
        {

            auto azi_j = 0 + j * ( (azi_max - azi_min) / (n - 1) ); 

            for (size_t k = 0; k < n_incl; ++k)
            {
                auto incl = 0; //Default

                if (n_incl > 1) {
                    auto incl = 0 + k * ( (incl_range[1] - incl_range[0]) / (n_incl - 1) );
                }
                print(std::cout << "Begin sim (ijk): (" << i << "," << j << "," << k << ")" << '\n');
                
                try {
                Particle p1{Mijk, Rsijk}; //R = 20GM_sun/c^2
                Particle p2{Mijk, Rsijk};
                Particle p3{Mijk, Rsijk};
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

                /*--------------------------------------------------Etc-----------------------------------------------------------*/
                 args.rtol = 1e-12; //Try to minimize number of failed sims
                 int ci = -1;
                 int cj = -1;

                /*------------------------------------STOP CONDITIONS--------------------------------------------------*/
                // Time out stop condition
                Scalar t_end = 6 * time_to_periapsis(orbit::group(p1, p2), p3);
                
                // to do: args.rtol <-- vary this
                args.add_stop_condition(t_end);
                    
                //Collision stop condition
                auto collision_detect = [&ci, &cj](auto& ptc, auto h) {
                    size_t particle_num = ptc.number();
                    for (size_t i = 0; i < particle_num; ++i) {
                        for (size_t j = i + 1; j < particle_num; ++j) {
                            if (distance(ptc.pos(i), ptc.pos(j)) < ptc.radius(i) + ptc.radius(i)) {
                                ci = i;
                                cj = j;
                                return true;
                            }
                        }
                    }
                    return false;
                };

                args.add_stop_condition(collision_detect);


                //Assume merger stop condition

                /*------------------------------------EXIT OPERATIONS--------------------------------------------------*/
                args.add_stop_point_operation([&ptc_res_file, &ordered_inputs_file, &b_i, &azi_j, &incl, &ci, &cj](auto &ptc, auto h)
                                {
                                    // print the end state of the system into file
                                    // time,pxyz,vxyz for each particle in the system
                                    ordered_inputs_file << '\n'
                                                        << b_i << "," << azi_j << "," << incl << "," << ci << "," << cj << "\n";
                                    ptc_res_file << ptc << '\n';
                                });
                solver.run(args);
                }
                catch(...){
                    //Flag these in case I want to plot/debug them in python
                    ptc_res_file << "-1," << b_i << "," << azi_j << "," << incl << ",0,0,0,0,0\n-2,0,0,0,0,0,0,0,0\n-3,0,0,0,0,0,0,0,0\n" << '\n';
                    print(std::cout << b_i << "," << azi_j << "," << incl << '\n');
                }
            }
        }
    }


    auto output_size = n*n*n;
      std::cout << "Complete! Output size: " << n << "^3 = " << output_size << std::endl;
    // g++ -std=c++17 -O3 -pthread simulations/topology/binary-single-3debug.cpp -o simulations/binary-single-3debug

    // simulations/binary-single-3demo
    return 0;
}
