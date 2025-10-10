
#include "../../src/spaceHub.hpp"
#include "../../src/rand-generator.hpp"
#include "../../src/taskflow/taskflow.hpp"
#include <typeinfo>
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace orbit; // save writing orbit::
using namespace force;
/*----------------------------------------------------------------------------------------------------------------*/
using f = Interactions<NewtonianGrav, PN2p5>;
using Solver = methods::AR_Chain_Plus<f, particles::SizeParticles>;

using Particle = Solver::Particle;

using Scalar = Solver::Scalar;

void job(std::vector<std::array<Scalar, 3>> &combinations, size_t n_start, size_t n_stop)
{

    std::fstream ptc_res_file("simulation_results/mt_dump/fig2_ptc_" + std::to_string(n_start) + ".txt", std::ios::out);
    std::fstream output_summary_file("simulation_results/mt_dump/fig2_ics_" + std::to_string(n_start) + ".txt", std::ios::out);

    if (n_start == 0)
    {
        ptc_res_file << "time,id,mass,px,py,pz,vx,vy,vz" << '\n';
        output_summary_file << "b,azi,incl,ci,cj,exitcon" << '\n';
    }
    // IC statics
    Scalar Mijk = 10_Ms;
    Scalar Rsijk = 4.2450051e-5_Rs;
    Scalar v_inf = 16.314_kms; // Yes

    Scalar ab = 1e-2_AU;

    Scalar b_i;
    Scalar azi_j;
    Scalar incl;

    print(std::cout << "Begin job i=" << n_start << '\n');

    for (size_t i = n_start; i < n_stop; i++)
    {
        auto bfi = combinations[i];
        b_i = bfi[0];
        azi_j = bfi[1];
        incl = bfi[2];

        auto r_start = sqrt(pow((20 * ab), 2) + pow(b_i, 2));

        char exit_condition = 't';
        // print(std::cout << "Begin sim (bfi): (" << b_i << "," << azi_j << "," << incl << ")" << '\n');

        
            Particle p1{Mijk, Rsijk}; // R = 20GM_sun/c^2
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
            args.rtol = 1e-12; // Try to minimize number of failed sims
            int ci = -1;
            int cj = -1;

            /*------------------------------------STOP CONDITIONS--------------------------------------------------*/
            // Time out stop condition
            
            Scalar ttp = time_to_periapsis(orbit::group(p1, p2), p3);
            Scalar t_end = 60 * ttp;

            // to do: args.rtol <-- vary this
            

            // Collision stop condition
            auto collision_detect = [&ci, &cj, &exit_condition](auto &ptc, auto h)
            {
                size_t particle_num = ptc.number();
                for (size_t i = 0; i < particle_num; ++i)
                {
                    for (size_t j = i + 1; j < particle_num; ++j)
                    {

                        if (distance(ptc.pos(i), ptc.pos(j)) < ptc.radius(i) + ptc.radius(i))
                        {
                            ci = i;
                            cj = j;
                            print(std::cout << "collision" << "\n");
                            exit_condition = 'c';
                            return true;
                        }
                    }
                }
                return false;
            };

           auto particle_ejected = [&ttp](auto &ptc, auto h)
            {
                //print(std::cout << ptc.time() << "<" << ttp <<'\n');
                if (ptc.time() >  4 * ttp)
                {
                    print(std::cout << ptc.time() << ">" << ttp <<'\n');
                    size_t particle_num = ptc.number();
                    for (size_t i = 0; i < particle_num; ++i)
                    {
                        for (size_t j = i + 1; j < particle_num; ++j)
                        {

                            auto total_mass = 30_Ms;
                            auto u = total_mass * 1;
                            auto dv = ptc.vel(i) - ptc.vel(j);
                            auto dr = ptc.pos(i) - ptc.pos(j);
                           
                            auto escalar = orbit::calc_eccentricity(u, dr, dv);
                            print(std::cout << "e: " << escalar << '\n');
                            if (escalar >= 1)
                            {
                                print(std::cout << "ejection: " << escalar << '\n');
                                return true;
                            }
                        }
                    }
                }
                return false;
            };
        

            args.add_stop_condition(t_end);

            args.add_stop_condition(particle_ejected);

            args.add_stop_condition(collision_detect);

            // Assume merger stop condition

            /*------------------------------------EXIT OPERATIONS--------------------------------------------------*/
            /*NEGATIVE TIMESTAMP CODES:
            -1: Exception caught while running
            -2: Printout of b, f, i, (exit condition?)
            */

            args.add_stop_point_operation([&ptc_res_file, &output_summary_file, &b_i, &azi_j, &incl, &ci, &cj, &exit_condition](auto &ptc, auto h)
                                          {
                                    // print the end state of the system into file
                                    // time,pxyz,vxyz for each particle in the system
                                    output_summary_file << '\n'
                                                        << b_i << "," << azi_j << "," << incl << "," << ci << "," << cj << ',' << exit_condition << "\n";
                                    ptc_res_file << ptc << '\n'; 
                                });
        try
        {
            solver.run(args);
        }
        catch (...)
        {
            
            print(std::cout << b_i << "," << azi_j << "," << incl << '\n');

            print(std::cout << "failed" << "\n");
        }
        print(std::cout, "job ", i, " finished\n");
    }
}

int main(int argc, char **argv)
{

    Scalar b_min = 4_AU; // [-0.04, 0.04] Corresponds to [-4, 4] in the rescaled b parameter from paper for a0=10^-4
    Scalar b_max = -4_AU;
    Scalar azi_max = 2 * consts::pi;
    Scalar azi_min = 0 * consts::pi;

    Scalar incl_range[2] = {0, 1 * consts::pi};
    int n_incl = 1;

    // pre-fill column headers

    size_t n = 10; // b-f grid size is n X n

    std::vector<std::array<Scalar, 3>> combinations;

    for (size_t i = 0; i < n; ++i)
    {
        Scalar b_i = b_min + i * ((b_max - b_min) / (n - 1));

        for (size_t j = 0; j < n; ++j)
        {

            Scalar azi_j = 0 + j * ((azi_max - azi_min) / (n - 1));

            for (size_t k = 0; k < n_incl; ++k)
            {
                Scalar incl = 0; // Default

                if (n_incl > 1)
                {
                    Scalar incl = 0 + k * ((incl_range[1] - incl_range[0]) / (n_incl - 1));
                }

                std::array<Scalar, 3> ijk_arr = {b_i, azi_j, incl};
                combinations.push_back(ijk_arr);
            }
        }
    }

    tf::Executor executor;
    tools::Timer timer; // Leaving in for fun
    timer.start();

    size_t jobs_per_file = 500;
    size_t num_combos = combinations.size();
    for (size_t i = 0; i < num_combos; i += jobs_per_file)
    {

        std::vector<std::array<Scalar, 3>> combos = combinations;
        executor.silent_async(job, combos, i, (i + jobs_per_file));
    }
    executor.wait_for_all();

    auto output_size = n * n * n_incl;
    std::cout << "Complete! Output size: " << n << "^2 * " << n_incl << " = " << output_size << std::endl;
    // g++ -std=c++17 -O3 -pthread simulations/topology/binary-single-3debug.cpp -o simulations/binary-single-3debug

    // simulations/binary-single-3demo
    return 0;
}
