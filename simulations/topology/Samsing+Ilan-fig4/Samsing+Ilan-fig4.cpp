
#include "../../../SpaceHub/src/spaceHub.hpp"
#include "../../../SpaceHub/src/rand-generator.hpp"
#include "../../../SpaceHub/src/taskflow/taskflow.hpp"
#include <typeinfo>
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace orbit; // save writing orbit::
using namespace force;
/*----------------------------------------------------------------------------------------------------------------*/
using f = Interactions<NewtonianGrav, PN2p5>;
using Solver = methods::DefaultMethod<f, particles::SizeParticles>;

using Particle = Solver::Particle;

using Scalar = Solver::Scalar;

void job(std::vector<std::array<Scalar, 3>> &combinations, size_t n_start, size_t n_stop)
{

    std::fstream ptc_res_file("simulations/topology/Samsing+Ilan-fig4/results/fig4_ptc_" + std::to_string(n_start) + ".txt", std::ios::out);
    std::fstream errlog_file("simulations/topology/Samsing+Ilan-fig4/results/logs/errlog_" + std::to_string(n_start) + ".txt", std::ios::out);
    

    if (n_start == 0)
    {
        ptc_res_file << "time,id,mass,px,py,pz,vx,vy,vz" << '\n';
        
    }
    // IC statics
    Scalar Mijk = 10_Ms;
    Scalar Rsijk = 4.2450051e-5_Rs;
    Scalar v_inf = 115.35562_kms; //Calculated in notebook

    Scalar ab = 1e-4_AU;

    Scalar b_i;
    Scalar azi_j;
    Scalar incl;

    //print(std::cout << "Begin job i=" << n_start << '\n');

    for (size_t i = n_start; i < n_stop; i++)
    {
        auto bfi = combinations[i];
        b_i = bfi[0];
        azi_j = bfi[1];
        incl = bfi[2];
        

        auto r_start = (20 * ab);

       int exit_condition = 0;
        /*
        -1: Error
        0: time stop
        1: collision
        2: GW inspiral (TODO)
        */

        try
        {
            Particle p1{Mijk, Rsijk}; // R = 20GM_sun/c^2
            Particle p2{Mijk, Rsijk};
            Particle p3{Mijk, Rsijk};
            /*--------------------------------------------------New-----------------------------------------------------------*/

            // Create binary elliptic orbit with M1=M2=1Ms, semimajor axis=5.0, e=0. incl, Omega, omega, nu
            /*In 2D, long. of asc. node, argument of periapsis, or true anomaly could be used to define the azimuth angle between
            the incident and binary orbit, however since I am not sure how spacehub defines longitude of asc node when incl=0,
            we will set Omega and nu to zero and vary the arg of periapsis
            --> Or this could be set in the incident orbit? Probably easier here*/
            Scalar incl_offset = 0_deg;
            Scalar Torb;
            if (b_i < 0){
                b_i = -b_i;
                incl_offset = 180_deg;
                
            }
            
            
            auto binary_orb = Elliptic(p1.mass, p2.mass, ab, 0.0, incl_offset, 0.0, azi_j, 0.0);
                move_particles(binary_orb, p2);
                Torb = period(binary_orb);
             // move p2 to the corresponding position/velocity of the orbit around origin(p1)

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
            Scalar t_end = 500 * Torb;
            
            Scalar E0  = 0;
            
            auto calc_e0 = [&E0](auto& ptc, auto step_size) {

                    int iid = 0;
                    int jid = 1;
                    auto const &m = ptc.mass();
                    auto const &p = ptc.pos();
                    auto const &v = ptc.vel();
                    decltype(p[iid]) dr = p[iid] - p[jid];
                    Scalar potential_eng = -(m[iid] * m[jid] / norm(dr)) * consts::G;
                    Scalar kinetic_eng_i = 0.5 * dot(v[iid], v[iid]) * m[iid];
                    Scalar kinetic_eng_j = 0.5 * dot(v[jid], v[jid]) * m[jid];
                    auto kinetic_eng = kinetic_eng_i + kinetic_eng_j;
                    auto total_energy = kinetic_eng + potential_eng;
                    E0 = total_energy;
                    return 0;
            };

            args.add_start_point_operation(calc_e0);


            auto gw_energy_stop = [&E0, &ci, &cj, &exit_condition](auto& ptc, auto step_size) {
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
                        if (total_energy < 100 * E0)
                        {
                            exit_condition = 2;
                            ci = i;
                            cj = j;
                            print(std::cout, "GW end");
                            return true;
                        }
                        
                    }
                }
                return false;

            };

            args.add_stop_condition(gw_energy_stop);
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
                            exit_condition = 1;
                            return true;
                        }
                    }
                }
                return false;
            };

       
        //print(std::cout << E0 << "\n");
        
        

        
        
        


            

            args.add_stop_condition(t_end);
            

            //args.add_stop_condition(particle_ejected);

            args.add_stop_condition(collision_detect);

            // Assume merger stop condition

            /*------------------------------------EXIT OPERATIONS--------------------------------------------------*/
            /*NEGATIVE TIMESTAMP CODES:
            -1: Exception caught while running
            -2: Printout of b, f, i, (exit condition?)
            */

            args.add_stop_point_operation([&ptc_res_file, &errlog_file, &i, &bfi, &azi_j, &incl, &ci, &cj, &exit_condition, &t_end](auto &ptc, auto h)
                                          {
                                    // print the end state of the system into file
                                    // time,pxyz,vxyz for each particle in the system
                                    print(std::cout << "Simulation #" << i << "complete \n");
                                    ptc_res_file << ptc << ptc.time() << ",-2," << bfi[0] << "," << azi_j << "," << incl << "," << exit_condition << "," << ci << "," << cj << ",0\n" << '\n'; });
            solver.run(args);
            
        }
        catch (...)
        {   
            print(std::cout << b_i << "," << azi_j << "," << incl);

            print(std::cout << " failed" << "\n");
        }
        //print(std::cout, "job ", i, " finished\n");
    }
}

int main(int argc, char **argv)
{

    Scalar b_min = -0.04_AU; // [-0.04, 0.04] Corresponds to [-4, 4] in the rescaled b parameter from paper for a0=10^-4
    Scalar b_max = 0.04_AU;
    Scalar azi_max = 2 * consts::pi;
    Scalar azi_min = 0 * consts::pi;

    Scalar incl_range[2] = {0, 1 * consts::pi};
    int n_incl = 1;

    // pre-fill column headers

    size_t n = 200; // b-f grid size is n X n

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

    size_t jobs_per_file = n;
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


