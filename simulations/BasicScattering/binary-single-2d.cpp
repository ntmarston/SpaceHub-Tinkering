
#include "../SpaceHub/src/spaceHub.hpp"
#include "../SpaceHub/src/rand-generator.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace orbit;  // save writing orbit::
/*----------------------------------------------------------------------------------------------------------------*/
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {
    Scalar v_inf = 13.5_kms;
    Scalar r_start = 100_AU;
    Scalar ab = 15_AU;
    std::fstream incident_orb_res_file("simulations/BasicScattering/results/2+1_incident_orbres-2d.txt", std::ios::out);
    std::fstream inner_orb_res_file("simulations/BasicScattering/results/2+1_inner_orbres-2d.txt", std::ios::out);
    std::fstream ptc_res_file("simulations/BasicScattering/results/2+1_ptcres-2d.txt", std::ios::out);
    // pre-fill column headers
    ptc_res_file << "time,id,mass,px,py,pz,vx,vy,vz" << '\n';
    incident_orb_res_file << "m1,m2,slr,e,i,Omega,omega,nu" << '\n';
    inner_orb_res_file << "m1,m2,slr,e,i,Omega,omega,nu" << '\n';
    Scalar b_min = 0.1_AU;
    Scalar b_max = 30_AU;
    
    auto incl = 0 * consts::pi;
    size_t n = 50; //sims per parameter will be n^2

    for(size_t i=0; i < n; ++i)
    {
        Scalar b_i = b_min + i*((b_max - b_min) / (n-1));

        for(size_t j=0; j < n; ++j){
            Particle p1{1_Ms}; //Create particles of equal mass at rest at origin
            Particle p2{1_Ms};
            Particle p3{1_Ms};
            
            auto azi_j = 0 + j*((2 * consts::pi - 0) / (n-1));
            /*--------------------------------------------------New-----------------------------------------------------------*/

            //Create binary elliptic orbit with M1=M2=1Ms, semimajor axis=5.0, e=0. incl, Omega, omega, nu
            /*In 2D, long. of asc. node, argument of periapsis, or true anomaly could be used to define the azimuth angle between 
             the incident and binary orbit, however since I am not sure how spacehub defines longitude of asc node when incl=0,
            we will set Omega and nu to zero and vary the arg of periapsis 
            --> Or this could be set in the incident orbit? Probably easier here*/
            auto binary_orb = Elliptic(p1.mass, p2.mass, ab, 0.0, 0.0, 0.0, azi_j, 0.0);

            move_particles(binary_orb, p2); //move p2 to the corresponding position/velocity of the orbit around origin(p1)

            move_to_COM_frame(p1, p2); //sets origin to the center of mass between p1/p2?

        
            //sig: incident_orbit(target_mass, incident_mass, v_inf, b_max, r_rel between target and incident)
            
            auto incident_orb = orbit::Hyperbolic(M_tot(p1, p2), p3.mass, v_inf, b_i, 0.0, 0.0, 0.0, r_start, orbit::Hyper::in);

            //Put p3 on the incident orbit trajectory
            move_particles(incident_orb, p3);

            //Move the origin to the com of all the particles? why? 
            // Would it not make more sense to stay in the reference frame of the binary?
            move_to_COM_frame(p1, p2, p3);
            /*----------------------------------------------------------------------------------------------------------------*/

            Solver solver{0, p1, p2, p3};

            Solver::RunArgs args;

            /*--------------------------------------------------New-----------------------------------------------------------*/
            // orbit::group(p1, p2) is the target object, p3 in the incident object
            Scalar t_end = 6 * time_to_periapsis(orbit::group(p1, p2), p3);

            args.add_stop_condition(t_end);
            /*----------------------------------------------------------------------------------------------------------------*/
            // [&] capture all variables in lambda by reference
            args.add_stop_point_operation([&ptc_res_file, &b_i, &azi_j, &incl](auto& ptc, auto h) {
                // print the end state of the system into file
                // time,pxyz,vxyz for each particle in the system
                ptc_res_file << '\n' << b_i << "," << azi_j << "," << incl << "\n";
                ptc_res_file << ptc << '\n';
          
            });

            args.add_stop_point_operation([&incident_orb_res_file, &incident_orb, &i, &b_i, &azi_j](auto& ptc, auto h) {
                // print the end state of the system into file
                // m1, m2, a(1-e^2), e, i, longOfAscNode, arg of periapsis, true anomaly
                print(incident_orb_res_file, i,',',b_i,',',azi_j,',', incident_orb, '\n');
            });

            args.add_stop_point_operation([&inner_orb_res_file, &binary_orb, &i, &b_i, &azi_j](auto& ptc, auto h) {
                // print the end state of the system into file
                // m1, m2, a(1-e^2), e, i, longOfAscNode, arg of periapsis, true anomaly
                print(inner_orb_res_file, i,',',b_i,',',azi_j,',', binary_orb, '\n');
            });
            
            if(i==25 && j==25)
            {
            auto writer = StepSlice(DefaultWriter("simulations/BasicScattering/results/2+1_sample-2d.txt"), 5);
            args.add_operation(writer);
            }

            solver.run(args);
             
        }
    }

    std::cout << "simulation of single binary scattering complete!\n";
    //g++ -std=c++17 -O3 -pthread simulations/binary-single-2d.cpp -o simulations/binary-single-2d

    //simulations/binary-single-2d
    return 0;
}
