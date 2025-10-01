
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;
/*--------------------------------------------------New-----------------------------------------------------------*/
using Scalar = Solver::Scalar;  // This is not explained in detail as far as I remember

int main(int argc, char** argv) {

    Scalar v_inf = 10_kms;
    Scalar b_max = 10_AU;
    Scalar r_start = 100_AU;  // drop the incident object 100 AU away from the scattered object. The trajectory from
                              // +inf to r_start will be calculated analytically

    std::fstream file("simulation_results/scattering_demo_finalstate.txt", std::ios::out);

    Particle p1{1_Ms};
    Particle p2{1_Ms};

    /*--------------------------------------------------New-----------------------------------------------------------*/
    // create the incident orbit (total mass of scattered object, total mass of incident object, velocity at
    // infinity, max impact parameter, relative distance at t = 0)

    //b and w are randomized automatically here by default
    //auto orb = scattering::incident_orbit(p1.mass, p2.mass, v_inf, b_max, r_start);

    //Manually set up incident orbit to avoid randomized parameters
    auto orb = orbit::Hyperbolic(p1.mass, p2.mass, v_inf, b_max, 0.0, 0.0, 0.0, r_start, orbit::Hyper::in);
    
    orbit::move_particles(orb, p2);

    orbit::move_to_COM_frame(p1, p2);
    /*----------------------------------------------------------------------------------------------------------------*/

    Solver solver{0, p1, p2};

    Solver::RunArgs args;

    /*--------------------------------------------------New-----------------------------------------------------------*/
    Scalar t_end = 3 * orbit::time_to_periapsis(p1, p2);  // calculate the time from start position to periapsis

    args.add_stop_condition(t_end);
    /*----------------------------------------------------------------------------------------------------------------*/

    args.add_stop_point_operation([&file, &orb](auto& ptc, auto h) {
        // print the end state of the system into file
        print(file, "orb:", orb, "\nptc:", ptc, '\n');
    });

    auto t_writer = TimeSlice(DefaultWriter("simulation_results/Scattering_demo_1+1.txt"), 0.0, t_end, 1000);
    args.add_operation(t_writer);

    solver.run(args);


    print(std::cout, "simulation with scattering tools complete!\n");

    return 0;
    //rm simulations/Scattering_demo_1+1

    //rm simulation_results/Scattering_demo_1+1.txt

    //rm simulation_results/scattering_demo_finalstate.txt

    //g++ -std=c++17 -O3 -pthread simulations/Scattering_demo_1+1.cpp -o simulations/Scattering_demo_1+1

    //simulations/Scattering_demo_1+1
}