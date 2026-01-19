
#include "../../src/spaceHub.hpp"
#include "../../src/rand-generator.hpp"
#include "../../src/taskflow/taskflow.hpp"
#include <typeinfo>
#include <mutex>
#include <fstream>
#include <random>
#include <future>
#include <chrono>
#include <ctime>
#include <iomanip>
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
using namespace orbit;

// Global mutex for thread-safe logging
std::mutex log_mutex;

// External timeout in seconds
constexpr int EXTERNAL_TIMEOUT_SECONDS = 120;

// Container for hung futures to prevent blocking on destruction
std::vector<std::future<void>> hung_futures;
std::mutex hung_futures_mutex;
/*--------------------------------------------------New-----------------------------------------------------------*/
// using the Newtonian gravity + first order Post-Newtonian correction
using f = Interactions<NewtonianGrav, StaticGasField>;

using Solver = methods::BS<f, particles::SizeParticles>;
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

void job(std::vector<std::array<Scalar, 4>> &combinations, size_t n_start, size_t n_stop) {
    for (size_t i = n_start; i < n_stop; i++)
    {
        auto inputs = combinations[i];
        
        
        auto mratio = inputs[0]; //     m2/m1 = q
        Scalar m1 = 1_Ms;
        Scalar m2 = mratio * m1;
     

        Particle p1{m1, 1_Rs};
        Particle p2{m2, 1_Rs};

        Scalar sma = inputs[1];
        auto ecc = inputs[2];
        auto inclination = inputs[3];
        auto longitude_of_ascending_node = 0_deg;
        auto argument_of_periapsis = 0_deg; 
        auto true_anomaly = 0_deg;

        auto orb = orbit::Elliptic(p1.mass, p2.mass, sma, ecc, inclination, longitude_of_ascending_node, argument_of_periapsis, true_anomaly);



        orbit::move_particles(orb, p2);

        orbit::move_to_COM_frame(p1, p2);

        Solver solver{0, p1, p2};

        Solver::RunArgs args;

        auto collision_detect = [](auto &ptc, auto h)
                {
                    //print(std::cout << "time:" << ptc.time() << "step:" << h << "\n");
                    size_t particle_num = ptc.number();
                    for (size_t i = 0; i < particle_num; ++i)
                    {
                        for (size_t j = i + 1; j < particle_num; ++j)
                        {

                            if (distance(ptc.pos(i), ptc.pos(j)) < ptc.radius(i) + ptc.radius(i))
                            {
                                return true;
                            }
                        }
                    }
                    return false;
                };

        args.add_stop_condition(collision_detect);


        auto Torb = period(orb);
        Scalar orbital_endtime = 500 * Torb;
        args.add_stop_condition(orbital_endtime);
        args.add_stop_condition(1000_year);

        // Add wallclock time stop condition (60 seconds)
        auto wallclock_limit = [start_time = std::chrono::steady_clock::now()](auto &ptc, auto h) {
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - start_time).count();
            return elapsed > 60;
        };
        args.add_stop_condition(wallclock_limit);

        //args.add_operation(StepSlice(DefaultWriter(("test/nick_test/KeplerSimple.dat")), 25));
        
        

        

    




        tools::Timer timer;
        timer.start();

        bool failed = false;
        bool external_timeout = false;
        std::string error_msg;

        // Run simulation with external timeout using std::async
        auto simulation_future = std::async(std::launch::async, [&]() {
            try {
                solver.run(args);
            } catch (const std::exception& e) {
                failed = true;
                error_msg = e.what();
            } catch (...) {
                failed = true;
                error_msg = "Unknown error";
            }
        });

        // Wait for simulation with timeout
        auto status = simulation_future.wait_for(std::chrono::seconds(EXTERNAL_TIMEOUT_SECONDS));

        if (status == std::future_status::timeout) {
            external_timeout = true;
            failed = true;
            error_msg = "External timeout - simulation hung";

            // Move the hung future to prevent blocking on destruction
            std::lock_guard<std::mutex> lock(hung_futures_mutex);
            hung_futures.push_back(std::move(simulation_future));
        }

        double elapsed_time = timer.get_time();

        // Log if simulation took longer than 60 seconds, failed, or hit external timeout
        if (elapsed_time > 60.0 || failed) {
            std::lock_guard<std::mutex> lock(log_mutex);

            if (external_timeout) {
                std::cout << "EXTERNAL TIMEOUT (>" << EXTERNAL_TIMEOUT_SECONDS << "s) - ";
            } else if (failed) {
                std::cout << "RUNTIME ERROR - ";
            } else {
                std::cout << "TIMEOUT (>60s) - ";
            }

            std::cout << "q=" << mratio
                      << ", a=" << sma / 1_AU << " AU"
                      << ", e=" << ecc
                      << ", i=" << inclination / 1_deg << " deg"
                      << ", time=" << elapsed_time << "s";

            if (failed && !external_timeout) {
                std::cout << ", error: " << error_msg;
            } else if (external_timeout) {
                std::cout << ", error: " << error_msg;
            }

            std::cout << std::endl;
        }


        
    }
}

int main(int argc, char **argv)
{

    Scalar massratio_min = 0.1;
    Scalar massratio_max = 10;
    Scalar sma_min = 1e-3_AU;
    Scalar sma_max = 10_AU;
    Scalar ecc_min = 0.1;
    Scalar ecc_max = 0.999;  // Must be < 1 for elliptic orbits (e=1 is parabolic)

    // Grid sizes for parameter sweep
    size_t n_q = 3;    // mass ratio
    size_t n_a = 3;    // semi-major axis
    size_t n_e = 5;    // eccentricity

    // Inclination test values to catch edge cases
    std::vector<Scalar> incl_values = {0_deg, 60_deg, 90_deg, 150_deg, 180_deg};

    std::vector<std::array<Scalar, 4>> combinations;

    // Populate combination array with parameter sweep
    // Using log-space for mass ratio and semi-major axis, linear for eccentricity
    for (size_t i_q = 0; i_q < n_q; i_q++) {
        Scalar q = massratio_min * std::pow(massratio_max / massratio_min, static_cast<Scalar>(i_q) / (n_q - 1));

        for (size_t i_a = 0; i_a < n_a; i_a++) {
            Scalar a = sma_min * std::pow(sma_max / sma_min, static_cast<Scalar>(i_a) / (n_a - 1));

            for (size_t i_e = 0; i_e < n_e; i_e++) {
                Scalar e = ecc_min + (ecc_max - ecc_min) * static_cast<Scalar>(i_e) / (n_e - 1);

                for (const auto& incl : incl_values) {
                    combinations.push_back({q, a, e, incl});
                }
            }
        }
    }

    std::cout << "Running parameter sweep with " << combinations.size() << " combinations..." << std::endl;
    std::cout << "Grid: " << n_q << " mass ratios × " << n_a << " SMAs × "
              << n_e << " eccentricities × " << incl_values.size() << " inclinations" << std::endl;

    tf::Executor executor;
    tools::Timer timer;
    timer.start();

    size_t jobs_per = 16;
    size_t num_combos = combinations.size();
    size_t batch_num = 0;
    size_t total_batches = (num_combos + jobs_per - 1) / jobs_per;

    for (size_t i = 0; i < num_combos; i += jobs_per)
    {
        std::vector<std::array<Scalar, 4>> combos = combinations;
        executor.silent_async(job, combos, i, std::min(i + jobs_per, num_combos));
        executor.wait_for_all();

        batch_num++;
        double elapsed = timer.get_time();

        // Get current system time
        auto now = std::chrono::system_clock::now();
        auto now_time_t = std::chrono::system_clock::to_time_t(now);
        auto now_tm = std::localtime(&now_time_t);

        std::cout << "Batch " << batch_num << "/" << total_batches
                  << " complete. Elapsed time: " << elapsed << "s"
                  << " [" << std::put_time(now_tm, "%Y-%m-%d %H:%M:%S") << "]"
                  << std::endl;
    }

    std::cout << "Complete! Tested " << combinations.size() << " total combinations in "
              << timer.get_time() << "s" << std::endl;
    // g++ -std=c++17 -O3 -pthread simulations/topology/binary-single-3debug.cpp -o simulations/binary-single-3debug

    // simulations/binary-single-3demo
    return 0;
}


// Commands to compile and run this simulation for copy+paste purposes:
// g++ -std=c++17 -O3 -pthread test/nick_test/ParamSweepTest.cpp -o test/nick_test/ParamSweepTest
// test/nick_test/ParamSweepTest

//nohup ./test/nick_test/ParamSweepTest > test/nick_test/sweep_results.txt 2>&1 &
//tail -f test/nick_test/sweep_results.txt