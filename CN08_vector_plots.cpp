#include "SpaceHub/src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;


using f = Interactions<NewtonianGrav, PN1>;
using Solver = methods::Sym6<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;
auto G = consts::G;
auto pi = consts::pi;
auto c = consts::C;


int main() {
    
    //===========Set-up/config========
    //Comment out one of these
    AGNDisk::init_from_file("SpaceHub/src/interaction/disk_tab/SG_01Edd.csv"); //AGN Scale
    //AGNDisk::init_from_file("SpaceHub/src/interaction/disk_tab/cn08_powerlaw_3col.csv"); //PPD scale

    Scalar M = 1e7_Ms;
    Scalar mp = 10_Ms;
    Scalar Rg = 2 * G * M / c / c;
    Scalar R = 1e3 * Rg;

    //============Generate grid============
    //todo iterate over inclination, eccentricity
    //todo from inclination, eccentricity, get everything else required for the vectors (everything that is an input into one of the helpers below)

    //=========V1: No switching between Type I and gas drag==============
    //=============Hold R & inclination constant, iterate over ecc===========

    //Hold R constant at ~1e3 Rg, i_h constant at 0.5, iterate with respect to e_tilde
    //todo for ecc in list of eccentricities
    
    Vec3<Scalar> accel_CN08{}; //contribution to acceleration from all 3 Type I "forces"
    AGNDisk::accel_ecc_damp();
    AGNDisk::accel_migration();
    AGNDisk::accel_inc_damp();
    Vec3<Scalar> accel_gasdrag{}; //contribution to acceleration from gas drag forces
    AGNDisk::accel_gas_drag();
    //todo append to some lists or something to store and plot curves — plot method might happen in python, still tbd. Just load into lists for now

    //=============Hold R & eccentricity constant, iterate over incl===========
    //todo for incl in list of inclinations
    Vec3<Scalar> accel_CN08{};
    AGNDisk::accel_ecc_damp();
    AGNDisk::accel_migration();
    AGNDisk::accel_inc_damp();
    Vec3<Scalar> accel_gasdrag{}; 
    AGNDisk::accel_gas_drag();
    //todo append to some lists or something to store and plot curves


    //=========V3: Switching between Type I and gas drag ENABLED==============

    //=============Hold R & inclination constant, iterate over ecc===========
    //todo for incl in list of inclinations
    Vec3<Scalar> accel_total{}; //Somehow need to get the result of the entire add_acc_to method

    //=============Hold R & eccentricity constant, iterate over incl===========
    //todo for ecc in list of eccentricities
    Vec3<Scalar> accel_total{}; //Somehow need to get the result of the entire add_acc_to method
    



}

// Build (HighFive headers vendored at project root; HDF5 C lib is system-installed):
// g++ -std=c++17 -O3 -pthread -I extern CN08_vector_plots.cpp -lhdf5 -o CN08_vector_plots.cpp (run from repo root)

// ./CN08_vector_plots.cpp
