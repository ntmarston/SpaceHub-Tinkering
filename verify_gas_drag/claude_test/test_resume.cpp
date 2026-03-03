// Test for nick-tools.hpp: read_last_snapshot for SizeParticle and TideParticle
// Prints the resumed state; compare against the last snapshot in each .dat file.
//
// Compile: g++ -std=c++17 -O3 -I../../SpaceHub/src test_resume.cpp -o test_resume
// Run:     ./test_resume

#include "../../SpaceHub/src/spaceHub.hpp"
#include "../../SpaceHub/src/nick-tools.hpp"

#include <iomanip>
#include <iostream>

using namespace hub;

int main() {
    std::cout << std::scientific << std::setprecision(16);

    // --- SizeParticles test ---
    {
        using Particle = particles::SizeParticle<Vec3<double>>;
        auto state = hub::nick::read_last_snapshot<Particle>("test_size_particles.dat");

        std::cout << "=== SizeParticles ===\n";
        std::cout << "resume_time = " << state.resume_time << "\n";
        for (size_t i = 0; i < state.particles.size(); ++i) {
            auto& p = state.particles[i];
            std::cout << "  particle " << i
                      << "  mass=" << p.mass
                      << "  radius=" << p.radius
                      << "\n  pos=(" << p.pos.x << ", " << p.pos.y << ", " << p.pos.z << ")"
                      << "\n  vel=(" << p.vel.x << ", " << p.vel.y << ", " << p.vel.z << ")\n";
        }
    }

    // --- TideParticles test ---
    {
        using Particle = particles::TideParticle<Vec3<double>>;
        auto state = hub::nick::read_last_snapshot<Particle>("test_tide_particles.dat");

        std::cout << "\n=== TideParticles ===\n";
        std::cout << "resume_time = " << state.resume_time << "\n";
        for (size_t i = 0; i < state.particles.size(); ++i) {
            auto& p = state.particles[i];
            std::cout << "  particle " << i
                      << "  mass=" << p.mass
                      << "  radius=" << p.radius
                      << "  k_AM=" << p.tide_apsidal_const
                      << "  tau_lag=" << p.tide_lag_time
                      << "\n  pos=(" << p.pos.x << ", " << p.pos.y << ", " << p.pos.z << ")"
                      << "\n  vel=(" << p.vel.x << ", " << p.vel.y << ", " << p.vel.z << ")\n";
        }
    }

    return 0;
}
