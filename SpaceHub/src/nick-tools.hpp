/*---------------------------------------------------------------------------*\
        .-''''-.         |
       /        \        |
      /_        _\       |  SpaceHub: The Open Source N-body Toolkit
     // \  <>  / \\      |
     |\__\    /__/|      |  Website:  https://yihanwangastro.github.io/SpaceHub/
      \    ||    /       |
        \  __  /         |  Copyright (C) 2019 Yihan Wang
         '.__.'          |
                         |  author: Nick Marston
---------------------------------------------------------------------
License
    This file is part of SpaceHub.
    SpaceHub is free software: you can redistribute it and/or modify it under
    the terms of the GPL-3.0 License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GPL-3.0 License
    for more details. You should have received a copy of the GPL-3.0 License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @file nick-tools.hpp
 *
 * Utilities for resuming N-body simulations from DefaultWriter output files.
 */
#pragma once

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "particles/drag-particles.hpp"
#include "particles/finite-size.hpp"
#include "particles/point-particles.hpp"
#include "particles/tide-particles.hpp"
#include "orbits/particle-manip.hpp"

namespace hub::nick {

    // Return type: the resume time and reconstructed particle vector.
    template <typename ParticleType>
    struct ResumeState {
        double resume_time;
        std::vector<ParticleType> particles;
    };

    // Primary template — left undefined on purpose.
    // Using an unsupported particle type produces a clear linker error.
    template <typename ParticleType>
    ParticleType parse_particle_line(const std::string& line, double& out_time);

    // PointParticle: time,id,mass,px,py,pz,vx,vy,vz
    template <>
    inline particles::PointParticle<Vec3<double>>
    parse_particle_line<particles::PointParticle<Vec3<double>>>(const std::string& line, double& out_time) {
        std::stringstream ss(line);
        char c;
        double time, mass, px, py, pz, vx, vy, vz;
        int id;
        ss >> time >> c >> id >> c >> mass >> c
           >> px >> c >> py >> c >> pz >> c
           >> vx >> c >> vy >> c >> vz;
        out_time = time;
        return particles::PointParticle<Vec3<double>>{mass, px, py, pz, vx, vy, vz};
    }

    // SizeParticle: time,id,mass,radius,px,py,pz,vx,vy,vz
    template <>
    inline particles::SizeParticle<Vec3<double>>
    parse_particle_line<particles::SizeParticle<Vec3<double>>>(const std::string& line, double& out_time) {
        std::stringstream ss(line);
        char c;
        double time, mass, radius, px, py, pz, vx, vy, vz;
        int id;
        ss >> time >> c >> id >> c >> mass >> c >> radius >> c
           >> px >> c >> py >> c >> pz >> c
           >> vx >> c >> vy >> c >> vz;
        out_time = time;
        return particles::SizeParticle<Vec3<double>>{mass, radius, px, py, pz, vx, vy, vz};
    }

    // TideParticle: time,id,mass,radius,k_AM,tau_lag,px,py,pz,vx,vy,vz
    template <>
    inline particles::TideParticle<Vec3<double>>
    parse_particle_line<particles::TideParticle<Vec3<double>>>(const std::string& line, double& out_time) {
        std::stringstream ss(line);
        char c;
        double time, mass, radius, k_AM, tau_lag, px, py, pz, vx, vy, vz;
        int id;
        ss >> time >> c >> id >> c >> mass >> c >> radius >> c
           >> k_AM >> c >> tau_lag >> c
           >> px >> c >> py >> c >> pz >> c
           >> vx >> c >> vy >> c >> vz;
        out_time = time;
        return particles::TideParticle<Vec3<double>>{mass, radius, k_AM, tau_lag, px, py, pz, vx, vy, vz};
    }

    // DragParticle: time,id,mass,px,py,pz,vx,vy,vz,sub_sonic_c,local_cs
    // Note: CSV column order differs from constructor order (mass, sub_sonic_c, cs, px...).
    template <>
    inline particles::DragParticle<Vec3<double>>
    parse_particle_line<particles::DragParticle<Vec3<double>>>(const std::string& line, double& out_time) {
        std::stringstream ss(line);
        char c;
        double time, mass, px, py, pz, vx, vy, vz, sub_sonic_c, local_cs;
        int id;
        ss >> time >> c >> id >> c >> mass >> c
           >> px >> c >> py >> c >> pz >> c
           >> vx >> c >> vy >> c >> vz >> c
           >> sub_sonic_c >> c >> local_cs;
        out_time = time;
        return particles::DragParticle<Vec3<double>>{mass, sub_sonic_c, local_cs, px, py, pz, vx, vy, vz};
    }

    // Read the last complete snapshot from a DefaultWriter output file.
    // Returns the resume time and the reconstructed particle vector.
    template <typename ParticleType>
    ResumeState<ParticleType> read_last_snapshot(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("nick::read_last_snapshot: cannot open file: " + filepath);
        }

        std::string line;
        std::getline(file, line);  // skip header

        std::vector<ParticleType> last_snapshot, current_snapshot;
        double current_time = 0.0, last_time = 0.0;

        while (std::getline(file, line)) {
            if (line.empty()) {
                if (!current_snapshot.empty()) {
                    last_snapshot = std::move(current_snapshot);
                    last_time = current_time;
                    current_snapshot.clear();
                }
                continue;
            }
            double t;
            current_snapshot.push_back(parse_particle_line<ParticleType>(line, t));
            current_time = t;
        }
        // Handle files that don't end with a blank line
        if (!current_snapshot.empty()) {
            last_snapshot = std::move(current_snapshot);
            last_time = current_time;
        }

        if (last_snapshot.empty()) {
            throw std::runtime_error("nick::read_last_snapshot: no data found in " + filepath);
        }

        return ResumeState<ParticleType>{last_time, std::move(last_snapshot)};
    }

    // Convenience wrapper: reads the last snapshot and moves particles to the COM frame.
    // Use this instead of read_last_snapshot when you want to resume cleanly.
    template <typename ParticleType>
    ResumeState<ParticleType> read_last_snapshot_com(const std::string& filepath) {
        auto state = read_last_snapshot<ParticleType>(filepath);
        orbit::move_to_COM_frame(state.particles);
        return state;
    }

}  // namespace hub::nick
