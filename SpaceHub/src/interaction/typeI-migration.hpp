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
 * @file typeI-migration.hpp
 *
 * Header file for Type I migration and eccentricity damping.
 */

#pragma once

#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "../dev-tools.hpp"
#include "../spacehub-concepts.hpp"
#include "../orbits/orbits.hpp"
using namespace hub::unit;

namespace hub::force
{
    // Set of values extracted for a specific R from the table
    struct DiskRow {
        double R, R_Rg, Tc, rho, P, cs, H, visc, Sigma, Q, grad_T, grad_Sigma, grad_P, gamma, f_thermal;
    };

    // Set of interpolated properties, returned by the interp_all() method
    struct DiskProps {
        double Sigma, H, rho, Tc, cs, grad_T, grad_Sigma, grad_P, gamma, f_thermal;
    };

    class DiskModel
    {
    public:
        constexpr static bool vel_dependent{true};

        static inline std::vector<DiskRow> disk_table; // A list of objects corresponding to single R-indexed rows in the disk file
        static inline bool initialized = false;

        // Force toggle flags (must be set before running solver)
        static inline bool eccDamping_CN06 = true;
        static inline bool migration_Jimenez = true;

        // Convenience aliases
        static inline bool& eccentricity_damping = eccDamping_CN06;
        static inline bool& migration = migration_Jimenez;

        static inline double Rmin;
        static inline double Rmax;

        // Load pre-tabulated disk CSV before running solver
        static void init_from_file(const std::string& filename) {
            load_disk_data(filename);
            if (!disk_table.empty()) {
                Rmin = disk_table.front().R;
                Rmax = disk_table.back().R;
            }
        }
        static void load_disk_data(const std::string& filename) {
            disk_table.clear();
            std::ifstream file(filename);
            if (!file) {
                throw std::runtime_error("DiskModel Error: Cannot open disk file (!file): " + filename);
            }

            std::string line;
            std::getline(file, line); // skip header (column order must match DiskRow fields)

            while (std::getline(file, line)) {
                std::stringstream ss(line);
                DiskRow row;
                char comma;
                row.gamma = 5.0/3.0;      // default if columns missing
                row.f_thermal = 1.0;
                ss >> row.R >> comma >> row.R_Rg >> comma >> row.Tc >> comma
                   >> row.rho >> comma >> row.P >> comma >> row.cs >> comma
                   >> row.H >> comma >> row.visc >> comma >> row.Sigma >> comma
                   >> row.Q >> comma >> row.grad_T >> comma >> row.grad_Sigma >> comma
                   >> row.grad_P >> comma >> row.gamma >> comma >> row.f_thermal;
                disk_table.push_back(row);
            }
            initialized = true;
        }

        // Linear interpolation of each disk property independently as f(R)
        static DiskProps interp_all(double R) {
            if (R <= Rmin || R >= Rmax) return {0, 0, 0, 0, 0, 0, 0, 0, 5.0/3.0, 1.0};

            auto const& t = disk_table;
            auto it = std::lower_bound(t.begin(), t.end(), R,
                [](const DiskRow& row, double r) { return row.R < r; });
            size_t i = (it - t.begin()) - 1;

            double frac = (R - t[i].R) / (t[i+1].R - t[i].R);
            auto lerp = [&](double f0, double f1) { return f0 + (f1 - f0) * frac; };

            return {lerp(t[i].Sigma, t[i+1].Sigma),
                    lerp(t[i].H, t[i+1].H),
                    lerp(t[i].rho, t[i+1].rho),
                    lerp(t[i].Tc, t[i+1].Tc),
                    lerp(t[i].cs, t[i+1].cs),
                    lerp(t[i].grad_T, t[i+1].grad_T),
                    lerp(t[i].grad_Sigma, t[i+1].grad_Sigma),
                    lerp(t[i].grad_P, t[i+1].grad_P),
                    lerp(t[i].gamma, t[i+1].gamma),
                    lerp(t[i].f_thermal, t[i+1].f_thermal)};
        }

        // Sub-keplerian disk velocity corrected for pressure gradient (Armitage eq. 2.30)
        template <typename Vec>
        static Vec disk_v(const Vec &r, double M, double n, double cs) {
            auto rr = sqrt(r.x * r.x + r.y * r.y);
            auto v_k = sqrt(consts::G * M / rr);
            auto v_k2 = v_k * v_k;
            auto cs2 = cs * cs;
            auto corr_factor = sqrt(1 - n * (cs2 / v_k2));
            auto v = v_k * corr_factor;
            return Vec{-v * r.y / rr, v * r.x / rr, 0};
        }
        
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

    };

    template <typename Particles>
    void DiskModel::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration)
    {
        if (!initialized) {
            throw std::runtime_error("DiskModel Error: Not initialized! Call DiskModel::init_from_file() first.");
        }

        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &v = particles.vel();
        auto const &m = particles.mass();
        auto const &r = particles.radius();

        

        for (size_t i = 1; i < num; ++i)
        {
            // Position relative to central mass
            auto dr = p[i] - p[0];
            auto dv = v[i] - v[0];
            double R_cyl = sqrt(dr.x * dr.x + dr.y * dr.y);
            double z = dr.z;
            auto u = consts::G * (m[0] + m[i]);
            auto [a_orb, ecc] = orbit::calc_a_e(u, dr, dv); //yihan built in method
            auto L_vec = cross(dr, dv);
            double incl = acos(L_vec.z / norm(L_vec));
            


            
            if (R_cyl <= Rmin || R_cyl > Rmax) continue; 

            auto props = interp_all(R_cyl);
            double Sigma = props.Sigma, H = props.H, rho = props.rho,
                   Tc = props.Tc, cs = props.cs, grad_T = props.grad_T,
                   grad_Sigma = props.grad_Sigma, grad_P = props.grad_P,
                   gamma = props.gamma, f_thermal = props.f_thermal;

            //double rho = rho_c * exp(-0.5 * (z * z) / (H * H));

            auto v_disk = disk_v(dr, m[0], grad_P, cs);
            auto v_rel = dv - v_disk;
            auto v2 = dot(v_rel, v_rel);
            auto vmag = sqrt(v2);

            if (vmag < 1e-10) continue;

            //Kept this in when I copied the structure over from disk-model.hpp, could probably delete
            auto r_eff = std::max(r[i], consts::G * m[i] / (v2 + cs * cs));
            auto Mach = vmag / cs;
            double Omega_k = sqrt(consts::G * m[0] / (R_cyl * R_cyl * R_cyl));

            double aspect_ratio = H/R_cyl;
            double e_tilde = ecc / aspect_ratio;


            //==========CN06 Eccentricity Damping (Eq. 17, 19)==========
            if (eccDamping_CN06) {
                //This one does not have extrema lining up properly, temp fix by changing coefficients
                //double Q_e = atan(-3.0 * e_tilde) * (2.0 / consts::pi) * 0.45 + 0.55;
                double Q_e = atan(-3.0 * e_tilde) * (2.0 / consts::pi) * 0.9 + 1.0;
                double t_e = (Q_e / 0.78) * (m[0] / m[i]) * (m[0] / (Sigma * a_orb * a_orb)) * pow(aspect_ratio, 4) * (1.0 + 0.25 * pow(e_tilde, 3)) / Omega_k;

                double vdotr = dot(dv, dr);
                double r2 = dot(dr, dr);
                auto accel_e = dr * (-2.0 * vdotr / (r2 * t_e));

                acceleration[i] += accel_e;
                acceleration[0] -= accel_e * (m[i] / m[0]); //I think this is how the scaling should work?
            }

            //==========Type I Migration Torque (Gilbaum+2025 Section 3.1 (uses JM17 lin_tot)==========
            if (migration_Jimenez) {
                double q = m[i] / m[0];
                double h = aspect_ratio;

                // --- jm_lin_tot ---
                double C_L = (-2.34 + 0.1 * grad_Sigma - 1.5 * grad_T) * f_thermal;
                double C_CR = (0.46 - 0.96 * grad_Sigma + 1.8 * grad_T) / gamma;
                double C_I = C_L + C_CR;

                // --- jm_lin_iso ---
                // double C_I = -1.36 - 0.54 * grad_Sigma - 0.5 * grad_T;
                

                // Normalizing torque (Eq. 12): Gamma0 = q^2 * Sigma * R^4 * Omega^2 * h^-3
                double Gamma0 = q * q * Sigma * R_cyl * R_cyl * R_cyl * R_cyl
                              * Omega_k * Omega_k / (h * h * h);

                // Type I torque (Eq. 13)
                double Gamma_I = C_I * h * Gamma0;

                // Torque -> tangential acceleration (Murray & Dermott T̄ component)
                double T_bar = Gamma_I / (m[i] * R_cyl);

                // Tbar * Rcyl/mag(Rcyl)
                auto a_mig = Vec3{-dr.y, dr.x, 0.0} * (T_bar / R_cyl); //Causes problems if this is not disabled for retrograde orbits!

                acceleration[i] += a_mig;
                acceleration[0] -= a_mig * (m[i] / m[0]);
            }

            //==================== Zhu+2019 model for inclined orbits =====================
            //-------------------- Zhu+2019 inclined migration rate -----------------------
            //TODO
            //--------------------- Zhu+2019 inclination damping --------------------------
            //TODO
            
        }
    }

} // namespace hub::force