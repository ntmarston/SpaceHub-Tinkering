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
 * @file disk-model.hpp
 *
 * Header file for tabulated disk model with one-time CSV loading.
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
using namespace hub::unit;

namespace hub::force
{
    // Set of values extracted for a specific R from the table
    struct DiskRow {
        double R, R_Rg, Tc, rho, P, cs, H, visc, Sigma, Q, grad_T, grad_Sigma, grad_P;
    };

    // Set of interpolated properties, returned by the interp_all() method
    struct DiskProps {
        double rho, cs, H, grad_P;
    };

    class DiskModel
    {
    public:
        constexpr static bool vel_dependent{true};

        static inline std::vector<DiskRow> disk_table; // A list of objects corresponding to single R-indexed rows in the disk file
        static inline bool initialized = false;

        // Force toggle flags (must be set before running solver)
        static inline bool enable_dynamical_friction = true;
        static inline bool enable_aerodynamic_drag = true;
        static inline bool enable_bondi_hoyle = true;

        static inline double Rmin;
        static inline double Rmax;
        static inline double log_Rmin;
        static inline double inv_dlogR;

        // Load pre-tabulated disk CSV before running solver
        static void init_from_file(const std::string& filename) {
            load_disk_data(filename);
            if (!disk_table.empty()) {
                Rmin = disk_table.front().R;
                Rmax = disk_table.back().R;
                log_Rmin = log(Rmin);
                double log_Rmax = log(Rmax);
                inv_dlogR = (disk_table.size() - 1) / (log_Rmax - log_Rmin);
            }
        }
        static void load_disk_data(const std::string& filename) {
            disk_table.clear();
            std::ifstream file(filename);
            if (!file) {
                throw std::runtime_error("DiskModel Error: Cannot open disk file (!file): " + filename);
            }

            std::string line;
            std::getline(file, line);
            {
                std::stringstream hss(line);
                std::string col;
                const std::vector<std::string> expected = {
                    "R", "R/Rg", "Tc", "rho", "P", "cs", "H", "visc",
                    "Sigma", "Q", "grad_T", "grad_Sigma", "grad_P"
                };
                for (const auto& exp : expected) {
                    if (!std::getline(hss, col, ',') || col != exp)
                        throw std::runtime_error("malformatted disktab file");
                }
            }
            while (std::getline(file, line)) {
                std::stringstream ss(line);
                DiskRow row;
                char comma;
                ss >> row.R >> comma >> row.R_Rg >> comma >> row.Tc >> comma
                   >> row.rho >> comma >> row.P >> comma >> row.cs >> comma
                   >> row.H >> comma >> row.visc >> comma >> row.Sigma >> comma
                   >> row.Q >> comma >> row.grad_T >> comma >> row.grad_Sigma >> comma
                   >> row.grad_P;
                disk_table.push_back(row);
            }
            initialized = true;
        }

        // Linear interpolation of each disk property independently as f(R)
        static DiskProps interp_all(double R) {
            if (R <= Rmin || R >= Rmax) return {0, 0, 0, 0};

            auto const& t = disk_table;
            // O(1) direct index: exploits perfectly log-spaced R grid.
            // Replaces std::lower_bound binary search (~12 comparisons).
            // To revert: replace these 3 lines with:
            //   auto it = std::lower_bound(t.begin(), t.end(), R,
            //       [](const DiskRow& row, double r) { return row.R < r; });
            //   size_t i = (it - t.begin()) - 1;
            double idx_f = (log(R) - log_Rmin) * inv_dlogR;
            size_t i = static_cast<size_t>(idx_f);
            if (i >= t.size() - 1) i = t.size() - 2;

            double frac = (R - t[i].R) / (t[i+1].R - t[i].R);
            auto lerp = [&](double f0, double f1) { return f0 + (f1 - f0) * frac; };

            return {lerp(t[i].rho, t[i+1].rho), lerp(t[i].cs, t[i+1].cs),
                    lerp(t[i].H, t[i+1].H), lerp(t[i].grad_P, t[i+1].grad_P)};
        }

        // Sub-keplerian disk velocity corrected for pressure gradient (Armitage eq. 2.30)
        template <typename Vec>
        static Vec disk_v(const Vec &r, double R_cyl, double M, double n, double cs) {
            auto v_disk_mag = sqrt(consts::G * M / R_cyl - n * cs * cs);
            return Vec{-v_disk_mag * r.y / R_cyl, v_disk_mag * r.x / R_cyl, 0};
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

        // Dynamical friction I(M) functions and transonic Hermite splice constants (logR=3.0 invariant)
        static auto I_sup = [](double M, double logR) { return (0.5 * log(1 - 1 / M / M) + logR) / M / M; };
        static auto I_sub = [](double M) { return (0.5 * log((1 + M) / (1 - M)) - M) / M / M; };
        static auto dIdM_sup = [](double M, double logR) {
            return (-2 * logR + 1 / (M * M - 1) - log(1 - 1 / M / M)) / M / M / M;
        };
        static auto dIdM_sub = [](double M) {
            return (M * M * M + (1 - M * M) * log((1 + M) / (1 - M)) - 2 * M) / (M * M * M * (M * M - 1));
        };

        static constexpr double logR = 3.0;
        static const double eps = 1.0 / std::exp(2.0 * logR / 3.0);
        static const double x1 = 1 - eps;
        static const double x2 = 1 + eps;
        static const double y1 = I_sub(x1);
        static const double y2 = I_sup(x2, logR);
        static const double k1 = dIdM_sub(x1);
        static const double k2 = dIdM_sup(x2, logR);
        static const double a = k1 * (x2 - x1) - (y2 - y1);
        static const double b = -k2 * (x2 - x1) + (y2 - y1);

        static auto tt = [](double M) { return (M - x1) / (x2 - x1); };
        static auto connect = [](double M) {
            double t = tt(M);
            return (1 - t) * y1 + y2 * t + (1 - t) * t * (t * b + (1 - t) * a);
        };

        for (size_t i = 1; i < num; ++i)
        {
            // Position relative to central mass
            auto dr = p[i] - p[0];
            auto dv = v[i] - v[0];
            double R_cyl = sqrt(dr.x * dr.x + dr.y * dr.y);
            double z = dr.z;
            
            if (R_cyl <= Rmin || R_cyl > Rmax) continue; 

            auto props = interp_all(R_cyl);
            double rho_c = props.rho, cs = props.cs, H = props.H, n = props.grad_P;

            double cs2 = cs * cs;
            double rho = rho_c * exp(-0.5 * (z * z) / (H * H));

            auto v_disk = disk_v(dr, R_cyl, m[0], n, cs);
            auto v_rel = dv - v_disk;
            auto v2 = dot(v_rel, v_rel);
            auto vmag = sqrt(v2);

            if (vmag < 1e-10) continue;

            auto r_eff = std::max(r[i], consts::G * m[i] / (v2 + cs2));
            auto Mach = vmag / cs;

            double f_total = 0;
            double I = 0;

            if (Mach >= 1 + eps) {
                I = (0.5 * log(1 - 1 / (Mach * Mach)) + logR) / (Mach * Mach);
            }
            else if ((0.1 < Mach) && (Mach < 1 - eps)) {
                I = (0.5 * log((1 + Mach) / (1 - Mach)) - Mach) / (Mach * Mach);
            }
            else if (Mach <= 0.1) {
                I = Mach / 3.0;
            }
            else {
                I = connect(Mach);
            }

            // Hoyle-Lyttleton common factor (shared by dynamical friction and Bondi-Hoyle)
            double f_HL = 4 * consts::pi * consts::G * consts::G * m[i] * m[i] * rho / cs2;

            if (enable_dynamical_friction) {
                f_total += I * f_HL;
            }

            if (enable_aerodynamic_drag) {
                double f_aero = consts::pi * r_eff * r_eff * rho * v2;
                f_total += f_aero;
            }

            if (enable_bondi_hoyle) {
                f_total += f_HL / (1 + Mach * Mach);
            }

            double inv_vmag = 1.0 / vmag;
            auto drag_acc = f_total * inv_vmag * v_rel;
            acceleration[i] -= drag_acc / m[i];
            acceleration[0] += drag_acc / m[0];
        }
    }

} // namespace hub::force