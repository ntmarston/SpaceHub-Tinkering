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
                ss >> row.R >> comma >> row.R_Rg >> comma >> row.Tc >> comma
                   >> row.rho >> comma >> row.P >> comma >> row.cs >> comma
                   >> row.H >> comma >> row.visc >> comma >> row.Sigma >> comma
                   >> row.Q >> comma >> row.grad_T >> comma >> row.grad_Sigma >> comma
                   >> row.grad_P;
                disk_table.push_back(row);
            }
            initialized = true;
        }

        // Catmull-Rom spline with Fritsch-Carlson monotonicity clamping
        // batch-interpolates all needed disk properties at radius R
        // Could probably replace with a simpler method
        static DiskProps interp_all(double R) {
            auto const& table = disk_table;
            auto it = std::lower_bound(table.begin(), table.end(), R,
                [](const DiskRow& row, double r) { return row.R < r; });

            size_t i = std::clamp<size_t>(it - table.begin(), 1, table.size() - 2);
            size_t i0 = (i > 1) ? i - 1 : 0;
            size_t i1 = i, i2 = i + 1;
            size_t i3 = std::min(i + 2, table.size() - 1);

            double x0 = table[i0].R, x1 = table[i1].R, x2 = table[i2].R, x3 = table[i3].R;
            double h = x2 - x1;
            double t = (R - x1) / h;
            double t2 = t * t, t3 = t2 * t;
            double h00 = 2*t3 - 3*t2 + 1, h10 = t3 - 2*t2 + t;
            double h01 = -2*t3 + 3*t2,     h11 = t3 - t2;

            auto spline = [&](double DiskRow::*field) {
                double y0 = table[i0].*field, y1 = table[i1].*field,
                       y2 = table[i2].*field, y3 = table[i3].*field;
                double m1 = (y2 - y0) / (x2 - x0);
                double m2 = (y3 - y1) / (x3 - x1);

                // Fritsch-Carlson monotonicity constraint
                double delta = (y2 - y1) / h;
                if (std::abs(delta) < 1e-30) {
                    m1 = m2 = 0.0;
                } else {
                    double alpha_fc = m1 / delta;
                    double beta_fc  = m2 / delta;
                    if (alpha_fc <= 0.0) m1 = 0.0;
                    if (beta_fc  <= 0.0) m2 = 0.0;
                    double r2 = alpha_fc * alpha_fc + beta_fc * beta_fc;
                    if (r2 > 9.0) {
                        double tau = 3.0 / std::sqrt(r2);
                        m1 = tau * alpha_fc * delta;
                        m2 = tau * beta_fc  * delta;
                    }
                }
                return h00*y1 + h10*h*m1 + h01*y2 + h11*h*m2;
            };

            return {spline(&DiskRow::rho), spline(&DiskRow::cs),
                    spline(&DiskRow::H), spline(&DiskRow::grad_P)};
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

            double rho = rho_c * exp(-0.5 * (z * z) / (H * H));

            auto v_disk = disk_v(dr, m[0], n, cs);
            auto v_rel = dv - v_disk;
            auto v2 = dot(v_rel, v_rel);
            auto vmag = sqrt(v2);

            if (vmag < 1e-10) continue;

            auto r_eff = std::max(r[i], consts::G * m[i] / (v2 + cs * cs));
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

            // M^2 factor is absorbed into I(M) to avoid mixed mach/velocity dependence
            if (enable_dynamical_friction) {
                double f_dyn = I * 4 * consts::pi * consts::G * consts::G * m[i] * m[i] * rho / (cs * cs);
                f_total += f_dyn;
            }

            if (enable_aerodynamic_drag) {
                double f_aero = consts::pi * r_eff * r_eff * rho * v2;
                f_total += f_aero;
            }

            if (enable_bondi_hoyle) {
                double f_HL = 4 * consts::pi * consts::G * consts::G * m[i] * m[i] * rho / (cs * cs);
                double f_BH = f_HL / (1 + Mach * Mach);
                f_total += f_BH;
            }

            acceleration[i] -= f_total * v_rel / vmag / m[i];
            acceleration[0] += f_total * v_rel / vmag / m[0];
        }
    }

} // namespace hub::force