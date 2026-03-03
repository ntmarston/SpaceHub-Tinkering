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
#include <cstdlib>
#include <stdexcept>
#include <iostream>

// TEST HEADER: Modified disk-model.hpp with #ifdef guards for isolation testing
// Define ONE of: TEST_NO_REACTION, TEST_INDIVIDUAL_INTERP, TEST_Q_DENSITY
#include "../../SpaceHub/src/dev-tools.hpp"
#include "../../SpaceHub/src/spacehub-concepts.hpp"
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

    //Implementation
    //===================================Initialization, declaration, helpers, etc=====================
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

        static inline double alpha;
        static inline double M;
        static inline double Mdot;
        static inline double Rmin;
        static inline double Rmax;

        //----------Called before running solver to load pre-tabulated disk file to memory--------
        static void init_from_file(const std::string& filename) {
            load_disk_data(filename);
            if (!disk_table.empty()) {
                Rmin = disk_table.front().R;
                print(std::cout << "Rmin:" << Rmin << "\n");
                Rmax = disk_table.back().R;
                print(std::cout << "Rmax:" << Rmin << "\n");
            }
        }
        /**
         * Generate a new disk model by calling the Python disktab CLI,
         * then load the resulting CSV.
         *
         * @param M_msun        Central BH mass in solar masses
         * @param alpha_val     Shakura-Sunyaev viscosity parameter
         * @param mdot_edd_frac Accretion rate as fraction of Eddington
         * @param epsilon       Radiative efficiency (default 0.1)
         * @param r_min_rg      Minimum radius in gravitational radii (default 10)
         * @param r_max_rg      Maximum radius in gravitational radii (default 1e9)
         * @param n_points      Number of radial grid points (default 500)
         */
        static void init_new_model(double M_msun, double alpha_val, double mdot_edd_frac,
                                   double epsilon = 0.1,
                                   double r_min_rg = 10.0, double r_max_rg = 1e9,
                                   int n_points = 500)
        {
            std::string output_path = "src/interaction/disk_tab/generated_disk.csv";

            std::ostringstream cmd;
            cmd << "python3 disktab/disktab.py"
                << " --M_msun " << M_msun
                << " --alpha " << alpha_val
                << " --mdot_edd_frac " << mdot_edd_frac
                << " --epsilon " << epsilon
                << " --r_min_rg " << r_min_rg
                << " --r_max_rg " << r_max_rg
                << " --n_points " << n_points
                << " --output " << output_path;

            int ret = std::system(cmd.str().c_str()); //capture return code
            if (ret != 0) {
                throw std::runtime_error(
                    "DiskModel Error: Python disk generation failed (exit code "
                    + std::to_string(ret) + "). Command: " + cmd.str());
            }

            M = M_msun;
            alpha = alpha_val;

            init_from_file(output_path);
        }
        
        //-----------Helper Function to load data from csv----------------------
        static void load_disk_data(const std::string& filename) {
            disk_table.clear();
            std::ifstream file(filename);
            if (!file) {
                throw std::runtime_error("DiskModel Error: Cannot open disk file (!file): " + filename);
            }

            std::string line;
            std::getline(file, line); // skip header
    //Note this is dependent on the ordering of the columns in the csv file, so if that ever changes this will need to be updated

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

        //-------------Interpolation helper methods-------------------
        static double interp(double R, double DiskRow::*field) { 
            
            
            //NOTE: this interpolation has worked in all cases thus far, however it is a mess and I will need to re-write it at some point
            //Spline interpolation using catmull-rom method for calculating slopes
            //Slightly unconventional, but it handles both log and linear spacing relatively well
            auto const& table = disk_table;
            auto it = std::lower_bound(table.begin(), table.end(), R,
                [](const DiskRow& row, double r) { return row.R < r; });

            size_t i = std::clamp<size_t>(it - table.begin(), 1, table.size() - 2); //I still don't understand this but stack overflow suggested it and it works
            size_t i0 = (i > 1) ? i - 1 : 0;
            size_t i1 = i;
            size_t i2 = i + 1;
            size_t i3 = std::min(i + 2, table.size() - 1);

            double x0 = table[i0].R, x1 = table[i1].R, x2 = table[i2].R, x3 = table[i3].R;
            double y0 = table[i0].*field, y1 = table[i1].*field,
                   y2 = table[i2].*field, y3 = table[i3].*field;

            double m1 = (y2 - y0) / (x2 - x0);
            double m2 = (y3 - y1) / (x3 - x1);
            double h = x2 - x1;

            // Fritsch-Carlson monotonicity constraint (same as interp_all)
            double delta = (y2 - y1) / h;
            if (std::abs(delta) < 1e-30) {
                m1 = m2 = 0.0;
            } else {
                double alpha = m1 / delta;
                double beta  = m2 / delta;
                if (alpha <= 0.0) m1 = 0.0;
                if (beta  <= 0.0) m2 = 0.0;
                double r2 = alpha * alpha + beta * beta;
                if (r2 > 9.0) {
                    double tau = 3.0 / std::sqrt(r2);
                    m1 = tau * alpha * delta;
                    m2 = tau * beta  * delta;
                }
            }

            double t = (R - x1) / h;
            double t2 = t * t, t3 = t2 * t;

            double h00 = 2*t3 - 3*t2 + 1;
            double h10 = t3 - 2*t2 + t;
            double h01 = -2*t3 + 3*t2;
            double h11 = t3 - t2;

            return h00*y1 + h10*h*m1 + h01*y2 + h11*h*m2;
        }

        // Batch interpolation: single binary search for all needed properties -> replaces interp() to reduce number of interpolation calls from 5 to 1 per time step
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

            // Monotone cubic Hermite interpolation (Fritsch-Carlson?).
            // Catmull-Rom slopes are clamped so the interpolant never overshoots
            // the bracket [y1,y2], preventing negative density/cs/H at opacity
            // transitions where the table has multi-order-of-magnitude jumps.
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
                    double alpha = m1 / delta;
                    double beta  = m2 / delta;
                    if (alpha <= 0.0) m1 = 0.0;
                    if (beta  <= 0.0) m2 = 0.0;
                    double r2 = alpha * alpha + beta * beta;
                    if (r2 > 9.0) {
                        double tau = 3.0 / std::sqrt(r2);
                        m1 = tau * alpha * delta;
                        m2 = tau * beta  * delta;
                    }
                }
                return h00*y1 + h10*h*m1 + h01*y2 + h11*h*m2;
            };

#ifdef TEST_GRADP_NO_FC
            // B5: Plain Catmull-Rom for grad_P (no Fritsch-Carlson clamping)
            // Tests whether FC non-smoothness on grad_P is the root cause
            auto spline_plain = [&](double DiskRow::*field) {
                double y0 = table[i0].*field, y1 = table[i1].*field,
                       y2 = table[i2].*field, y3 = table[i3].*field;
                double m1 = (y2 - y0) / (x2 - x0);
                double m2 = (y3 - y1) / (x3 - x1);
                return h00*y1 + h10*h*m1 + h01*y2 + h11*h*m2;
            };
            return {spline(&DiskRow::rho), spline(&DiskRow::cs),
                    spline(&DiskRow::H), spline_plain(&DiskRow::grad_P)};
#else
            return {spline(&DiskRow::rho), spline(&DiskRow::cs),
                    spline(&DiskRow::H), spline(&DiskRow::grad_P)};
#endif
        }


        //----------Calculate corrected sub-keplerian disk velocity--------------
         template <typename Vec>
        static Vec disk_v(Vec &r, double M, double n, double cs) {
            auto rr = sqrt(r.x * r.x + r.y * r.y);
            auto v_k = sqrt(consts::G * M / rr);
            auto v_k2 = v_k * v_k;
            auto cs2 = cs * cs;
            auto corr_factor = sqrt(1-n*(cs2/v_k2));
            auto v = v_k * corr_factor;
            return Vec{-v * r.y / rr, v * r.x / rr, 0};
        }
        
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration); // Still unclear how this C++ syntax works

    };

    template <typename Particles>
    void DiskModel::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration)
    {
        if (!initialized) {
            throw std::runtime_error("DiskModel Error: Not initialized! Call DiskModel::load_disk_data() first.");
        }

        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &v = particles.vel();
        auto const &m = particles.mass();
        auto const &r = particles.radius();


        //===========Dynamical friction stuff===============
        auto I_sup = [](double M, double logR) { return (0.5 * log(1 - 1 / M / M) + logR) / M / M; };
        auto I_sub = [](double M) { return (0.5 * log((1 + M) / (1 - M)) - M) / M / M; };
        auto dIdM_sup = [](double M, double logR) {
            return (-2 * logR + 1 / (M * M - 1) - log(1 - 1 / M / M)) / M / M / M;
        };
        auto dIdM_sub = [](double M) {
            return (M * M * M + (1 - M * M) * log((1 + M) / (1 - M)) - 2 * M) / (M * M * M * (M * M - 1));
        };

        constexpr double logR = 3.0;
        const double eps = 1 / std::exp(2.0 * logR / 3.0);
        const double x1 = 1 - eps;
        const double x2 = 1 + eps;
        const double y1 = I_sub(x1);
        const double y2 = I_sup(x2, logR);
        const double k1 = dIdM_sub(x1);
        const double k2 = dIdM_sup(x2, logR);
        const double a = k1 * (x2 - x1) - (y2 - y1);
        const double b = -k2 * (x2 - x1) + (y2 - y1);

        auto tt = [&](double M) { return (M - x1) / (x2 - x1); };
        auto connect = [&](double M) {
            double t = tt(M);
            return (1 - t) * y1 + y2 * t + (1 - t) * t * (t * b + (1 - t) * a);
        };


        //==============Calculate and apply the forces to each particle=====================
        for (size_t i = 1; i < num; ++i)
        {


            // Calc position relative to the central mass
            auto dr = p[i] - p[0];
            auto dv = v[i] - v[0];
            double R_cyl = sqrt(dr.x * dr.x + dr.y * dr.y);
            double z = dr.z;
            
            if (R_cyl <= Rmin || R_cyl > Rmax) {
                continue; // Skip particles outside the disk model range
            }

            // --- DIAGNOSTIC: NaN check on inputs ---
            //if (std::isnan(dr.x) || std::isnan(dr.y) || std::isnan(dr.z) ||
            //    std::isnan(dv.x) || std::isnan(dv.y) || std::isnan(dv.z)) {
            //    std::cerr << "[DIAG] NaN in input! dr=(" << dr.x << "," << dr.y << "," << dr.z << ")"
            //              << " dv=(" << dv.x << "," << dv.y << "," << dv.z << ")\n"
            //              << "Rmin: " << Rmin << ")\n";
            //}
            //std::cerr << "[DIAG] i=" << i << " t: " << particles.time() << " R_cyl=" << R_cyl <<  " Rmin: " << Rmin << " z=" << z << "\n";
            // --- END DIAGNOSTIC ---

            

            // Interpolate disk properties for this R
#ifdef TEST_INDIVIDUAL_INTERP
            // B3: Use individual interp() calls instead of interp_all()
            double rho_c = interp(R_cyl, &DiskRow::rho);
            double cs = interp(R_cyl, &DiskRow::cs);
            double H = interp(R_cyl, &DiskRow::H);
            double n = interp(R_cyl, &DiskRow::grad_P);
#else
            auto props = interp_all(R_cyl);
            double rho_c = props.rho, cs = props.cs, H = props.H, n = props.grad_P;
#endif

#ifdef TEST_Q_DENSITY
            // B4: Restore Q-dependent vertical density
            double Q_val = interp(R_cyl, &DiskRow::Q);
            double rho;
            if (std::abs(1 - Q_val) < 0.01) {
                rho = rho_c;  // Self-Reg zone: uniform density
            } else {
                rho = rho_c * exp(-0.5 * (z * z) / (H * H));  // Standard zone: Gaussian
            }
#else
            double rho = rho_c * exp(-0.5 * (z * z) / (H * H));
#endif

            

            // --- DIAGNOSTIC: disk properties ---
            //std::cerr << "[DIAG] rho_c=" << rho_c << " cs=" << cs
            //          << " H=" << H << " n=" << n
            //          << " rho(z)=" << rho << "\n";
            // --- END DIAGNOSTIC ---

            /*
            double rho = 0;
            if (abs(1-Q) < 0.01) {
                rho = rho_c;
            }
            else {
                // calculate local density at height z
                rho = rho_c * exp(-0.5 * (z * z) / (H * H));
            }
            */
            

            auto v_disk = disk_v(dr, m[0], n, cs);
            auto v_rel = dv - v_disk;
            auto v2 = dot(v_rel, v_rel);
            auto vmag = sqrt(v2);

            // --- DIAGNOSTIC: velocity and Mach ---
            //std::cerr << "[DIAG] v_disk=(" << v_disk.x << "," << v_disk.y << "," << v_disk.z << ")"
            //          << " dv=(" << dv.x << "," << dv.y << "," << dv.z << ")"
            //          << " v_rel=(" << v_rel.x << "," << v_rel.y << "," << v_rel.z << ")"
            //          << " vmag=" << vmag << "\n";
            // --- END DIAGNOSTIC ---

            // Skip drag when relative velocity is effectively zero (e.g. circular prograde orbit)
            if (vmag < 1e-10) {
                continue;
            }

            auto r_eff = std::max(r[i], consts::G * m[i] / (v2 + cs * cs));
            auto Mach = vmag / cs;

#ifdef TEST_MACH_GUARD
            // B6 FIX: Skip when Mach << 1 (orbit nearly co-rotating with gas)
            // Physically: I(M) = M/3 at M<<1, so force is negligible at M << 1e-3
            // This prevents extremely slow integration after orbital circularization
            // when the sub-Keplerian correction defines a non-zero floor velocity
            if (Mach < 1e-4) {
                continue;
            }
#endif
            
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

            // --- DIAGNOSTIC: Mach regime and I ---
            //std::cerr << "[DIAG] Mach=" << Mach << " I=" << I
            //          << " regime=" << (Mach >= 1+eps ? "supersonic" : Mach > 1-eps ? "connect" : Mach > 0.1 ? "subsonic" : "low") << "\n";
            // --- END DIAGNOSTIC ---

            // Factor of M^2 is consolidated into the I() function to parameterize in terms of M only
            if (enable_dynamical_friction) {
                double f_dyn = I * 4 * consts::pi * consts::G * consts::G * m[i] * m[i] * rho / (cs * cs);
                f_total += f_dyn;
            }

            if (enable_aerodynamic_drag) {
                double f_aero = consts::pi * r_eff * r_eff * rho * v2; //Working as of 11/26/2025, do not break again
                f_total += f_aero;
            }

            if (enable_bondi_hoyle) {
                double f_HL = 4 * consts::pi * consts::G * consts::G * m[i] * m[i] * rho / (cs * cs);
                double f_BH = f_HL * (pow(Mach, 2) / (1 + pow(Mach, 2) ) ) / pow(Mach, 2);
                f_total += f_BH;
            }

            // --- DIAGNOSTIC: force totals ---
            //std::cerr << "[DIAG] f_total=" << f_total
            //          << " acc_on_i=" << f_total/m[i]
            //          << " acc_on_0=" << f_total/m[0] << "\n";
            //std::cerr << "[DIAG] unit_v=(" << v_rel.x/vmag << "," << v_rel.y/vmag << "," << v_rel.z/vmag << ")\n";
            //std::cerr << "---\n";
            // --- END DIAGNOSTIC ---

            acceleration[i] -= f_total * v_rel / vmag / m[i]; // v_rel / vmag always produces a unit vector in theory
#ifndef TEST_NO_REACTION
            acceleration[0] += f_total * v_rel / vmag / m[0];
#endif
        }
    }

} // namespace hub::force