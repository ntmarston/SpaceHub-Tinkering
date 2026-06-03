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
 * 
 * THIS IS THE CURRENT VERSION OF THE DISK FORCE IMPLEMENTATIONS AS OF 4/29/2026
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
    class DiskMigration
    {
    public:
        struct DiskRow {
            double R, R_Rg, Tc, rho, P, cs, H, visc, Sigma, Q, grad_T, grad_Sigma, grad_P, gamma, f_thermal, v_disk;
        };

        struct DiskProps {
            double Sigma, H, rho, Tc, cs, grad_T, grad_Sigma, grad_P, gamma, f_thermal, v_disk;
        };

        constexpr static bool vel_dependent{true};

        static inline std::vector<DiskRow> disk_table;
        static inline bool initialized = false;

  


        // Force-disable overrides — set to true to globally suppress a force.
        // These are user-facing controls for testing and debugging only.
        // Normal operation is governed by the orbital-regime switching logic (classify()).
        //
        // To use in a simulation, set these after including headers but before running:
        //   force::DiskMigration::DISABLE_MIGRATION = true;
        //   force::DiskMigration::DISABLE_E_DAMPING = true;
        //   solver.run(args);  // uses updated settings
        static inline bool DISABLE_MIGRATION          = false;
        static inline bool DISABLE_I_DAMPING          = false;
        static inline bool DISABLE_E_DAMPING          = false;
        static inline bool DISABLE_DYNAMICAL_FRICTION = false;
        static inline bool DISABLE_AERODYNAMIC_DRAG   = false;
        static inline bool DISABLE_BONDI_HOYLE        = false;
        static inline bool ignore_dynamical_friction  = false;  // suppresses the dyn-fric regime flag in classify()
        static inline bool always_use_dynamical_friction = false; // forces DF-only path; bypasses regime switching & all migration/damping
        static inline bool use_cn08_calibration = false; //use 2.7+1.1\beta/2 instead of the JM17 C_I calibration
        static inline bool force_damping = false; //Always use i,e damping and migration, even if out of disk bounds
        static inline bool mute_diagnostics = false;
        //static inline bool use_ppd_powerlaw = false;
        // details
        static inline double ecc_tol = 0.0001; //  below which is considered circular
        static inline double incl_tol = 0.001_deg; // below which is considered in-plane
       

        

        static inline double Rmin;
        static inline double Rmax;
        static inline std::string disk_filename;
        static inline bool diagnostics_printed = false;

        /*
        ==================================================================================================================================================
                                            ========================HELPER METHOD IMPLEMENTATIONS======================
        ==================================================================================================================================================
        */
        static void print_diagnostics() {
            if (mute_diagnostics) return;
            auto on_off = [](bool disabled) { return disabled ? "OFF" : "ON"; };
            std::cout << "\n=== DiskMigration Configuration ===\n";
            std::cout << "  Disk file         : " << disk_filename << "\n";
            std::cout << "  Migration         : " << on_off(DISABLE_MIGRATION)          << "\n";
            std::cout << "  E-damping         : " << on_off(DISABLE_E_DAMPING)          << "\n";
            std::cout << "  I-damping         : " << on_off(DISABLE_I_DAMPING)          << "\n";
            std::cout << "  Dynamical fric.   : " << on_off(DISABLE_DYNAMICAL_FRICTION) << "\n";
            std::cout << "  Aerodynamic drag  : " << on_off(DISABLE_AERODYNAMIC_DRAG)   << "\n";
            std::cout << "  Bondi-Hoyle       : " << on_off(DISABLE_BONDI_HOYLE)        << "\n";
            std::cout << "  Torque cal.       : " << (use_cn08_calibration ? "CN08 (2.7+1.1*beta/2)" : "JM17 C_I") << "\n";
            std::cout << "===================================\n\n";
        }

        // Load pre-tabulated disk CSV before running solver
        static void init_from_file(const std::string& filename) {
            disk_filename = filename;
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
                throw std::runtime_error("DiskMigration Error: Cannot open disk file (!file): " + filename);
            }

            std::string line;
            std::getline(file, line); // skip header (column order must match DiskRow fields)

            while (std::getline(file, line)) {
                std::stringstream ss(line);
                DiskRow row;
                char comma;
                row.gamma = 5.0/3.0;      // default if columns missing
                row.f_thermal = 1.0;
                row.v_disk = 0.0;
                ss >> row.R >> comma >> row.R_Rg >> comma >> row.Tc >> comma
                   >> row.rho >> comma >> row.P >> comma >> row.cs >> comma
                   >> row.H >> comma >> row.visc >> comma >> row.Sigma >> comma
                   >> row.Q >> comma >> row.grad_T >> comma >> row.grad_Sigma >> comma
                   >> row.grad_P;
                // gamma, f_thermal, v_disk are optional columns — only read if present.
                // Must save defaults before attempting extraction: operator>> overwrites
                // the target to 0 on parse failure (e.g. if the next column is a string
                // like "zone"), destroying the default set above.
                double g_default = row.gamma, f_default = row.f_thermal, vd_default = row.v_disk;
                if (ss >> comma >> row.gamma) {
                    ss >> comma >> row.f_thermal;
                    if (ss.fail()) row.f_thermal = f_default;
                } else {
                    row.gamma = g_default;
                }
                if (ss >> comma >> row.v_disk) { /* parsed */ }
                else row.v_disk = vd_default;
                disk_table.push_back(row);
            }
            initialized = true;
        }

        // Linear interpolation of each disk property independently as f(R)
        static DiskProps interp_all(double R) {
            if (std::isnan(R) || R <= Rmin || R >= Rmax) return {0, 0, 0, 0, 0, 0, 0, 0, 5.0/3.0, 1.0, 0.0};

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
                    lerp(t[i].f_thermal, t[i+1].f_thermal),
                    lerp(t[i].v_disk, t[i+1].v_disk)};
        }

        struct OrbitalRegime {
            /*Structure for auto-switching*/
            bool retrograde, embedded, in_plane, dynamical_friction;
        };


        static OrbitalRegime classify(double ecc, double incl, double sin_i, double aspect_ratio) {
            bool retrograde         = (incl > consts::pi / 2.0);
            bool embedded           = (sin_i < aspect_ratio);
            bool in_plane           = (incl < incl_tol);
            bool dynamical_friction = retrograde || ((!embedded || (ecc >= 0.3)) && !ignore_dynamical_friction);
            //Override for testing: bypass disk-bounds and dyn-fric checks, but
            //leave in_plane based on actual inclination so i-damping turns off
            //naturally once incl < incl_tol (otherwise it fights numerical
            //precision indefinitely once inclination is damped below ~rtol).
            if (force_damping){
                retrograde = false;
                embedded = true;
                dynamical_friction = false;
            }
            return {retrograde, embedded, in_plane, dynamical_friction};
        }

        // Sub-keplerian disk velocity corrected for pressure gradient (Armitage eq. 2.30)
        template <typename Vec>
        static Vec disk_v(const Vec &r, double M, double n, double cs) {
            // Add an override method if we want to do disk turbulence as a stochastic velocity fluctuation
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
    void DiskMigration::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration)
    {
        if (!initialized) {
            throw std::runtime_error("DiskMigration Error: Not initialized! Call DiskMigration::init_from_file() first.");
        }
        if (!diagnostics_printed) {
            print_diagnostics();
            diagnostics_printed = true;
        }

        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &v = particles.vel();
        auto const &m = particles.mass();
        auto const &r = particles.radius();

        

        //TODO these
        // if (ensemble)
        // if ecc < ecc_tol -> circular=true
        // if incl < incl_tol -> in_plane=true
        // if circular && in_plane -> Jimenez
        // elif circular && !in_plane -> Zhu
        // elif !circular && in_plane -> CN06
        // else -> dynamical friction (need to call the disk-model add_acc_to method?)


        for (size_t i = 1; i < num; ++i)
        {
            // Position relative to central mass
            auto dr = p[i] - p[0];
            auto dv = v[i] - v[0];
            double R_cyl = sqrt(dr.x * dr.x + dr.y * dr.y);
            double z = dr.z;
            auto u = consts::G * (m[0] + m[i]);
            auto [a_orb, ecc] = orbit::calc_a_e(u, dr, dv); //yihan built in method
            auto h_vec = cross(dr, dv);
            double h_mag = norm(h_vec);
            double incl = acos(h_vec.z / h_mag);
            double sin_i = sqrt(h_vec.x * h_vec.x + h_vec.y * h_vec.y) / h_mag;

            //==================================================================================
            // Dynamical friction functions and statics
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
        //====================================================================

            
            if (R_cyl <= Rmin || R_cyl > Rmax) continue; //Out of disk condition


            //===============================LOAD DISK PROPERTIES FROM CSV================================
            auto props = interp_all(R_cyl);
            double Sigma = props.Sigma, H = props.H, rho = props.rho,
                   Tc = props.Tc, cs = props.cs, grad_T = props.grad_T,
                   grad_Sigma = props.grad_Sigma, grad_P = props.grad_P,
                   gamma = props.gamma, f_thermal = props.f_thermal;

            
            
            //=========================CALCULATE SIMPLE PROPERTIES==================================
            

            double v_disk_speed = props.v_disk;
            auto v_disk = typename Particles::Vector{-v_disk_speed * dr.y / R_cyl,
                                                      v_disk_speed * dr.x / R_cyl, 0.0};
            //double rho = rho_c * exp(-0.5 * (z * z) / (H * H)); //Gaussian density profile
            auto v_rel = dv - v_disk;
            auto v2 = dot(v_rel, v_rel);
            auto cs2 = props.cs * props.cs;
            auto vmag = sqrt(v2);

            if (vmag < 1e-10) continue;
            double Omega_k = sqrt(consts::G * m[0] / (R_cyl * R_cyl * R_cyl));
            double Omega_CN08 = v_disk_speed / R_cyl;
            double aspect_ratio = H / R_cyl;
            double e_tilde = ecc / aspect_ratio;
            double i_h = incl / aspect_ratio;
            double r2 = dot(dr, dr);
            double r_mag = sqrt(r2);
            double q = m[i] / m[0];
            double t_wave = (m[0] / m[i]) * (m[0] / Sigma / a_orb / a_orb) * pow(aspect_ratio, 4) / Omega_k;

            //==================================================================================================================================================
            //                  ======FORCE ENSEMBLE SWITCHING MECHANISM======
            //==================================================================================================================================================
   
                auto [retrograde, embedded, in_plane, dynamical_friction] = classify(ecc, incl, sin_i, aspect_ratio);
                if (always_use_dynamical_friction) {
                    retrograde = false; embedded = false; in_plane = false; dynamical_friction = true;
                }

                
                
            

            

            //==================================================================================================================================================
            //                                    ==========CN06 Eccentric co-planar orbits==========
            // The paper actually gives the acceleration vectors for this one (eqs 18/19 in CN06), which is very convenient -> TODO re-derive for higher eccentricity orbits
            // Intended case: orbits which are not circular, but have negligible inclination (they are essentially co-planar with the midplane of the disk)
            // =========================================LOW ECCENTRICITY APPROXIMATION================================================================
            
            // =========================================DECOUPLED ECCENTRICITY DAMPING================================================================
            
            //==> TODO: look at the 3D Cresswell paper and see if this needs a correction for embedded vs in-plane cases
            
            if (!DISABLE_E_DAMPING && (embedded && !dynamical_friction)) {

                //-------------CN08 Eccentricity Damping-----------------
                
                
                //double Q_e = atan(-20.0 *(e_tilde-1)) * (2.0 / consts::pi) * 0.45 + 0.55; //Fitting function to handle issue covered in Fairbairn & Rafikov
                //double t_e = (Q_e / 0.78) * (m[0] / m[i]) * (m[0] / (Sigma * a_orb * a_orb)) * pow(aspect_ratio, 4) * (1.0 + 0.25 * pow(e_tilde, 3)) / Omega_k;
                double ecc_cn08 = (1 - 0.14 * pow(e_tilde, 2) + 0.06 * pow(e_tilde, 3) + 0.18 * e_tilde * pow(i_h, 2));
                double t_e = (t_wave / 0.780) * ecc_cn08;
                
                
                double vdotr = dot(dv, dr);
                auto r_hat = dr * (1.0 / r_mag);
                

                /*
                These are energy/semi-major axis conserving, not angular momentum-conserving. 
                This prevents the e-damping force from reducing the semi-major axis at all, 
                is only valid if the reduction of semi-major axis from eccentricity is incorporated in the migration component
                
                double T_bar = ecc * ecc * h_mag / (r_mag * (1.0 - ecc * ecc) * t_e); //INJECTS angular momentum at each step to keep sma constant
                double R_bar = -T_bar * h_mag / vdotr;
                
                auto accel_e = r_hat * R_bar + cross(h_vec * (1.0 / h_mag), r_hat) * T_bar;
                */
                
                auto accel_e = r_hat * (-2.0 * vdotr / t_e / r_mag);

                acceleration[i] += accel_e;
                acceleration[0] -= accel_e * q;


            }


            ////==================================================================================================================================================
            //                          ==========Type I Migration Torque (JM17 lin_tot with ?)==========
            // intended case: embedded orbits with eccentricity less than 1.1H/r 
            
            //==================================================================================================================================================
            if (!DISABLE_MIGRATION && (embedded && !dynamical_friction)) {

                // --- jm_lin_tot ---
                double C_L = (-2.34 + 0.1 * grad_Sigma - 1.5 * grad_T) * f_thermal;
                double C_CR = (0.46 - 0.96 * grad_Sigma + 1.8 * grad_T) / gamma;
                double C_I = C_L + C_CR;
                

                if (use_cn08_calibration){
                    // CN08 Eq.14: a_m = -u/t_m, so inward migration needs t_m > 0. Here Γ_I = C_I·h·Γ₀ with Γ₀ > 0, so inward needs C_I < 0 — sign flips relative to the CN08 timescale formula.
                    C_I = -(2.7 + 1.1 * grad_Sigma)/2.0;
                }
                

                // Normalizing torque (Eq. 12): Gamma0 = q^2 * Sigma * R^4 * Omega^2 * h^-3
                double Gamma0 = q * q * Sigma * R_cyl * R_cyl * R_cyl * R_cyl
                              * Omega_k * Omega_k / (aspect_ratio * aspect_ratio * aspect_ratio); //IS THIS THE RIGHT OMEGA???

                // Type I torque (Eq. 13)
                double Gamma_I = C_I * aspect_ratio * Gamma0;

                double P_e = (1.0 + pow(ecc / (2.25 * aspect_ratio), 1.2) + pow(ecc / (2.84 * aspect_ratio), 6.0))
                           / (1.0 - pow(ecc / (2.02 * aspect_ratio), 4.0));
                //$$P(e) = \frac{1 + \left( \frac{e}{2.25 H/r} \right)^{1.2} + \left( \frac{e}{2.84 H/r} \right)^{6}}{1 - \left( \frac{e}{2.02 H/r} \right)^{4}}$$
                
                double f_cn08 = P_e + (P_e/std::abs(P_e)) * (0.070 * i_h + 0.085 * pow(i_h, 4) - 0.080 * e_tilde * i_h * i_h);

                double Gamma_CN08 = Gamma_I / f_cn08;

                double mu = m[i] * m[0] / (m[i] + m[0]);
                auto a_mig = cross(h_vec, dr) * (Gamma_CN08 / (mu * h_mag * r2));

                acceleration[i] += a_mig;
                acceleration[0] -= q * a_mig;
            }


            //==================================================================================================================================================
            //                              ==================== CN08 inclination damping =====================
            // Intended case: mildly inclined orbits (such that at least 90% of the orbit remains embedded in the disk)  
            //==================================================================================================================================================

            if (!DISABLE_I_DAMPING && (!in_plane && embedded) && !dynamical_friction) {
                //TODO need to add a condition in this loop to defer to dynamical friction if branch2 is activated
                
                if (incl < 1e-10) continue; //Should not get here anyway if auto-switching is enabled

                double incl_cn08 = (1 - 0.3 * pow(i_h, 2) + 0.24 * pow(i_h, 3) + 0.14 * pow(e_tilde, 2) * i_h);
                double tau_I = (t_wave/0.544) * incl_cn08;
                double tau_I_inv = 1/tau_I;

                // CN08 Eq. 16: a_i = -(v_z / t_i) k_hat, with k_hat = disk-midplane normal (z_hat).
                typename Particles::Vector accel_inc{0.0, 0.0, -dv.z * tau_I_inv};

                acceleration[i] += accel_inc;
                acceleration[0] -= accel_inc * q;

                /* ---- Previous Gauss-inversion form kept for reference ----
                // Line of nodes: n = Z_hat x h_vec = (-h_vec.y, h_vec.x, 0)
                auto n_vec = typename Particles::Vector{-h_vec.y, h_vec.x, 0.0};
                double n_mag = sqrt(n_vec.x * n_vec.x + n_vec.y * n_vec.y);
                auto n_hat = n_vec * (1.0 / n_mag);

                // R22: N_bar = |r x v| / (r_vec . n_hat) * (-I * tau_I_inv)
                double r_dot_nhat = dr.x * n_hat.x + dr.y * n_hat.y + dr.z * n_hat.z;

                if (H > 1e-10) {
                // Clamp |r_dot_nhat| to H (scale height) to prevent divergence at nodes
                // Below H, the 2D disk-planet interaction formalism breaks down
                double r_dot_nhat_clamped = (r_dot_nhat >= 0)
                    ? std::max(r_dot_nhat, H) : std::min(r_dot_nhat, -H);
                    double N_bar = (h_mag / r_dot_nhat_clamped) * (-incl * tau_I_inv);

                    // N_bar acts along orbit normal: w_hat = h_vec / |h_vec|
                    auto w_hat = h_vec * (1.0 / h_mag);
                    auto accel_inc_old = w_hat * N_bar;

                    acceleration[i] += accel_inc_old;
                    acceleration[0] -= accel_inc_old * q;
                }
                 */
            }

            
            if (dynamical_friction){

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

                if (!DISABLE_DYNAMICAL_FRICTION) {
                    f_total += I * f_HL;
                }

                if (!DISABLE_AERODYNAMIC_DRAG && !always_use_dynamical_friction) {
                    double f_aero = consts::pi * r_eff * r_eff * rho * v2;
                    f_total += f_aero;
                }

                if (!DISABLE_BONDI_HOYLE && !always_use_dynamical_friction) {
                    f_total += f_HL / (1 + Mach * Mach);
                }

                double inv_vmag = 1.0 / vmag;
                auto drag_acc = f_total * inv_vmag * v_rel;
                acceleration[i] -= drag_acc / m[i];
                acceleration[0] += drag_acc / m[0];
            }
        }
    }

} // namespace hub::force