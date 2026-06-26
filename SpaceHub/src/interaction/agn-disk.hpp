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
 * THIS IS THE CURRENT VERSION OF THE DISK FORCE IMPLEMENTATIONS AS OF 6/10/2026
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
    class AGNDisk
    {

        public:
        //Structs
            struct DiskRow {
            double R, R_Rg, Tc, rho, P, cs, H, visc, Sigma, Q, grad_T, grad_Sigma, grad_P, gamma, f_thermal, v_disk;
        };

        struct DiskProps {
            double Sigma, H, rho, Tc, cs, grad_T, grad_Sigma, grad_P, gamma, f_thermal, v_disk;
        };


        constexpr static bool vel_dependent{true};

        static inline std::vector<DiskRow> disk_table;
        static inline bool initialized = false;
        static inline bool diagnostics_printed = false;
        //=================Debug options==========
        static inline bool keplerian_disk_vel = false; //use keplerian velocity instead of sub-kep disk gas velocity

        /*==================OPTIONS===================*/
        static inline bool mute_diagnostics = false;
        static inline bool use_JM17_calibration = false; //use the JM17 C_I calibration instead of 2.7+1.1\beta/2 
        static inline bool secondary_gas_forces = false; // includes contribution from bondi-hoyle-littleton drag and aerodynamic drag

        static inline double type_i_eccentricity_max = 0.0; // >0 overrides the default e < 4*H/R cutoff with this absolute eccentricity
        static inline auto type_i_inclination_max = 0.0_deg; // >0 overrides the default sin(i) < 1.5*H/R cutoff with this absolute inclination angle (radians)
         

        //static inline double migration_max_eccentricity = 0.3; //deprecated, but kept for future testing
        //static inline double migration_max_inclination = 8_deg; //deprecated, but kept for future testing

        static inline double Rmin;
        static inline double Rmax;
        static inline std::string disk_filename;
       

         /*
        ==================================================================================================================================================
                                            ========================HELPER METHOD IMPLEMENTATIONS======================
        ==================================================================================================================================================
        */
        static void print_diagnostics() {
            if (mute_diagnostics) return;
            auto on_off = [](bool disabled) { return disabled ? "OFF" : "ON"; };
            std::cout << "\n=== AGNDisk Class Configuration ===\n";
            std::cout << "  Disk file         : " << disk_filename << "\n";
            std::cout << "  Torque cal.       : " << (use_JM17_calibration ? "JM17 C_I" : "CN08 (2.7+1.1*beta/2)") << "\n";
            std::cout << "===================================\n\n";
        }


        static void load_disk_data(const std::string& filename) {

            // TODO get rid of all the column-guessing stuff and just add a line that reads the column headers and compares it to a hard coded string of what is correct. Throw error if they don't match
            disk_table.clear();
            std::ifstream file(filename);
            if (!file) {
                throw std::runtime_error("AGNDisk Class Error: Cannot open disk file (!file): " + filename);
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
                // gamma, f_thermal, v_disk are "optional" columns
                // Save defaults before attempting load from file: operator >> overwrites
                
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

        // Load pre-tabulated disk CSV before running solver
        static void init_from_file(const std::string& filename) {
            disk_filename = filename;
            load_disk_data(filename);
            if (!disk_table.empty()) {
                Rmin = disk_table.front().R;
                Rmax = disk_table.back().R;
            }
        }

        // Linear interpolation of each disk property f independently as f(R)
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


         // Sub-keplerian disk velocity corrected for pressure gradient (Armitage eq. 2.30)
        template <typename Vec> //My understanding is that this needs to be here because it is called by add_acc_to which uses the template. 
        // It was not compiling before I added this line to make the signatures match
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

        template <typename Vec>
        static Vec accel_ecc_damp(double e_tilde, double i_h, double t_wave,
                                  double vdotr, double r_mag, const Vec &r_hat) {
            double ecc_cn08 = (1 - 0.14 * pow(e_tilde, 2) + 0.06 * pow(e_tilde, 3) + 0.18 * e_tilde * pow(i_h, 2));
            double t_e = (t_wave / 0.780) * ecc_cn08;
            return (-2.0 * vdotr / t_e / r_mag / r_mag) * r_mag * r_hat; //CN08 Eq. 15, todo add dr to sig
        }

        template <typename Vec>
        static Vec accel_inc_damp(double e_tilde, double i_h, double t_wave, double dvz) {
            double incl_corr_factor = (1 - 0.3 * pow(i_h, 2) + 0.24 * pow(i_h, 3) + 0.14 * pow(e_tilde, 2) * i_h);
            double t_i = (t_wave / 0.544) * incl_corr_factor;
            // CN08 Eq. 16: a_i = -(v_z / t_i) k_hat, with k_hat = direction normal to disk plane (here z_hat).
            double accel_i_scalar = -dvz / t_i;
            return Vec{0.0, 0.0, accel_i_scalar};
        }

        template <typename Vec>
        static Vec accel_migration(double grad_Sigma, double grad_T, double f_thermal, double gamma,
                                   double Sigma, double R_cyl, double Omega_k, double aspect_ratio,
                                   double ecc, double e_tilde, double i_h, double q,
                                   double m_i, double m_0, double h_mag, double r2,
                                   const Vec &h_vec, const Vec &dr, const Vec &dv) {
            double C_I;
            if (use_JM17_calibration) {
            double C_L = (-2.34 + 0.1 * grad_Sigma - 1.5 * grad_T) * f_thermal;
            double C_CR = (0.46 - 0.96 * grad_Sigma + 1.8 * grad_T) / gamma;
            C_I = C_L + C_CR;
            }
            else {
                C_I = -(2.7 + 1.1 * grad_Sigma)/2.0; //CN08 calibration
            }

           
            // Normalizing torque (GGS24 Eq. 12): Gamma0 = q^2 * Sigma * R^4 * Omega^2 * h^-3.
            // Omega is the (Keplerian) orbital angular velocity, sqrt(GM/R^3) = Omega_k. Correct as written.
            double Gamma0 = q * q * Sigma * R_cyl * R_cyl * R_cyl * R_cyl
                        * Omega_k * Omega_k / (aspect_ratio * aspect_ratio * aspect_ratio);

            // Net Type I torque (GGS24 Eq. 13): Gamma_I = C_I * (H/R) * Gamma0 = C_I * q^2*Sigma*R^4*Omega^2*h^-2.
            // This is the standard T&W net torque, for C_I = -(2.7+1.1*beta)/2 it equals the CN08 torque.
            // Applied below as a_m = -v/t_m, it reproduces CN08 Eq. 13 exactly: t_m = L/|Gamma_I| = 2*t_wave/(2.7+1.1*beta)*h^-2.
            
            // NB: CN08 and JM17 have different definitions of timescales, which are off by a factor of 2
            double Gamma_I = C_I * aspect_ratio * Gamma0;

            double P_e = (1.0 + pow(ecc / (2.25 * aspect_ratio), 1.2) + pow(ecc / (2.84 * aspect_ratio), 6.0))
                    / (1.0 - pow(ecc / (2.02 * aspect_ratio), 4.0));

            double f_cn08 = P_e + (P_e/std::abs(P_e)) * (0.070 * i_h + 0.085 * pow(i_h, 4) - 0.080 * e_tilde * i_h * i_h);

            double Gamma_CN08 = Gamma_I / f_cn08;

            double mu = m_i * m_0 / (m_i + m_0);
            double j = mu * h_mag;
            return  - dv * abs(Gamma_CN08) / j;
        }

        // Gas drag acceleration on body i: gas dynamical friction (+ optional aerodynamic drag
        // and Bondi-Hoyle-Littleton drag when secondary_gas_forces is enabled). \mathcal{I} taken from Yihan's code
        template <typename Vec>
        static Vec accel_gas_drag(double rho, double cs, double m_i, double r_eff,
                                  double vmag, double v2, double I_factor, double Mach,
                                  const Vec &v_rel) {
            if (vmag < 1e-10) return Vec{0.0, 0.0, 0.0};   // manual guard: prevent ostriker's built in singularity
            double f_HL = 4 * consts::pi * consts::G * consts::G * m_i * m_i * rho / cs /cs;
            double f_total = I_factor * f_HL;              // ostriker gas dynamical friction
            if (secondary_gas_forces) {
                f_total += consts::pi * r_eff * r_eff * rho * v2;   // aerodynamic drag
                f_total += f_HL / (1 + Mach * Mach);                // Bondi-Hoyle-Littleton drag
            }
            return (-f_total / vmag / m_i) * v_rel;       
        }

        
        // Embedded (Type I) regime test, angle-space: i < i_max AND e < e_max.
        // i_max / e_max default to 1.5*(H/R) and 4*(H/R); the type_i_*_max statics override when > 0.
        // Because incl in [0, pi] and i_max is a small prograde angle, retrograde / i>90 deg orbits
        // are always excluded (they can never be in the Type I regime).
        // todo NOTE THIS CURRENTLY USES THE SMALL ANGLE APPROXIMATION WITH [sin(i) ~ i] AND MIGHT NEED TO BE CHANGED LATER
        static bool is_typeI_regime(double incl, double ecc, double aspect_ratio) {
            double i_max = (type_i_inclination_max > 0.0) ? type_i_inclination_max : 1.5 * aspect_ratio;
            double e_max = (type_i_eccentricity_max > 0.0) ? type_i_eccentricity_max : 4.0 * aspect_ratio;
            return (incl < i_max) && (ecc < e_max);
        }



        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);
    };



    template <typename Particles>
    void AGNDisk::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration)
    {
        if (!initialized) {
            throw std::runtime_error("AGNDisk Class Error: Not initialized! Call AGNDisk::init_from_file() first.");
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


        for (size_t i = 1; i < num; ++i)
        {
            
            auto dr = p[i] - p[0]; // Position relative to the central mass
            auto dv = v[i] - v[0]; // Velocity relative to the central mass
            double vdotr = dot(dv, dr);
            double r2 = dot(dr, dr);
            double r_mag = sqrt(r2);
            auto r_hat = dr * (1.0 / r_mag);
            double R_cyl = sqrt(dr.x * dr.x + dr.y * dr.y);
            double z = dr.z;
            auto u = consts::G * (m[0] + m[i]); //gravitational parameter
            auto [a_orb, ecc] = orbit::calc_a_e(u, dr, dv); //yihan built in method
            auto h_vec = cross(dr, dv); //specific angular momentum
            double h_mag = norm(h_vec); //magnitude of specific angular momentum
            double incl = acos(h_vec.z / h_mag); //inclination


            //==================================================================================
            // Dynamical friction functions and statics - credit Yihan

            //These are velocity dependent and therefore need to be in the for loop
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

            auto props = interp_all(R_cyl);

            double Sigma = props.Sigma, H = props.H, rho = props.rho,
                    Tc = props.Tc, cs = props.cs, grad_T = props.grad_T,
                    grad_Sigma = props.grad_Sigma, grad_P = props.grad_P,
                    gamma = props.gamma, f_thermal = props.f_thermal;


            //=========================CALCULATE PROPERTIES==================================
                

            double v_disk_speed = props.v_disk;
            auto v_disk = typename Particles::Vector{-v_disk_speed * dr.y / R_cyl,
                                                    v_disk_speed * dr.x / R_cyl, 0.0};
            //double rho = rho_c * exp(-0.5 * (z * z) / (H * H)); //Gaussian density profile
            auto v_rel = dv - v_disk;
            auto vrel2 = dot(v_rel, v_rel);
            auto cs2 = props.cs * props.cs;
            auto vrel_mag = sqrt(vrel2);

            double Omega_k = sqrt(consts::G * m[0] / (R_cyl * R_cyl * R_cyl));
            double Omega_CN08 = v_disk_speed / R_cyl;
            double aspect_ratio = H / R_cyl;
            double e_tilde = ecc / aspect_ratio;
            double i_h = incl / aspect_ratio;

            double mass_ratio = m[i] / m[0]; // q << 1
            double t_wave = (m[0] / m[i]) * (m[0] / Sigma / a_orb / a_orb) * pow(aspect_ratio, 4) / Omega_k;

            bool in_typeI_regime = is_typeI_regime(incl, ecc, aspect_ratio);
            //CN08 incorporating gradual switch to dynamical friction at e,i >> h/r allows looser criteria than strict h/r cut-offs


            if (in_typeI_regime) {

                /*
                ========================== Type I Migration & Damping =========================
                - Eccentricity & Inclination damping and corrections taken from Cresswell & Nelson 2008
                - Optional Jimenez & Masset hydrodynamic calibration for AGN scales
                
                */

                /*====================== Eccentricity Damping (CN08) ======================*/

                auto accel_e = accel_ecc_damp(e_tilde, i_h, t_wave, vdotr, r_mag, r_hat);

                acceleration[i] += accel_e;
                //acceleration[0] -= accel_e * mass_ratio; 

                /*====================== Inclination Damping (CN08) ======================*/


                auto accel_inc = accel_inc_damp<typename Particles::Vector>(e_tilde, i_h, t_wave, dv.z);

                acceleration[i] += accel_inc;
                //acceleration[0] -= accel_inc * mass_ratio;

                /*====================== Migration Torque (CN08 or JM17) ======================*/
                auto a_mig = accel_migration(grad_Sigma, grad_T, f_thermal, gamma, Sigma, R_cyl,
                                                Omega_k, aspect_ratio, ecc, e_tilde, i_h, mass_ratio,
                                                m[i], m[0], h_mag, r2, h_vec, dr, dv);

                acceleration[i] += a_mig;
                //acceleration[0] -= mass_ratio * a_mig;

                
            }

            else // if not in type I enabled regime
            {
                if (dr.z > H)
                {
                    rho = 0.0;
                }
                auto r_eff = std::max(r[i], consts::G * m[i] / (vrel2 + cs2));
                auto Mach = vrel_mag / cs;

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

                auto a_drag = accel_gas_drag(rho, cs, m[i], r_eff, vrel_mag, vrel2, I, Mach, v_rel);
                acceleration[i] += a_drag;
                //acceleration[0] -= a_drag * mass_ratio;

            }
            
        
        }
    } //end add_acc_to

    



}