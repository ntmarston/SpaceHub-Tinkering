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
    class DiskMigration
    {
    public:
        struct DiskRow {
            double R, R_Rg, Tc, rho, P, cs, H, visc, Sigma, Q, grad_T, grad_Sigma, grad_P, gamma, f_thermal;
        };

        struct DiskProps {
            double Sigma, H, rho, Tc, cs, grad_T, grad_Sigma, grad_P, gamma, f_thermal;
        };

        constexpr static bool vel_dependent{true};

        static inline std::vector<DiskRow> disk_table;
        static inline bool initialized = false;

        // Force toggle flags (must be set before running solver)
        // Eventually re-work to automatically determine which to use based on orbital parameters, but keep manual overrides for testing purposes
        static inline bool LoweccDamping_CN06 = false;
        static inline bool CN06_ECC_DECOUPLED = false;
        static inline bool CN06_MIG_DECOUPLED = false;
        static inline bool migration_Jimenez = false;
        static inline bool inclined_zhu = false;
        static inline bool ensemble = false;


        // Overrides
        bool enable_migration = true;
        bool enable_i_damping = true;
        bool enable_e_damping = true;
        bool ignore_dynamical_friction = false;

        // details
        double ecc_tol = 0.01; //  below which is considered circular
        double incl_tol = 0.001_deg; // TODO below which is considered in-plane
        double e_max = 0.75; // Maximum eccentricty before switching to dynamical friction
        // Look through CN06 to find a better value of e_max

        

        static inline double Rmin;
        static inline double Rmax;

        /*
        ==================================================================================================================================================
                                            ========================HELPER METHOD IMPLEMENTATIONS======================
        ==================================================================================================================================================
        */
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
                ss >> row.R >> comma >> row.R_Rg >> comma >> row.Tc >> comma
                   >> row.rho >> comma >> row.P >> comma >> row.cs >> comma
                   >> row.H >> comma >> row.visc >> comma >> row.Sigma >> comma
                   >> row.Q >> comma >> row.grad_T >> comma >> row.grad_Sigma >> comma
                   >> row.grad_P;
                // gamma and f_thermal are optional columns — only read if present.
                // Must save defaults before attempting extraction: operator>> overwrites
                // the target to 0 on parse failure (e.g. if the next column is a string
                // like "zone"), destroying the default set above.
                double g_default = row.gamma, f_default = row.f_thermal;
                if (ss >> comma >> row.gamma) {
                    ss >> comma >> row.f_thermal;
                    if (ss.fail()) row.f_thermal = f_default;
                } else {
                    row.gamma = g_default;
                }
                disk_table.push_back(row);
            }
            initialized = true;
        }

        // Linear interpolation of each disk property independently as f(R)
        static DiskProps interp_all(double R) {
            if (std::isnan(R) || R <= Rmin || R >= Rmax) return {0, 0, 0, 0, 0, 0, 0, 0, 5.0/3.0, 1.0};

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
    void DiskMigration::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration)
    {
        if (!initialized) {
            throw std::runtime_error("DiskMigration Error: Not initialized! Call DiskMigration::init_from_file() first.");
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
            double incl = acos(h_vec.z / norm(h_vec));
            double sin_i = sqrt(h_vec.x * h_vec.x + h_vec.y * h_vec.y) / norm(h_vec); //(2.132)^ + (2.133)^2, solve for SinI



            
            if (R_cyl <= Rmin || R_cyl > Rmax) continue; //Out of disk condition


            //===============================LOAD DISK PROPERTIES FROM CSV================================
            auto props = interp_all(R_cyl);
            double Sigma = props.Sigma, H = props.H, rho = props.rho,
                   Tc = props.Tc, cs = props.cs, grad_T = props.grad_T,
                   grad_Sigma = props.grad_Sigma, grad_P = props.grad_P,
                   gamma = props.gamma, f_thermal = props.f_thermal;

            
            
            //=========================CALCULATE SIMPLE PROPERTIES==================================
            

            auto v_disk = disk_v(dr, m[0], grad_P, cs);
            //double rho = rho_c * exp(-0.5 * (z * z) / (H * H)); //Gaussian density profile
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

            //==================================================================================================================================================
            //                  ======FORCE ENSEMBLE SWITCHING MECHANISM======
            //==================================================================================================================================================
                    //TODO these
            // if (ensemble)
            // if ecc < ecc_tol -> circular=true
            // if incl < incl_tol -> in_plane=true
            // if circular && in_plane -> Jimenez
            // elif circular && !in_plane -> Zhu
            // elif !circular && in_plane -> CN06
            // else -> dynamical friction (need to call the disk-model add_acc_to method?)
              //Default everything to false 
                bool in_plane = false;
                
                bool embedded = false;

                bool dynamical_friction = false;

                bool over_eccentric = false;


                if (incl < incl_tol) {
                    in_plane = true;
                } //else remains false
            
                if (sini < aspect_ratio){
                    embedded = true;
                }//else remains false

                if (ecc > e_max){
                    over_eccentric = true; //Pertains to the eccentricity damping cut-off; migration threshold is 1.1H/R
                } //else remains false


                if ((embedded || over_eccentric) && !ignore_dynamical_friction){
                    //Conditions not met for type I migration, switch to dynamical friction model
                    dynamical_friction = true;
                }// else remain false, proceed to migration calculations
                
            

            

            //==================================================================================================================================================
            //                                    ==========CN06 Eccentric co-planar orbits==========
            // The paper actually gives the acceleration vectors for this one (eqs 18/19 in CN06), which is very convenient -> TODO re-derive for higher eccentricity orbits
            // Intended case: orbits which are not circular, but have negligible inclination (they are essentially co-planar with the midplane of the disk)
            // =========================================LOW ECCENTRICITY APPROXIMATION================================================================
            
            // =========================================DECOUPLED ECCENTRICITY DAMPING================================================================
            
            //==> TODO: look at the 3D Cresswell paper and see if this needs a correction for embedded vs in-plane cases
            
            if (enable_e_damping && (embedded && !dynamical_friction)) {

                //-------------CN06 Eccentricity Damping (Eq. 17, 19)-----------------
                //This one does not have extrema lining up properly, temp fix by changing coefficients
                
                double Q_e = atan(-20.0 *(e_tilde-1)) * (2.0 / consts::pi) * 0.45 + 0.55; //Fitting function to handle issue covered in Fairbairn & Rafikov
                double t_e = (Q_e / 0.78) * (m[0] / m[i]) * (m[0] / (Sigma * a_orb * a_orb)) * pow(aspect_ratio, 4) * (1.0 + 0.25 * pow(e_tilde, 3)) / Omega_k;

                double vdotr = dot(dv, dr);
                if (std::abs(vdotr) < 1e-9) {
                    vdotr = std::copysign(1e-9, vdotr);
                }

                double r2 = dot(dr, dr);
                
                double L = norm(h_vec);
                double r_mag = sqrt(r2);
                double T_bar = ecc * ecc * L / (r_mag * (1.0 - ecc * ecc) * t_e);
                double R_bar = -T_bar * L / vdotr;
                auto r_hat = dr * (1.0 / r_mag);
                auto accel_e = r_hat * R_bar + cross(h_vec * (1.0 / L), r_hat) * T_bar;

                acceleration[i] += accel_e;
                acceleration[0] -= accel_e * (m[i] / m[0]); //I think this is how the scaling should work?


            }


            ////==================================================================================================================================================
            //                          ==========Type I Migration Torque (JM17 lin_tot with PL00 eccentric correction factor)==========
            // intended case: embedded orbits with eccentricity less than 1.1H/r 
            
            //==================================================================================================================================================
            if (enable_migration && (embedded && !dynamical_friction)) {

                double q = m[i] / m[0];

                if (ecc >= 1.0999 * aspect_ratio){

                    //PL00 has a singularity at 1.1H/R, and migration stalls beyond this threshold
                    continue;
                }
                // --- jm_lin_tot ---
                double C_L = (-2.34 + 0.1 * grad_Sigma - 1.5 * grad_T) * f_thermal;
                double C_CR = (0.46 - 0.96 * grad_Sigma + 1.8 * grad_T) / gamma;
                double C_I = C_L + C_CR;

                

                // Normalizing torque (Eq. 12): Gamma0 = q^2 * Sigma * R^4 * Omega^2 * h^-3
                double Gamma0 = q * q * Sigma * R_cyl * R_cyl * R_cyl * R_cyl
                              * Omega_k * Omega_k / (aspect_ratio * aspect_ratio * aspect_ratio);

                // Type I torque (Eq. 13)
                double Gamma_I = C_I * aspect_ratio * Gamma0;


                double Gamma_PL00 = Gamma_I *
                    (1.0 - pow(ecc * R_cyl / (1.1 * H), 4)) /
                    (1.0 + pow(ecc * R_cyl / (1.3 * H), 5));

                double mu = m[i] * m[0] / (m[i] + m[0]);
                double h_mag = norm(h_vec);
                double r2 = dot(dr, dr);
                auto a_mig = cross(h_vec, dr) * (Gamma_PL00 / (mu * h_mag * r2));

                acceleration[i] += a_mig;
                acceleration[0] -= (m[i] / m[0]) * a_mig;
            }


            //==================================================================================================================================================
            //                              ==================== Zhu+2019 inclination damping =====================
            // Intended case: mildly inclined orbits (such that at least 90% of the orbit remains embedded in the disk)  
            //==================================================================================================================================================

            //===> TODO read Cresswell 3d paper and see if this needs alterations

            //===> TODO look more closely at when branch2 is selected, determine how this can be offloaded to the swithcing mechanism above

            if (enable_i_damping && (!in_plane && embedded) && !dynamical_friction) {
                //TODO need to add a condition in this loop to defer to dynamical friction if branch2 is activated
                
                if (incl < 1e-10) continue; //Should not get here anyway if auto-switching is enabled
                if(false){ 
                double q_mratio = m[i] / m[0];
                double h_ar = aspect_ratio;  // Alias so I can keep notation consistent (This is H/R)
                double I = incl;              // inclination, already computed above
                double sinI2 = sin(I / 2.0); // Precomputing saves time later
                double sinI  = sin(I);
                double alpha_s = grad_Sigma;  // = -d(ln Sigma)/d(ln r) = Zhu's alpha_s

                //-------------------- Zhu+2019 inclined migration rate (eqs 14-15) ----------------------
                
                // t_mig^{-1} from Zhu eq 15
                double t_mig_inv = Omega_k * q_mratio * (Sigma * R_cyl * R_cyl / m[0])
                                 / (h_ar * h_ar);

                // tau_mig^{-1} from Zhu eq 14 (min of linear theory + dynamical friction)
                double mig_branch1 = (2.7 + 1.1 * alpha_s) * t_mig_inv;
                double mig_branch2 = 8.8 * h_ar * h_ar / (sinI2 * sinI) * t_mig_inv;
                double tau_mig_inv = std::min(mig_branch1, mig_branch2);


                double T_bar_zhu = -0.5 * R_cyl * Omega_k * tau_mig_inv;

                auto a_mig_zhu = typename Particles::Vector{-dr.y, dr.x, 0.0} * (T_bar_zhu / R_cyl);

                acceleration[i] += a_mig_zhu;
                acceleration[0] -= a_mig_zhu * (m[i] / m[0]);


                //Disabled, but left here for now for testing purposes

                }
                //--------------------- Zhu+2019 inclination damping (eqs 16-17) ------
                // t_inc^{-1} from Zhu eq 17
                double t_inc_inv = Omega_k * q_mratio * (Sigma * R_cyl * R_cyl / m[0])
                                 / (h_ar * h_ar * h_ar * h_ar);

                // tau_I^{-1} from R21 (min of linear and dynamical friction branches)
                double inc_branch1 = 0.544 * t_inc_inv;
                double inc_branch2 = 1.46 * h_ar*h_ar*h_ar*h_ar / (sinI2*sinI2*sinI2 * I) * t_inc_inv;
                double tau_I_inv = std::min(inc_branch1, inc_branch2);

                // Line of nodes: n = Z_hat x h_vec = (-h_vec.y, h_vec.x, 0)
                auto n_vec = typename Particles::Vector{-h_vec.y, h_vec.x, 0.0};
                double n_mag = sqrt(n_vec.x * n_vec.x + n_vec.y * n_vec.y);

                auto n_hat = n_vec * (1.0 / n_mag);

                // R22: N_bar = |r x v| / (r_vec . n_hat) * (-I * tau_I_inv)
                double h_mag = norm(h_vec);  // specific angular momentum magnitude
                double r_dot_nhat = dr.x * n_hat.x + dr.y * n_hat.y + dr.z * n_hat.z;

                // Clamp |r_dot_nhat| to H (scale height) to prevent divergence at nodes
                // Below H, the 2D disk-planet interaction formalism breaks down
                double r_dot_nhat_clamped = (r_dot_nhat >= 0)
                    ? std::max(r_dot_nhat, H) : std::min(r_dot_nhat, -H);

                if (H > 1e-10) {
                    double N_bar = (h_mag / r_dot_nhat_clamped) * (-I * tau_I_inv);

                    // N_bar acts along orbit normal: w_hat = h_vec / |h_vec|
                    auto w_hat = h_vec * (1.0 / h_mag);
                    auto accel_inc = w_hat * N_bar;

                    acceleration[i] += accel_inc;
                    acceleration[0] -= accel_inc * (m[i] / m[0]);
                }
            
            }
        }
    }

} // namespace hub::force