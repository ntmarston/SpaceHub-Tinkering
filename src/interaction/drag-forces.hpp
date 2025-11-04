/*---------------------------------------------------------------------------*\
        .-''''-.         |
       /        \        |
      /_        _\       |  SpaceHub: The Open Source N-body Toolkit
     // \  <>  / \\      |
     |\__\    /__/|      |  Website:  https://yihanwangastro.github.io/SpaceHub/
      \    ||    /       |
        \  __  /         |  Copyright (C) 2019 Yihan Wang
         '.__.'          |
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
 * @file drag-forces.hpp
 *
 * Header file.
 */
#pragma once

#include <cmath>

#include "../dev-tools.hpp"
#include "../spacehub-concepts.hpp"
using namespace hub::unit;

namespace hub::force
{
    class StaticGasField
    {
    public:
        constexpr static bool vel_dependent{true};
        // Type members
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

        // static double rho;
        // static double cs;
        
    };

    template <typename Particles>
    void StaticGasField::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration)
    {
        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &v = particles.vel();
        auto const &m = particles.mass();
        auto const &r = particles.radius();

        
        // Taken from Yihan's Alpha disk model
        auto I_sup = [](double M, double logR)
        { return (0.5 * log(1 - 1 / M / M) + logR) / M / M; };
        auto I_sub = [](double M)
        { return (0.5 * log((1 + M) / (1 - M)) - M) / M / M; };


        // Taken from Yihan's Alpha disk model blindly
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

        auto tt = [&](double M)
        { return (M - x1) / (x2 - x1); };

        auto connect = [&](double M) //Taylor approx?
        {
            double t = tt(M);
            return (1 - t) * y1 + y2 * t + (1 - t) * t * (t * b + (1 - t) * a);
        };

        for (size_t i = 1; i < num; ++i)
        {
            
            auto dr = p[i] - p[0];
            auto dv = v[i] - v[0];
            auto v_rel = dv;
            auto rho = 1; // assuming kg/m3
            auto v2 = dot(v_rel, v_rel);
            auto vmag = norm(v_rel);
            auto grav_rad = consts::G * m[i]/(pow(vmag, 2));
            auto r_eff = std::max(grav_rad, r[i]); //NOTE: Should always assume the gravitational radius >> physical radius? Valid for compact objects.
            //Note ^ G ~= 1558.5 in the cursed units system (AU, Msun, 1yr=2pi)
            
            auto cs = 0.58625_kms;
            auto Mach = vmag / cs;

            double I = 0;

            double f_total = 0;

            if (Mach > 1 + eps){
                I = (0.5 * log(1 - 1 / (Mach * Mach)) + logR) / (Mach * Mach);
            }
            else if ((0.1 < Mach) && (Mach < 1 - eps)){
                I = (0.5 * log((1 + Mach) / (1 - Mach)) - Mach) / (Mach * Mach);
            }
            else if (Mach < 0.1){
                I = Mach / 3.0;
            }
            else{
                I = connect(Mach);
            }

            auto vesc = sqrt(2*consts::G*m[i] / r[i]);
            

            double f_dyn = I * 4 * consts::pi * consts::G * consts::G * m[i] * m[i] * rho / pow(vmag, 3);
            double f_aero = consts::pi * r_eff * r_eff * rho * vmag;
            double f_bhdrag = 0; //4 * consts::pi * consts::G * consts::G * m[i] * m[i] * rho / pow(vmag, 3);
            
            if (pow(vmag, 2)/pow(vesc,2) <= 1){
                f_total = f_dyn + f_bhdrag;
            }
            else{
                f_total = f_aero + f_bhdrag;
            }
            

            acceleration[i] -= f_total * v_rel / m[i];
            acceleration[0] += f_total * v_rel / m[0];
        }
    }

} // namespace hub::force
