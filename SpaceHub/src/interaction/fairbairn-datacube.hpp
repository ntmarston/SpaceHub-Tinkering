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
 * @file fairbairn-datacube.hpp
 *
 * Header file for accessing Fairbairn & Rafikov (2025) pre-computed
 * planet-disk interaction data cubes. Provides 4D multilinear interpolation
 * over (q, p, h, e) parameter space for migration torques and damping rates.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace hub::force {

class FairbairnDataCube {
   public:
    static inline bool initialized = false;

    /**
     * Load Fairbairn data from a directory (circular or eccentric) containing parameters.out
     * and data files (e.g. tau_a_inv.out, t_net.out, tau_e_inv.out).
     */
    static void init(const std::string& data_dir) {
        load_parameters(data_dir + "/parameters.out");
        Nq = q_grid.size();
        Np = p_grid.size();
        Nh = h_grid.size();
        Ne = e_grid.size();
        size_t total = Nq * Np * Nh * Ne;

        tau_a_inv_data = load_data_file(data_dir + "/tau_a_inv.out", total);
        t_net_data = load_data_file(data_dir + "/t_net.out", total);
        tau_e_inv_data = load_data_file(data_dir + "/tau_e_inv.out", total);

        initialized = true;
    }

    //------------------Query specific elements------------------

    /**
     * Get tau_a^{-1} (semi-major axis damping rate), dimensionless.
     * Normalized by tau_0^{-1} = F_{J,0} / (M_p * Omega_p * a^2).
     *
     * @param p  Surface density exponent (Sigma ~ r^{-p})
     * @param q  Temperature exponent (T ~ r^{-q})
     * @param h  Disk aspect ratio H/r
     * @param e  Orbital eccentricity
     */
    static double get_tau_a_inv(double p, double q, double h, double e) {
        check_init(); //Make sure class is initialized
        return interp4d(tau_a_inv_data, q, p, h, e);
    }

    /**
     * Query T_net (net torque), dimensionless.
     * Normalized by F_{J,0} = Sigma_p * a^4 * Omega_p^2 * h_p^{-3} * (M_p/M_*)^2.
     */
    static double get_T_net(double p, double q, double h, double e) {
        check_init();
        return interp4d(t_net_data, q, p, h, e);
    }

    /**
     * Query tau_e^{-1} (eccentricity damping rate), dimensionless.
     * Normalized by 100 * tau_0^{-1}.
     */
    static double get_tau_e_inv(double p, double q, double h, double e) {
        check_init();
        return interp4d(tau_e_inv_data, q, p, h, e);
    }


    //----------Initialization/declarations-----------
   private:
    // initialize vectors for grid axes
    static inline std::vector<double> q_grid;
    static inline std::vector<double> p_grid;
    static inline std::vector<double> h_grid;
    static inline std::vector<double> e_grid;

    // Grid dimensions
    static inline size_t Nq, Np, Nh, Ne;

    // Flat data arrays, row-major with shape (Nq, Np, Nh, Ne)
    static inline std::vector<double> tau_a_inv_data;
    static inline std::vector<double> t_net_data;
    static inline std::vector<double> tau_e_inv_data;

    static void check_init() {
        if (!initialized) {
            throw std::runtime_error(
                "FairbairnDataCube not initialized. Call FairbairnDataCube::init() first.");
        }
    }

    // ---- File I/O ----

    //Loads the parameters (which q,p,h, and e values are defined in the directory)
    // from parameters.out
    static void load_parameters(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file) {
            throw std::runtime_error("FairbairnDataCube: Cannot open " + filepath);
        }
        std::string line;
        std::getline(file, line);
        q_grid = parse_line(line);
        std::getline(file, line);
        p_grid = parse_line(line);
        std::getline(file, line);
        h_grid = parse_line(line);
        std::getline(file, line);
        e_grid = parse_line(line);
    }


    //I do not understand c++ File I/O but claude does
    static std::vector<double> parse_line(const std::string& line) {
        std::vector<double> vals;
        std::istringstream ss(line);
        double v;
        while (ss >> v) {
            vals.push_back(v);
        }
        return vals;
    }

    // Load data from one of the data files (e.g. tau_a_inv.out)
    static std::vector<double> load_data_file(const std::string& filepath, size_t expected) {
        std::ifstream file(filepath);
        if (!file) {
            throw std::runtime_error("FairbairnDataCube: Cannot open " + filepath);
        }
        std::vector<double> data;
        data.reserve(expected);
        double v;
        while (file >> v) { //reads one double at a time until the end of the file
            data.push_back(v);
        }
        if (data.size() != expected) {
            throw std::runtime_error("FairbairnDataCube: Expected " + std::to_string(expected) +
                                     " values in " + filepath + ", got " +
                                     std::to_string(data.size()));
        }
        return data;
    }

    // ---- Indexing ----

    static size_t flat_index(size_t iq, size_t ip, size_t ih, size_t ie) {
        return ie + Ne * (ih + Nh * (ip + Np * iq));
    }

    // ---- Interpolation ----
    // ---- This is a placeholder for now, need to re-do all of the interpolation at some point ----

    /**
     * Find the lower bracketing index for x in a sorted grid.
     * Clamped to [0, grid.size()-2] so that idx and idx+1 are always valid.
     */
    static size_t find_interval(const std::vector<double>& grid, double x) {
        auto it = std::upper_bound(grid.begin(), grid.end(), x);
        if (it == grid.begin()) return 0;
        size_t idx = static_cast<size_t>(it - grid.begin()) - 1;
        return std::min(idx, grid.size() - 2);
    }

    /**
     * Compute the fractional weight within the interval [grid[i], grid[i+1]].
     * NOT clamped to [0,1] — values outside that range produce linear extrapolation,
     * matching scipy RegularGridInterpolator with fill_value=None.
     */
    static double weight(const std::vector<double>& grid, size_t i, double x) {
        return (x - grid[i]) / (grid[i + 1] - grid[i]);
    }

    static double lerp(double a, double b, double t) { return a + t * (b - a); }

    /**
     * 4D multilinear interpolation via nested lerps.
     * Axes are (q, p, h, e) matching the data layout.
     * Extrapolates linearly outside grid bounds.
     */
    static double interp4d(const std::vector<double>& data, double q_val, double p_val,
                           double h_val, double e_val) {
        size_t iq = find_interval(q_grid, q_val);
        double tq = weight(q_grid, iq, q_val);
        size_t ip = find_interval(p_grid, p_val);
        double tp = weight(p_grid, ip, p_val);
        size_t ih = find_interval(h_grid, h_val);
        double th = weight(h_grid, ih, h_val);
        size_t ie = find_interval(e_grid, e_val);
        double te = weight(e_grid, ie, e_val);

        // Interpolate along e (innermost), collapsing 16 corners to 8
        double c00 = lerp(data[flat_index(iq, ip, ih, ie)],
                          data[flat_index(iq, ip, ih, ie + 1)], te);
        double c01 = lerp(data[flat_index(iq, ip, ih + 1, ie)],
                          data[flat_index(iq, ip, ih + 1, ie + 1)], te);
        double c10 = lerp(data[flat_index(iq, ip + 1, ih, ie)],
                          data[flat_index(iq, ip + 1, ih, ie + 1)], te);
        double c11 = lerp(data[flat_index(iq, ip + 1, ih + 1, ie)],
                          data[flat_index(iq, ip + 1, ih + 1, ie + 1)], te);
        double c20 = lerp(data[flat_index(iq + 1, ip, ih, ie)],
                          data[flat_index(iq + 1, ip, ih, ie + 1)], te);
        double c21 = lerp(data[flat_index(iq + 1, ip, ih + 1, ie)],
                          data[flat_index(iq + 1, ip, ih + 1, ie + 1)], te);
        double c30 = lerp(data[flat_index(iq + 1, ip + 1, ih, ie)],
                          data[flat_index(iq + 1, ip + 1, ih, ie + 1)], te);
        double c31 = lerp(data[flat_index(iq + 1, ip + 1, ih + 1, ie)],
                          data[flat_index(iq + 1, ip + 1, ih + 1, ie + 1)], te);

        // Interpolate along h, collapsing 8 to 4
        double d0 = lerp(c00, c01, th);
        double d1 = lerp(c10, c11, th);
        double d2 = lerp(c20, c21, th);
        double d3 = lerp(c30, c31, th);

        // Interpolate along p, collapsing 4 to 2
        double e0 = lerp(d0, d1, tp);
        double e1 = lerp(d2, d3, tp);

        // Interpolate along q
        return lerp(e0, e1, tq);
    }
};

}  // namespace hub::force
