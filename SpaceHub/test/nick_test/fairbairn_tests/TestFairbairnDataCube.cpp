/**
 * TestFairbairnDataCube.cpp
 *
 * Standalone test for the FairbairnDataCube 4D interpolation engine.
 * Outputs CSV files for comparison against the Python reference implementation
 * in disktab/accessDataCube.ipynb.
 *
 * Compile:
 *   g++ -std=c++17 -O3 TestFairbairnDataCube.cpp -o TestFairbairnDataCube
 *
 * Run from the repository root:
 *   test/nick_test/TestFairbairnDataCube
 */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../../../src/interaction/fairbairn-datacube.hpp"

using hub::force::FairbairnDataCube;

// ---------------------------------------------------------------------------
// Minimal disk CSV reader (standalone, does not depend on DiskModel)
// ---------------------------------------------------------------------------
struct DiskRow {
    double R, R_Rg, Tc, rho, P, cs, H, visc, Sigma, Q, grad_T, grad_Sigma, grad_P;
};

std::vector<DiskRow> load_disk_csv(const std::string& filename) {
    std::vector<DiskRow> table;
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open disk file: " + filename);
    }
    std::string line;
    std::getline(file, line);  // skip header
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        DiskRow row;
        char comma;
        ss >> row.R >> comma >> row.R_Rg >> comma >> row.Tc >> comma >> row.rho >> comma >>
            row.P >> comma >> row.cs >> comma >> row.H >> comma >> row.visc >> comma >>
            row.Sigma >> comma >> row.Q >> comma >> row.grad_T >> comma >> row.grad_Sigma >>
            comma >> row.grad_P;
        table.push_back(row);
    }
    return table;
}

// Catmull-Rom spline interpolation (same algorithm as disk-model.hpp)
double disk_interp(const std::vector<DiskRow>& table, double R, double DiskRow::*field) {
    auto it = std::lower_bound(table.begin(), table.end(), R,
                               [](const DiskRow& row, double r) { return row.R < r; });
    size_t i = std::clamp<size_t>(it - table.begin(), 1, table.size() - 2);
    size_t i0 = (i > 1) ? i - 1 : 0;
    size_t i1 = i;
    size_t i2 = i + 1;
    size_t i3 = std::min(i + 2, table.size() - 1);

    double x0 = table[i0].R, x1 = table[i1].R, x2 = table[i2].R, x3 = table[i3].R;
    double y0 = table[i0].*field, y1 = table[i1].*field, y2 = table[i2].*field,
           y3 = table[i3].*field;

    double m1 = (y2 - y0) / (x2 - x0);
    double m2 = (y3 - y1) / (x3 - x1);
    double h = x2 - x1;
    double t = (R - x1) / h;
    double t2 = t * t, t3 = t2 * t;

    double h00 = 2 * t3 - 3 * t2 + 1;
    double h10 = t3 - 2 * t2 + t;
    double h01 = -2 * t3 + 3 * t2;
    double h11 = t3 - t2;

    return h00 * y1 + h10 * h * m1 + h01 * y2 + h11 * h * m2;
}

// ---------------------------------------------------------------------------
// Test 1: Query at exact grid points
// ---------------------------------------------------------------------------
void test_grid_points(const std::string& outfile) {
    std::ofstream out(outfile);
    out << std::setprecision(15);
    out << "q,p,h,e,tau_a_inv,T_net,tau_e_inv\n";

    // Access the grid axes through the public query functions at grid points
    // We know the grids from parameters.out:
    // q: 0.0 0.25 0.5 0.75
    // p: 0.0 0.5 1.0 1.5
    // h: 0.06 0.08 0.10
    // e: 0.0 0.01 ... 0.12
    double qs[] = {0.0, 0.25, 0.5, 0.75};
    double ps[] = {0.0, 0.5, 1.0, 1.5};
    double hs[] = {0.06, 0.08, 0.10};
    double es[] = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
                   0.07, 0.08, 0.09, 0.10, 0.11, 0.12};

    for (double q : qs) {
        for (double p : ps) {
            for (double h : hs) {
                for (double e : es) {
                    double tai = FairbairnDataCube::get_tau_a_inv(p, q, h, e);
                    double tn = FairbairnDataCube::get_T_net(p, q, h, e);
                    double tei = FairbairnDataCube::get_tau_e_inv(p, q, h, e);
                    out << q << "," << p << "," << h << "," << e << "," << tai << "," << tn
                        << "," << tei << "\n";
                }
            }
        }
    }
    std::cout << "  Grid points -> " << outfile << "\n";
}

// ---------------------------------------------------------------------------
// Test 2: Query at interior (midpoint) locations
// ---------------------------------------------------------------------------
void test_interior_points(const std::string& outfile) {
    std::ofstream out(outfile);
    out << std::setprecision(15);
    out << "q,p,h,e,tau_a_inv,T_net,tau_e_inv\n";

    double qs[] = {0.125, 0.375, 0.625};
    double ps[] = {0.25, 0.75, 1.25};
    double hs[] = {0.07, 0.09};
    double es[] = {0.005, 0.015, 0.025, 0.035, 0.045, 0.055,
                   0.065, 0.075, 0.085, 0.095, 0.105, 0.115};

    for (double q : qs) {
        for (double p : ps) {
            for (double h : hs) {
                for (double e : es) {
                    double tai = FairbairnDataCube::get_tau_a_inv(p, q, h, e);
                    double tn = FairbairnDataCube::get_T_net(p, q, h, e);
                    double tei = FairbairnDataCube::get_tau_e_inv(p, q, h, e);
                    out << q << "," << p << "," << h << "," << e << "," << tai << "," << tn
                        << "," << tei << "\n";
                }
            }
        }
    }
    std::cout << "  Interior points -> " << outfile << "\n";
}

// ---------------------------------------------------------------------------
// Test 3: Query outside grid bounds (extrapolation)
// ---------------------------------------------------------------------------
void test_extrapolation(const std::string& outfile) {
    std::ofstream out(outfile);
    out << std::setprecision(15);
    out << "q,p,h,e,tau_a_inv,T_net,tau_e_inv\n";

    // Outside each axis individually
    struct TestPoint {
        double q, p, h, e;
    };
    TestPoint pts[] = {
        // Below q
        {-0.1, 0.5, 0.08, 0.05},
        // Above q
        {0.9, 0.5, 0.08, 0.05},
        // Below p (negative p, common in radiation-dominated zone)
        {0.375, -0.5, 0.08, 0.05},
        {0.375, -1.0, 0.08, 0.05},
        {0.375, -1.2349, 0.08, 0.05},
        // Above p
        {0.375, 2.0, 0.08, 0.05},
        // Below h
        {0.375, 0.5, 0.04, 0.05},
        {0.375, 0.5, 0.0066, 0.05},
        // Above h
        {0.375, 0.5, 0.15, 0.05},
        // Above e (beyond grid max of 0.12)
        {0.375, 0.5, 0.08, 0.15},
        {0.375, 0.5, 0.08, 0.20},
        // Multiple axes outside simultaneously (matching notebook Test 4 radii)
        {0.3751, -1.2349, 0.1467, 0.05},  // R/Rg=50
        {0.6556, 1.4940, 0.0066, 0.05},   // R/Rg=5000
        {0.7474, 1.4994, 0.0663, 0.05},   // R/Rg=5e5
    };

    for (auto& pt : pts) {
        double tai = FairbairnDataCube::get_tau_a_inv(pt.p, pt.q, pt.h, pt.e);
        double tn = FairbairnDataCube::get_T_net(pt.p, pt.q, pt.h, pt.e);
        double tei = FairbairnDataCube::get_tau_e_inv(pt.p, pt.q, pt.h, pt.e);
        out << pt.q << "," << pt.p << "," << pt.h << "," << pt.e << "," << tai << "," << tn
            << "," << tei << "\n";
    }
    std::cout << "  Extrapolation -> " << outfile << "\n";
}

// ---------------------------------------------------------------------------
// Test 4: End-to-end (disk table -> parameter extraction -> data cube query)
// ---------------------------------------------------------------------------
void test_end_to_end(const std::vector<DiskRow>& disk_table, const std::string& outfile) {
    std::ofstream out(outfile);
    out << std::setprecision(15);
    out << "R_Rg,p,q,h,e,tau_a_inv,T_net,tau_e_inv\n";

    // Compute Rg from the first row to convert R/Rg -> R(AU)
    double Rg = disk_table.front().R / disk_table.front().R_Rg;

    double test_R_Rg[] = {50.0, 5000.0, 500000.0};
    double test_e = 0.05;

    for (double R_Rg : test_R_Rg) {
        double R_AU = R_Rg * Rg;

        double p = disk_interp(disk_table, R_AU, &DiskRow::grad_Sigma);
        double q = disk_interp(disk_table, R_AU, &DiskRow::grad_T);
        double H = disk_interp(disk_table, R_AU, &DiskRow::H);
        double h = H / R_AU;

        double tai = FairbairnDataCube::get_tau_a_inv(p, q, h, test_e);
        double tn = FairbairnDataCube::get_T_net(p, q, h, test_e);
        double tei = FairbairnDataCube::get_tau_e_inv(p, q, h, test_e);

        out << R_Rg << "," << p << "," << q << "," << h << "," << test_e << "," << tai << ","
            << tn << "," << tei << "\n";

        std::cout << "    R/Rg=" << R_Rg << ": p=" << std::setprecision(6) << p << " q=" << q
                  << " h=" << h << " -> tau_a_inv=" << tai << " T_net=" << tn << "\n";
        std::cout << std::setprecision(15);
    }
    std::cout << "  End-to-end -> " << outfile << "\n";
}

// ---------------------------------------------------------------------------
// Test 5: Eccentricity sweep at fixed (p, q, h) from 3 disk regions
// ---------------------------------------------------------------------------
void test_ecc_sweep(const std::vector<DiskRow>& disk_table, const std::string& outfile) {
    std::ofstream out(outfile);
    out << std::setprecision(15);
    out << "region,R_Rg,p,q,h,e,tau_a_inv,T_net,tau_e_inv\n";

    double Rg = disk_table.front().R / disk_table.front().R_Rg;

    struct Region {
        const char* name;
        double R_Rg;
    };
    Region regions[] = {
        {"Radiation-dom", 50.0},
        {"Gas-dom", 5000.0},
        {"Q1-SelfReg", 500000.0},
    };

    for (auto& reg : regions) {
        double R_AU = reg.R_Rg * Rg;
        double p = disk_interp(disk_table, R_AU, &DiskRow::grad_Sigma);
        double q = disk_interp(disk_table, R_AU, &DiskRow::grad_T);
        double H = disk_interp(disk_table, R_AU, &DiskRow::H);
        double h = H / R_AU;

        // Sweep e = 0.01 to 0.15 in 15 steps (matching notebook)
        for (int i = 0; i < 15; ++i) {
            double e = 0.01 + i * (0.14 / 14.0);
            double tai = FairbairnDataCube::get_tau_a_inv(p, q, h, e);
            double tn = FairbairnDataCube::get_T_net(p, q, h, e);
            double tei = FairbairnDataCube::get_tau_e_inv(p, q, h, e);
            out << reg.name << "," << reg.R_Rg << "," << p << "," << q << "," << h << "," << e
                << "," << tai << "," << tn << "," << tei << "\n";
        }
    }
    std::cout << "  Eccentricity sweep -> " << outfile << "\n";
}

// ---------------------------------------------------------------------------
int main() {
    std::string fairbairn_dir = "SpaceHub/src/interaction/tables/FairbairnData/eccentric";
    std::string disk_csv = "SpaceHub/src/interaction/disk_tab/disk_test_log.csv";
    std::string output_dir = "SpaceHub/test/nick_test/fairbairn_tests/fairbairn_output";

    // Create output directory
    std::system(("mkdir -p " + output_dir).c_str());

    // Initialize data cube
    std::cout << "Loading Fairbairn data cube from: " << fairbairn_dir << "\n";
    FairbairnDataCube::init(fairbairn_dir);
    std::cout << "  Loaded successfully.\n\n";

    // Load disk table
    std::cout << "Loading disk table from: " << disk_csv << "\n";
    auto disk_table = load_disk_csv(disk_csv);
    std::cout << "  Loaded " << disk_table.size() << " rows.\n\n";

    // Run tests
    std::cout << "Running tests:\n";
    test_grid_points(output_dir + "/grid_points.csv");
    test_interior_points(output_dir + "/interior_points.csv");
    test_extrapolation(output_dir + "/extrapolation.csv");

    std::cout << "\n  End-to-end results:\n";
    test_end_to_end(disk_table, output_dir + "/end_to_end.csv");

    std::cout << "\n";
    test_ecc_sweep(disk_table, output_dir + "/ecc_sweep.csv");

    std::cout << "\nAll tests complete. Output in " << output_dir << "/\n";
    return 0;
}
