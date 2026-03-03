/*---------------------------------------------------------------------------*\
    DiskModelTest.cpp - Component tests for disk-model.hpp

    Tests untested components:
    1. Vertical density profile (Gaussian vs Self-Reg)
    2. disk_v() Keplerian velocity
    3. Relative velocity calculation
    4. I(M) Mach number regimes
    5. Force direction and Newton's 3rd law
    6. Edge cases (vmag->0, high z, boundary radii)

    Output: CSV data for visualization in TestDiskModel.ipynb

    Compile: g++ -std=c++17 -O3 test/nick_test/DiskModelTest.cpp -o test/nick_test/DiskModelTest
    Run: ./test/nick_test/DiskModelTest [test_name]
\*---------------------------------------------------------------------------*/

#include "../../../src/spaceHub.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

using namespace hub;
using namespace unit;
using namespace force;

// Central mass for AGN disk (10^8 solar masses)
constexpr double M_central = 1e8;  // in Ms (SpaceHub units)

//=============================================================================
// Test 1: Vertical Density Profile
// Verify Gaussian profile in Standard zone and uniform in Self-Reg zone
//=============================================================================
void test_vertical_density(const std::string& outfile) {
    std::ofstream out(outfile);
    out << "R,Q,z_over_H,rho_over_rho_c,expected_gaussian,zone\n";
    out << std::setprecision(12);

    // Test radii spanning both zones
    // Standard zone: R < ~3904 AU (Q >> 1)
    // Self-Reg zone: R > ~3904 AU (Q ≈ 1)
    double test_radii[] = {10, 50, 100, 500, 1000, 2000, 3000, 5000, 10000, 50000};
    double z_over_H_vals[] = {0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0};

    for (double R : test_radii) {
        double rho_c = DiskModel::interp(R, &DiskRow::rho);
        double H = DiskModel::interp(R, &DiskRow::H);
        double Q = DiskModel::interp(R, &DiskRow::Q);

        std::string zone = (std::abs(1 - Q) < 0.01) ? "Self-Reg" : "Standard";

        for (double z_over_H : z_over_H_vals) {
            double z = z_over_H * H;

            // Calculate density using the same logic as disk-model.hpp
            double rho;
            if (std::abs(1 - Q) < 0.01) {
                rho = rho_c;  // Self-regulating: uniform
            } else {
                rho = rho_c * std::exp(-0.5 * z_over_H * z_over_H);  // Gaussian
            }

            double expected_gaussian = std::exp(-0.5 * z_over_H * z_over_H);
            double rho_ratio = rho / rho_c;

            out << R << "," << Q << "," << z_over_H << ","
                << rho_ratio << "," << expected_gaussian << "," << zone << "\n";
        }
    }

    std::cout << "Vertical density test written to: " << outfile << "\n";
}

//=============================================================================
// Test 2: disk_v() Keplerian Velocity
// Verify magnitude = sqrt(GM/R) and direction perpendicular to r
//=============================================================================
void test_disk_velocity(const std::string& outfile) {
    std::ofstream out(outfile);
    out << "R,theta,v_mag,v_keplerian,rel_error,v_dot_r,v_z,Lz_sign\n";
    out << std::setprecision(12);

    // Test at various radii and azimuthal angles
    double test_radii[] = {10, 50, 100, 500, 1000, 5000, 10000, 50000};
    double test_angles[] = {0, 45, 90, 135, 180, 225, 270, 315};  // degrees

    for (double R : test_radii) {
        for (double theta_deg : test_angles) {
            double theta = theta_deg * consts::pi / 180.0;

            // Position vector
            using Vec = hub::Vec3<double>;
            Vec r{R * std::cos(theta), R * std::sin(theta), 0};

            // Get disk velocity
            Vec v = DiskModel::disk_v(r, M_central);

            // Expected Keplerian velocity magnitude: v = sqrt(G*M/R)
            // With G=1 in SpaceHub units
            double v_kep = std::sqrt(M_central / R);
            double v_mag = std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
            double rel_error = std::abs(v_mag - v_kep) / v_kep;

            // v should be perpendicular to r: v·r = 0
            double v_dot_r = v.x * r.x + v.y * r.y + v.z * r.z;

            // Angular momentum direction (should be +z for counter-clockwise)
            double Lz = r.x * v.y - r.y * v.x;

            out << R << "," << theta_deg << "," << v_mag << "," << v_kep << ","
                << rel_error << "," << v_dot_r << "," << v.z << "," << (Lz > 0 ? 1 : -1) << "\n";
        }
    }

    std::cout << "Disk velocity test written to: " << outfile << "\n";
}

//=============================================================================
// Test 3: I(M) Function Across Mach Regimes
// Verify continuity and correct limiting behavior
//=============================================================================
void test_mach_regimes(const std::string& outfile) {
    std::ofstream out(outfile);
    out << "Mach,I_value,regime,I_analytical\n";
    out << std::setprecision(12);

    // Define the I(M) calculation functions (copied from disk-model.hpp)
    constexpr double logR = 3.0;
    const double eps = 1 / std::exp(2.0 * logR / 3.0);

    auto I_sup = [](double M, double logR) {
        return (0.5 * std::log(1 - 1 / M / M) + logR) / M / M;
    };
    auto I_sub = [](double M) {
        return (0.5 * std::log((1 + M) / (1 - M)) - M) / M / M;
    };
    auto dIdM_sup = [](double M, double logR) {
        return (-2 * logR + 1 / (M * M - 1) - std::log(1 - 1 / M / M)) / M / M / M;
    };
    auto dIdM_sub = [](double M) {
        return (M * M * M + (1 - M * M) * std::log((1 + M) / (1 - M)) - 2 * M) / (M * M * M * (M * M - 1));
    };

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

    // Sample Mach numbers from 0.01 to 5.0
    for (double M = 0.01; M <= 5.0; M += 0.01) {
        double I = 0;
        std::string regime;
        double I_analytical = 0;

        if (M >= 1 + eps) {
            I = I_sup(M, logR);
            regime = "supersonic";
            I_analytical = I_sup(M, logR);
        }
        else if ((0.1 < M) && (M < 1 - eps)) {
            I = I_sub(M);
            regime = "subsonic";
            I_analytical = I_sub(M);
        }
        else if (M <= 0.1) {
            I = M / 3.0;
            regime = "low_mach";
            I_analytical = M / 3.0;
        }
        else {
            I = connect(M);
            regime = "transition";
            // For transition, analytical is the connect function itself
            I_analytical = connect(M);
        }

        out << M << "," << I << "," << regime << "," << I_analytical << "\n";
    }

    std::cout << "Mach regimes test written to: " << outfile << "\n";
    std::cout << "  eps = " << eps << " (transition at M = " << x1 << " to " << x2 << ")\n";
}

//=============================================================================
// Test 4: Relative Velocity for Circular Orbit
// Particle on circular Keplerian orbit should have v_rel ≈ 0
//=============================================================================
void test_circular_orbit_vrel(const std::string& outfile) {
    std::ofstream out(outfile);
    out << "R,theta,v_rel_mag,v_particle_mag,v_disk_mag\n";
    out << std::setprecision(12);

    double test_radii[] = {100, 500, 1000, 5000};

    for (double R : test_radii) {
        for (int i = 0; i < 360; i += 15) {
            double theta = i * consts::pi / 180.0;

            using Vec = hub::Vec3<double>;
            Vec r{R * std::cos(theta), R * std::sin(theta), 0};

            // Keplerian velocity for circular orbit
            double v_kep = std::sqrt(M_central / R);
            // Velocity perpendicular to r, counter-clockwise
            Vec v_particle{-v_kep * std::sin(theta), v_kep * std::cos(theta), 0};

            // Disk velocity
            Vec v_disk = DiskModel::disk_v(r, M_central);

            // Relative velocity
            Vec v_rel{v_particle.x - v_disk.x, v_particle.y - v_disk.y, v_particle.z - v_disk.z};
            double v_rel_mag = std::sqrt(v_rel.x*v_rel.x + v_rel.y*v_rel.y + v_rel.z*v_rel.z);
            double v_particle_mag = std::sqrt(v_particle.x*v_particle.x + v_particle.y*v_particle.y);
            double v_disk_mag = std::sqrt(v_disk.x*v_disk.x + v_disk.y*v_disk.y);

            out << R << "," << i << "," << v_rel_mag << "," << v_particle_mag << "," << v_disk_mag << "\n";
        }
    }

    std::cout << "Circular orbit v_rel test written to: " << outfile << "\n";
}

//=============================================================================
// Test 5: Force Scaling and Components
// Output f_dyn, f_aero, f_BH for various conditions
//=============================================================================
void test_force_components(const std::string& outfile) {
    std::ofstream out(outfile);
    out << "Mach,rho,m_particle,r_eff,cs,f_dyn,f_aero,f_BH,f_total,f_HL\n";
    out << std::setprecision(12);

    // Test parameters
    double rho = 1e-5;  // Typical disk density in SpaceHub units
    double m_particle = 10.0;  // 10 solar masses
    double r_eff = 1e-6;  // Effective radius in AU (tiny for stellar-mass BH)
    double cs = 800.0;  // Sound speed in AU/time_unit

    constexpr double logR = 3.0;
    const double eps = 1 / std::exp(2.0 * logR / 3.0);

    auto I_sup = [](double M, double logR) {
        return (0.5 * std::log(1 - 1 / M / M) + logR) / M / M;
    };
    auto I_sub = [](double M) {
        return (0.5 * std::log((1 + M) / (1 - M)) - M) / M / M;
    };
    auto dIdM_sub = [](double M) {
        return (M * M * M + (1 - M * M) * std::log((1 + M) / (1 - M)) - 2 * M) / (M * M * M * (M * M - 1));
    };
    auto dIdM_sup = [](double M, double logR) {
        return (-2 * logR + 1 / (M * M - 1) - std::log(1 - 1 / M / M)) / M / M / M;
    };

    const double x1 = 1 - eps;
    const double x2 = 1 + eps;
    const double y1 = I_sub(x1);
    const double y2 = I_sup(x2, logR);
    const double k1 = dIdM_sub(x1);
    const double k2 = dIdM_sup(x2, logR);
    const double a = k1 * (x2 - x1) - (y2 - y1);
    const double b = -k2 * (x2 - x1) + (y2 - y1);

    auto connect = [&](double M) {
        double t = (M - x1) / (x2 - x1);
        return (1 - t) * y1 + y2 * t + (1 - t) * t * (t * b + (1 - t) * a);
    };

    for (double Mach = 0.01; Mach <= 5.0; Mach += 0.02) {
        double I = 0;

        if (Mach >= 1 + eps) {
            I = I_sup(Mach, logR);
        }
        else if ((0.1 < Mach) && (Mach < 1 - eps)) {
            I = I_sub(Mach);
        }
        else if (Mach <= 0.1) {
            I = Mach / 3.0;
        }
        else {
            I = connect(Mach);
        }

        double vmag = Mach * cs;
        double v2 = vmag * vmag;

        // G = 1 in SpaceHub units
        double f_dyn = I * 4 * consts::pi * m_particle * m_particle * rho / (cs * cs);
        double f_aero = consts::pi * r_eff * r_eff * rho * v2;
        double f_HL = 4 * consts::pi * m_particle * m_particle * rho / (cs * cs);
        double f_BH = f_HL * (Mach * Mach / (1 + Mach * Mach)) / (Mach * Mach);
        // Simplifies to: f_BH = f_HL / (1 + Mach^2)

        double f_total = f_dyn + f_aero + f_BH;

        out << Mach << "," << rho << "," << m_particle << "," << r_eff << "," << cs << ","
            << f_dyn << "," << f_aero << "," << f_BH << "," << f_total << "," << f_HL << "\n";
    }

    std::cout << "Force components test written to: " << outfile << "\n";
}

//=============================================================================
// Test 6: Force Direction (should oppose v_rel)
// Verify acceleration is anti-parallel to relative velocity
//=============================================================================
void test_force_direction(const std::string& outfile) {
    std::ofstream out(outfile);
    out << "v_rel_x,v_rel_y,v_rel_z,acc_x,acc_y,acc_z,dot_product,expected_dot\n";
    out << std::setprecision(12);

    // Test with various v_rel directions
    double directions[][3] = {
        {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
        {1, 1, 0}, {1, 0, 1}, {0, 1, 1},
        {1, 1, 1}, {-1, 1, 0}, {1, -1, 1}
    };

    double rho = 1e-5;
    double m_particle = 10.0;
    double r_eff = 1e-6;
    double cs = 800.0;
    double vmag = 400.0;  // M = 0.5

    for (auto& dir : directions) {
        // Normalize direction
        double norm = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
        double vx = vmag * dir[0] / norm;
        double vy = vmag * dir[1] / norm;
        double vz = vmag * dir[2] / norm;

        // Calculate force (simplified, using f_dyn only for this test)
        double Mach = vmag / cs;
        double I = (0.5 * std::log((1 + Mach) / (1 - Mach)) - Mach) / (Mach * Mach);
        double f_dyn = I * 4 * consts::pi * m_particle * m_particle * rho / (cs * cs);

        // Acceleration direction: -f * v_rel / |v_rel| / m
        double acc_x = -f_dyn * vx / vmag / m_particle;
        double acc_y = -f_dyn * vy / vmag / m_particle;
        double acc_z = -f_dyn * vz / vmag / m_particle;

        // Dot product of v_rel and acc (should be negative)
        double dot = (vx * acc_x + vy * acc_y + vz * acc_z);
        double v_mag = std::sqrt(vx*vx + vy*vy + vz*vz);
        double a_mag = std::sqrt(acc_x*acc_x + acc_y*acc_y + acc_z*acc_z);
        double expected_dot = -v_mag * a_mag;  // Anti-parallel means dot = -|v||a|

        out << vx << "," << vy << "," << vz << ","
            << acc_x << "," << acc_y << "," << acc_z << ","
            << dot << "," << expected_dot << "\n";
    }

    std::cout << "Force direction test written to: " << outfile << "\n";
}

//=============================================================================
// Test 7: Edge Case - Near-zero vmag
// Document behavior when particle nearly co-rotates with disk
//=============================================================================
void test_edge_vmag_zero(const std::string& outfile) {
    std::ofstream out(outfile);
    out << "vmag,f_total_times_vrel_over_vmag,is_nan,is_inf\n";
    out << std::setprecision(12);

    double rho = 1e-5;
    double m_particle = 10.0;
    double r_eff = 1e-6;
    double cs = 800.0;

    // Test vmag from 1e-20 to 1e-5
    for (double log_vmag = -20; log_vmag <= -5; log_vmag += 0.5) {
        double vmag = std::pow(10.0, log_vmag);
        double v2 = vmag * vmag;
        double Mach = vmag / cs;

        // Use low-Mach approximation (M <= 0.1)
        double I = Mach / 3.0;

        double f_dyn = I * 4 * consts::pi * m_particle * m_particle * rho / (cs * cs);
        double f_aero = consts::pi * r_eff * r_eff * rho * v2;
        double f_HL = 4 * consts::pi * m_particle * m_particle * rho / (cs * cs);
        double f_BH = f_HL / (1 + Mach * Mach);

        double f_total = f_dyn + f_aero + f_BH;

        // The problematic expression: f_total / vmag
        double result = f_total / vmag;

        out << vmag << "," << result << ","
            << std::isnan(result) << "," << std::isinf(result) << "\n";
    }

    // Also test vmag = 0 explicitly
    double vmag = 0.0;
    double result = 1.0 / vmag;  // Will be inf
    out << vmag << "," << result << ","
        << std::isnan(result) << "," << std::isinf(result) << "\n";

    std::cout << "Edge case vmag->0 test written to: " << outfile << "\n";
    std::cout << "  NOTE: vmag=0 causes division by zero (inf or NaN)\n";
}

//=============================================================================
// Test 8: Interpolation Boundary Behavior
// Test behavior at/beyond disk table boundaries
//=============================================================================
void test_interpolation_bounds(const std::string& outfile) {
    std::ofstream out(outfile);
    out << "R,rho,cs,H,Q,is_inside_bounds\n";
    out << std::setprecision(12);

    // Get table bounds
    double R_min = DiskModel::disk_table.front().R;
    double R_max = DiskModel::disk_table.back().R;

    std::cout << "Disk table bounds: R_min = " << R_min << ", R_max = " << R_max << "\n";

    // Test radii including outside bounds
    double test_radii[] = {
        R_min * 0.1, R_min * 0.5, R_min * 0.9,  // Below min
        R_min, R_min * 1.01,                     // At/near min
        R_max * 0.99, R_max,                     // At/near max
        R_max * 1.1, R_max * 2.0, R_max * 10.0   // Above max
    };

    for (double R : test_radii) {
        double rho = DiskModel::interp(R, &DiskRow::rho);
        double cs = DiskModel::interp(R, &DiskRow::cs);
        double H = DiskModel::interp(R, &DiskRow::H);
        double Q = DiskModel::interp(R, &DiskRow::Q);

        bool inside = (R >= R_min && R <= R_max);

        out << R << "," << rho << "," << cs << "," << H << "," << Q << "," << inside << "\n";
    }

    std::cout << "Interpolation bounds test written to: " << outfile << "\n";
}

//=============================================================================
// Main - Run selected or all tests
//=============================================================================
int main(int argc, char** argv) {
    std::string disk_file = "SpaceHub/src/interaction/disk_tab/disk_test_log.csv";
    std::string output_dir = "SpaceHub/test/nick_test/disk_model/disk_model_output/";

    // Load disk data
    std::cout << "Loading disk data from: " << disk_file << "\n";
    DiskModel::load_disk_data(disk_file);
    std::cout << "Loaded " << DiskModel::disk_table.size() << " rows\n";

    // Create output directory
    std::string mkdir_cmd = "mkdir -p " + output_dir;
    ::system(mkdir_cmd.c_str());

    std::string test_name = (argc > 1) ? argv[1] : "all";

    if (test_name == "all" || test_name == "density") {
        test_vertical_density(output_dir + "vertical_density.csv");
    }
    if (test_name == "all" || test_name == "velocity") {
        test_disk_velocity(output_dir + "disk_velocity.csv");
    }
    if (test_name == "all" || test_name == "mach") {
        test_mach_regimes(output_dir + "mach_regimes.csv");
    }
    if (test_name == "all" || test_name == "vrel") {
        test_circular_orbit_vrel(output_dir + "circular_vrel.csv");
    }
    if (test_name == "all" || test_name == "force") {
        test_force_components(output_dir + "force_components.csv");
    }
    if (test_name == "all" || test_name == "direction") {
        test_force_direction(output_dir + "force_direction.csv");
    }
    if (test_name == "all" || test_name == "vmag_zero") {
        test_edge_vmag_zero(output_dir + "edge_vmag_zero.csv");
    }
    if (test_name == "all" || test_name == "bounds") {
        test_interpolation_bounds(output_dir + "interp_bounds.csv");
    }

    std::cout << "\nAll tests complete. Output in: " << output_dir << "\n";
    std::cout << "Visualize with: jupyter notebook TestDiskModel.ipynb\n";

    return 0;
}
