/*---------------------------------------------------------------------------*\
    TestDiskVelocity.cpp - Verify disk_v() produces physically reasonable
    sub-Keplerian gas velocities.

    Checks:
    1. |v| is close to Keplerian (within 10%)
    2. |v| differs from Keplerian by more than machine precision
    3. v is perpendicular to r (tangential flow)
    4. v is prograde (positive Lz)
    5. v has no z-component (midplane model)

    Compile: g++ -std=c++17 -O3 -pthread test/nick_test/TestDiskVelocity.cpp -o test/nick_test/TestDiskVelocity
    Run:     ./test/nick_test/TestDiskVelocity   (from project root)
\*---------------------------------------------------------------------------*/

#include "../../../src/spaceHub.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

using namespace hub;
using namespace unit;
using namespace force;

static int failures = 0;

void check(const std::string& name, bool condition) {
    if (condition) {
        std::cout << "  PASS: " << name << "\n";
    } else {
        std::cout << "  FAIL: " << name << "\n";
        failures++;
    }
}

int main() {
    std::string disk_file = "SpaceHub/src/interaction/disk_tab/disk_test_log.csv";
    std::cout << "Loading disk data from: " << disk_file << "\n";
    DiskModel::init_from_file(disk_file);
    std::cout << "Loaded " << DiskModel::disk_table.size() << " rows\n";
    std::cout << "R range: [" << DiskModel::Rmin << ", " << DiskModel::Rmax << "]\n\n";

    constexpr double M_central = 1e8;

    double test_radii[] = {10.0, 100.0, 1000.0, 10000.0};
    double test_angles[] = {0.0, 90.0, 180.0, 270.0};

    using Vec = hub::Vec3<double>;

    for (double R : test_radii) {
        double cs = DiskModel::interp(R, &DiskRow::cs);
        double n  = DiskModel::interp(R, &DiskRow::grad_P);
        double v_kep = std::sqrt(consts::G * M_central / R);

        double ratio = n * cs * cs / (v_kep * v_kep);

        std::cout << std::setprecision(8);
        std::cout << "=== R = " << R << " AU ===\n";
        std::cout << "  cs = " << cs << ", n = " << n
                  << ", v_kep = " << v_kep << "\n";
        std::cout << "  n*cs^2/v_k^2 = " << ratio
                  << ", expected correction = " << std::sqrt(1.0 - ratio) << "\n";

        for (double theta_deg : test_angles) {
            double theta = theta_deg * consts::pi / 180.0;
            Vec r{R * std::cos(theta), R * std::sin(theta), 0.0};

            Vec v = DiskModel::disk_v(r, M_central, n, cs);

            double v_mag = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
            double rel_diff = std::abs(v_mag - v_kep) / v_kep;

            double v_dot_r = v.x * r.x + v.y * r.y + v.z * r.z;
            double cos_angle = v_dot_r / (v_mag * R);

            double Lz = r.x * v.y - r.y * v.x;

            std::string label = "R=" + std::to_string((int)R)
                              + " theta=" + std::to_string((int)theta_deg);

            std::cout << "  --- " << label << " ---\n";
            std::cout << "    |v| = " << v_mag << ", rel_diff = " << rel_diff << "\n";

            check(label + " close to Keplerian (< 10%)", rel_diff < 0.10);
            check(label + " sub-Keplerian (> 1e-12)",    rel_diff > 1e-12);
            check(label + " perpendicular to r",         std::abs(cos_angle) < 1e-12);
            check(label + " prograde (Lz > 0)",          Lz > 0);
            check(label + " v.z == 0",                   v.z == 0.0);
        }
        std::cout << "\n";
    }

    std::cout << "================================\n";
    if (failures == 0) {
        std::cout << "ALL CHECKS PASSED\n";
    } else {
        std::cout << failures << " CHECK(S) FAILED\n";
    }
    std::cout << "================================\n";

    return (failures > 0) ? 1 : 0;
}
