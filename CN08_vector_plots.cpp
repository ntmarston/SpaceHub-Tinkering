#include <array>
#include <cmath>
#include <fstream>
#include "SpaceHub/src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace force;

using Type        = Types<double>;
using Particles   = particles::SizeParticles<Type>;
using Particle    = Particles::Particle;
using Vector      = Particles::Vector;
using VectorArray = Particles::VectorArray;

const double G  = consts::G;
const double c  = consts::C;
const double pi = consts::pi;

const double M_smbh = 1e7_Ms;  // central SMBH mass (placeholder until matched disk file is ready)
const double m_body = 10_Ms;   // embedded body mass


// Per-component RMS of the disk acceleration over an orbit, sqrt(<a_c^2>).
// Samples uniformly in mean anomaly (= uniform in time), so this is a true time-RMS.
// We report RMS, not the signed time-average, because the eccentricity- and
// inclination-damping accelerations are odd in true anomaly (a_R ~ -v_r, a_N ~ -v_z):
// their signed orbit-average is ~0 even though they strongly damp e and i (the damping
// enters via the phase-correlated <a_R sin nu>, not <a_R>). RMS is positive-definite,
// so the radial/normal forcing magnitude stays visible. Uses whatever regime toggles
// are currently set on AGNDisk.
std::array<double, 3> avg_accel_RTN(double a, double e, double inc) {
    const int N = 200;
    double sumR = 0, sumT = 0, sumN = 0;
    for (int k = 0; k < N; ++k) {
        double mean_anom = -pi + 2 * pi * k / N;
        double nu = orbit::M_anomaly_to_T_anomaly(mean_anom, e);

        Particle central{M_smbh, 0.0}, body{m_body, 0.0};
        auto orb = orbit::Elliptic(M_smbh, m_body, a, e, inc, 0.0, 0.0, nu);
        orbit::move_particles(orb, body);
        orbit::move_to_COM_frame(central, body);

        Particles ptcs{0.0, orbit::group(central, body)};
        VectorArray acc(ptcs.number(), Vector{0, 0, 0});
        AGNDisk::add_acc_to(ptcs, acc);

        Vector dr = ptcs.pos()[1] - ptcs.pos()[0];
        Vector dv = ptcs.vel()[1] - ptcs.vel()[0];
        Vector r_hat = dr / norm(dr);
        Vector n_hat = cross(dr, dv) / norm(cross(dr, dv));
        Vector t_hat = cross(n_hat, r_hat);

        double aR = dot(acc[1], r_hat);
        double aT = dot(acc[1], t_hat);
        double aN = dot(acc[1], n_hat);
        sumR += aR * aR;
        sumT += aT * aT;
        sumN += aN * aN;
    }
    return {std::sqrt(sumR / N), std::sqrt(sumT / N), std::sqrt(sumN / N)};
}


int main() {
    AGNDisk::init_from_file("SpaceHub/src/interaction/disk_tab/SG_01Edd_1e7Msun.csv");
    AGNDisk::mute_diagnostics = true;
    AGNDisk::use_JM17_calibration = true;  // Type I migration torque calibration: JM17 (true, default here) vs CN08 (false)

    double Rg = 2 * G * M_smbh / c / c;
    double a  = 1e4 * Rg;                       // semi-major axis held at ~1e3 Rg
    double h0 = AGNDisk::interp_all(a).H / a;   // aspect ratio at R, maps e/h,i/h -> physical e,i

    const int    n         = 60;
    const double e_h_max   = 6.0;   // brackets the e/h = 4 switch threshold
    const double i_h_max   = 3.0;   // brackets the i/h = 1.5 switch threshold
    const double i_h_fixed = 0.0;   // held constant during the e/h sweep
    const double e_h_fixed = 0.0;   // held constant during the i/h sweep

    // forced/default toggle triplets for the three fields we evaluate per grid point
    auto force_typeI    = [] { AGNDisk::type_i_eccentricity_max = AGNDisk::type_i_inclination_max = 1e30; };
    auto force_gas_drag = [] { AGNDisk::type_i_eccentricity_max = AGNDisk::type_i_inclination_max = 1e-30; };
    auto natural_switch = [] { AGNDisk::type_i_eccentricity_max = AGNDisk::type_i_inclination_max = 0.0; };

    const std::string OUT = "SpaceHub/test/migration_test/AnalyticalValidation/CN08/";

    // ===================== Sweep over e/h (i/h fixed) =====================
    std::ofstream v1_ecc(OUT + "v1_typeI_drag_ecc.csv");
    std::ofstream v2_ecc(OUT + "v2_switch_ecc.csv");
    v1_ecc << "e_h,typeI_R,typeI_T,typeI_N,drag_R,drag_T,drag_N\n";
    v2_ecc << "e_h,switch_R,switch_T,switch_N\n";

    for (int k = 0; k <= n; ++k) {
        double e_h = e_h_max * k / n;
        double e = e_h * h0, inc = i_h_fixed * h0;

        force_typeI();    auto ti = avg_accel_RTN(a, e, inc);
        force_gas_drag(); auto gd = avg_accel_RTN(a, e, inc);
        natural_switch(); auto sw = avg_accel_RTN(a, e, inc);

        v1_ecc << e_h << "," << ti[0] << "," << ti[1] << "," << ti[2] << ","
                              << gd[0] << "," << gd[1] << "," << gd[2] << "\n";
        v2_ecc << e_h << "," << sw[0] << "," << sw[1] << "," << sw[2] << "\n";
    }

    // ===================== Sweep over i/h (e/h fixed) =====================
    std::ofstream v1_inc(OUT + "v1_typeI_drag_inc.csv");
    std::ofstream v2_inc(OUT + "v2_switch_inc.csv");
    v1_inc << "i_h,typeI_R,typeI_T,typeI_N,drag_R,drag_T,drag_N\n";
    v2_inc << "i_h,switch_R,switch_T,switch_N\n";

    for (int k = 0; k <= n; ++k) {
        double i_h = i_h_max * k / n;
        double e = e_h_fixed * h0, inc = i_h * h0;

        force_typeI();    auto ti = avg_accel_RTN(a, e, inc);
        force_gas_drag(); auto gd = avg_accel_RTN(a, e, inc);
        natural_switch(); auto sw = avg_accel_RTN(a, e, inc);

        v1_inc << i_h << "," << ti[0] << "," << ti[1] << "," << ti[2] << ","
                              << gd[0] << "," << gd[1] << "," << gd[2] << "\n";
        v2_inc << i_h << "," << sw[0] << "," << sw[1] << "," << sw[2] << "\n";
    }

    return 0;
}

// Build (from repo root):
// g++ -std=c++17 -O3 -pthread CN08_vector_plots.cpp -o CN08_vector_plots
// ./CN08_vector_plots
