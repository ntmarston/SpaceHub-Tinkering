// Minimal scratch test for AGNDisk::add_acc_to.
// Not part of the final product — used to explore unit test patterns.
//
// Compile and run from repo root:
//   g++ -std=c++17 -O0 -pthread SpaceHub/test/agn-disk-test.cpp -o SpaceHub/test/agn-disk-test
//   SpaceHub/test/agn-disk-test

#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "../src/type-class.hpp"
#include "../src/particles/finite-size.hpp"
#include "../src/orbits/orbits.hpp"
#include "../src/orbits/particle-manip.hpp"
#include "../src/interaction/agn-disk.hpp"

using namespace hub;
using namespace hub::unit;
using namespace hub::force;

using Type        = hub::Types<double>;
using Particles   = hub::particles::SizeParticles<Type>;
using Particle    = Particles::Particle;
using Vector      = Particles::Vector;
using VectorArray = Type::VectorArray;

static double M, m_sbh;
static Particle central_smbh, sbh;
static bool verbose = true;
std::vector<bool> test_load_disk_data();
std::vector<bool> test_ecc_damp();
std::vector<bool> test_incl_damp();
std::vector<bool> test_gas_drag();
std::vector<bool> test_regime_switching();

static void setup_particles() {
    constexpr double Rs_per_Msun = 4.24e-6;
    M     = 1e8  * 1_Ms;
    m_sbh = 10.0 * 1_Ms;
    central_smbh = Particle{M,     M     * Rs_per_Msun * 1_Rs};
    sbh          = Particle{m_sbh, m_sbh * Rs_per_Msun * 1_Rs};
    
}

static void write_mock_csv(const std::string &path) {
/* Needed because AGNDisk will terminate if no disk file is initialized*/

    std::ofstream f(path);
    f << "R,R_Rg,Tc,rho,P,cs,H,visc,Sigma,Q,grad_T,grad_Sigma,grad_P,gamma,f_thermal\n";
    auto row = [&](double R, double H, double rho, double Tc, double cs, double Sigma) {
        f << R << "," << R*10 << "," << Tc << "," << rho << ","
          << rho*cs*cs << "," << cs << "," << H << ",0.01,"
          << Sigma << ",1.0,-0.5,-1.5,-2.5,1.66,1.0\n";
    };
    row(1.0, 0.05, 1e-8, 300.0, 1.0,  100.0);
    row(2.0, 0.10, 5e-9, 150.0, 0.5,   50.0);
    row(3.0, 0.15, 1e-9,  75.0, 0.25,  25.0);
}

int main() {
    // Leftover example exploration of add_acc_to — kept for reference.
    // write_mock_csv("/tmp/agn_disk_test.csv");
    // AGNDisk::init_from_file("/tmp/agn_disk_test.csv");
    // setup_particles();
    // auto orb = orbit::Elliptic(M, m_sbh, 2.0_AU, 0.0, 0.0_deg, 0.0_deg, 0.0_deg, 0.0_deg);
    //
    // orbit::move_particles(orb, sbh);
    // orbit::move_to_COM_frame(central_smbh, sbh);
    //
    // Particles ptcs{0.0, orbit::group(central_smbh, sbh)};
    // VectorArray acc(ptcs.number(), Vector{0, 0, 0});
    // AGNDisk::add_acc_to(ptcs, acc);

    

    struct NamedTest { const char *name; std::vector<bool> results; };
    std::vector<NamedTest> tests = {
        {"test_load_disk_data",  test_load_disk_data()},
        {"test_ecc_damp",        test_ecc_damp()},
        {"test_incl_damp",       test_incl_damp()},
        {"test_gas_drag",        test_gas_drag()},
        {"test_regime_switching", test_regime_switching()},
    };

    int total = 0, passed = 0;
    for (auto const &t : tests) {
        for (size_t i = 0; i < t.results.size(); ++i) {
            ++total;
            if (t.results[i]) ++passed;
            else std::cout << t.name << " check " << i << " FAILED\n";
        }
    }
    std::cout << "agn-disk-test: " << passed << "/" << total << " passed\n";
    return 0;
}


std::vector<bool> test_load_disk_data() {
    // !c pass a few different pre-generated csvs to load_disk_data and make sure the data loaded in matches hard-coded expectations
    // load_disk_data does no unit conversion: each loaded field must equal the raw CSV value exactly.

    // CSV A: every column present (incl. optional gamma, f_thermal, v_disk), two rows.
    std::ofstream("/tmp/agn_load_A.csv")
        << "R,R_Rg,Tc,rho,P,cs,H,visc,Sigma,Q,grad_T,grad_Sigma,grad_P,gamma,f_thermal,v_disk\n"
        << "1.0,10.0,300.0,1e-8,2.0,1.5,0.05,0.01,100.0,1.0,-0.5,-1.5,-2.5,1.4,0.9,0.7\n"
        << "2.0,20.0,150.0,5e-9,1.0,0.5,0.10,0.02,50.0,0.8,-0.6,-1.6,-2.6,1.66,1.0,0.3\n";

    AGNDisk::load_disk_data("/tmp/agn_load_A.csv");
    auto a = AGNDisk::disk_table;  // copy: the reload below overwrites the shared static table

    // CSV B: optional columns omitted (stops at grad_P), one row — exercises defaults + reload clears the table.
    std::ofstream("/tmp/agn_load_B.csv")
        << "R,R_Rg,Tc,rho,P,cs,H,visc,Sigma,Q,grad_T,grad_Sigma,grad_P\n"
        << "3.0,30.0,75.0,1e-9,0.5,0.25,0.15,0.03,25.0,0.6,-0.7,-1.7,-2.7\n";

    AGNDisk::load_disk_data("/tmp/agn_load_B.csv");
    const auto& b = AGNDisk::disk_table;

    if (verbose && false) {
        std::cout << "load_disk_data [CSV A]:\n";
        std::cout << "  size: " << a.size() << " expected 2\n";
        std::cout << "  a[0].R=" << a[0].R << " [AU] expected 1.0\n";
        std::cout << "  a[1].R=" << a[1].R << " [AU] expected 2.0\n";
        std::cout << "  a[0].R_Rg=" << a[0].R_Rg << " [dimensionless] expected 10.0\n";
        std::cout << "  a[0].Tc=" << a[0].Tc << " [K] expected 300.0\n";
        std::cout << "  a[0].rho=" << a[0].rho << " [Ms/AU^3] expected 1e-8\n";
        std::cout << "  a[0].P=" << a[0].P << " [Ms/(AU*(yr/2pi)^2)] expected 2.0\n";
        std::cout << "  a[0].cs=" << a[0].cs << " [AU/(yr/2pi)] expected 1.5\n";
        std::cout << "  a[0].H=" << a[0].H << " [AU] expected 0.05\n";
        std::cout << "  a[0].visc=" << a[0].visc << " [dimensionless] expected 0.01\n";
        std::cout << "  a[0].Sigma=" << a[0].Sigma << " [Ms/AU^2] expected 100.0\n";
        std::cout << "  a[0].Q=" << a[0].Q << " [dimensionless] expected 1.0\n";
        std::cout << "  a[0].grad_T=" << a[0].grad_T << " [dimensionless] expected -0.5\n";
        std::cout << "  a[0].grad_Sigma=" << a[0].grad_Sigma << " [dimensionless] expected -1.5\n";
        std::cout << "  a[0].grad_P=" << a[0].grad_P << " [dimensionless] expected -2.5\n";
        std::cout << "  a[0].gamma=" << a[0].gamma << " [dimensionless] expected 1.4\n";
        std::cout << "  a[0].f_thermal=" << a[0].f_thermal << " [dimensionless] expected 0.9\n";
        std::cout << "  a[0].v_disk=" << a[0].v_disk << " [AU/(yr/2pi)] expected 0.7\n";
        std::cout << "load_disk_data [CSV B]:\n";
        std::cout << "  size: " << b.size() << " expected 1\n";
        std::cout << "  b[0].R=" << b[0].R << " [AU] expected 3.0\n";
        std::cout << "  b[0].grad_P=" << b[0].grad_P << " [dimensionless] expected -2.7\n";
        std::cout << "  b[0].gamma=" << b[0].gamma << " [dimensionless] expected " << 5.0/3.0 << " (default)\n";
        std::cout << "  b[0].f_thermal=" << b[0].f_thermal << " [dimensionless] expected 1.0 (default)\n";
        std::cout << "  b[0].v_disk=" << b[0].v_disk << " [AU/(yr/2pi)] expected 0.0 (default)\n";
    }
    return {
        a.size() == 2,                                  // both rows loaded
        a[0].R == 1.0,        a[1].R == 2.0,
        a[0].R_Rg == 10.0,    a[0].Tc == 300.0,
        a[0].rho == 1e-8,     a[0].P == 2.0,
        a[0].cs == 1.5,       a[0].H == 0.05,
        a[0].visc == 0.01,    a[0].Sigma == 100.0,
        a[0].Q == 1.0,        a[0].grad_T == -0.5,
        a[0].grad_Sigma == -1.5, a[0].grad_P == -2.5,
        a[0].gamma == 1.4,    a[0].f_thermal == 0.9,
        a[0].v_disk == 0.7,                             // optional columns read when present

        b.size() == 1,                                  // reload cleared the previous table
        b[0].R == 3.0,        b[0].grad_P == -2.7,
        b[0].gamma == 5.0/3.0,                          // defaults applied when optional columns absent
        b[0].f_thermal == 1.0,
        b[0].v_disk == 0.0,
    };
}

std::vector<bool> test_init_from_file() {
    //todo
    return {true};
}

std::vector<bool> test_interp_all() {
    //make sure interpolated values are reasonably close to some manually computed examples
    return {true};
}

std::vector<bool> test_disk_v() {
    // feed known values and compare to known answer
    return {true};
}


/*Force tests:
Null-effect case: case wherein the acceleration contribution of the force being tested should be zero

*/
std::vector<bool> test_ecc_damp() {
    // Eccentricity-damping acceleration (CN08):  a_e = -2 (v.r)/(r t_e) r_hat.
    // SpaceHub natural units: length [AU], mass [Ms], time [yr/2pi], speed [AU/(yr/2pi)].
    // e_tilde = e/(H/r) and i_h = i/(H/r) are dimensionless.

    // null-effect 1: a circular, in-plane orbit. For a circle the velocity is perpendicular
    // to the radius, so v.r = 0 and the damping vanishes exactly (here e_tilde = i_h = 0).
    Vector r_hat{1.0, 0.0, 0.0};          // unit radial direction (dimensionless)
    auto a1 = AGNDisk::accel_ecc_damp(/*e_tilde=*/0.0, /*i_h=*/0.0,
                                      /*t_wave [yr/2pi]=*/1.0,   // irrelevant here: v.r = 0
                                      /*vdotr [AU^2/(yr/2pi)]=*/0.0,
                                      /*r_mag=*/2_AU, r_hat);
    bool null1 = (a1.norm() == 0.0);

    // null-effect 2: a near-empty disk. The wave time t_wave is inversely proportional to the
    // surface density Sigma, so as Sigma -> 0 the damping time t_e -> infinity and the damping
    // vanishes -- even for an eccentric orbit with v.r != 0. We use Sigma = 1e-16 (not 0) to
    // avoid the 1/Sigma divide-by-zero; the remaining quantities are physical values at this
    // radius and only set the overall scale of t_wave, not the (zero) null result.
    double Sigma = 1e-16;                            // surface density [Ms/AU^2], ~ zero
    double m0 = 1e8_Ms, mi = 10_Ms;                  // SMBH + stellar-mass BH
    double a_orb = 2_AU;                             // semi-major axis
    double aspect_ratio = 0.02;                      // H/r (dimensionless)
    double Omega_k = std::sqrt(consts::G * m0 / (a_orb * a_orb * a_orb));  // Keplerian freq [(yr/2pi)^-1]
    // wave time (CN08 / agn-disk.hpp): t_wave = (M*/m)(M*/(Sigma a^2))(H/r)^4 / Omega_k  [yr/2pi]
    double t_wave = (m0 / mi) * (m0 / Sigma / (a_orb * a_orb)) * std::pow(aspect_ratio, 4) / Omega_k;
    auto a2 = AGNDisk::accel_ecc_damp(/*e_tilde=*/1.0, /*i_h=*/0.0, t_wave,
                                      /*vdotr [AU^2/(yr/2pi)]=*/0.5,   // nonzero: eccentric orbit
                                      /*r_mag=*/2_AU, r_hat);
    bool null2 = (a2.norm() < 1e-15);     // damping is negligible

    if (verbose) {
        std::cout << "ecc_damp null1: |a|=" << a1.norm() << " [AU/(yr/2pi)^2] expected 0.0\n";
        std::cout << "ecc_damp null2: |a|=" << a2.norm() << " [AU/(yr/2pi)^2] expected < 1e-15\n";
    }
    return {null1, null2};
}

std::vector<bool> test_migration() {
    bool null1 = true; // null-effect 1: zero local surface density (must be 1e-16 or something to avoid divide by zero errors)
    bool null2 = true; // null-effect 2: Special case - e~2H/R (migration reversal threshold)
    return {null1, null2};
}

std::vector<bool> test_incl_damp() {
    // Inclination-damping acceleration (CN08):  a_i = -(v_z / t_i) z_hat.
    // e_tilde = e/(H/r), i_h = i/(H/r) are dimensionless; t_wave is a time [yr/2pi]; the
    // vertical velocity dvz is a speed [AU/(yr/2pi)]. The null result depends only on dvz,
    // so the t_wave and e_tilde values below are arbitrary placeholders.

    // null-effect 1: a non-inclined (midplane) orbit. i_h = 0 and the vertical velocity
    // dvz = 0, so a_i = 0 exactly.
    auto a1 = AGNDisk::accel_inc_damp<Vector>(/*e_tilde=*/0.5, /*i_h=*/0.0,
                                              /*t_wave [yr/2pi]=*/1.0, /*dvz [AU/(yr/2pi)]=*/0.0);
    bool null1 = (a1.norm() == 0.0);

    // null-effect 2: the anti-nodes of an inclined orbit. The orbit is genuinely inclined
    // (i_h != 0), but at maximum vertical excursion the vertical velocity is instantaneously
    // zero (dvz = 0), so a_i = 0 there as well.
    auto a2 = AGNDisk::accel_inc_damp<Vector>(/*e_tilde=*/0.5, /*i_h=*/2.0,
                                              /*t_wave [yr/2pi]=*/1.0, /*dvz [AU/(yr/2pi)]=*/0.0);
    bool null2 = (a2.norm() == 0.0);

    if (verbose) {
        std::cout << "incl_damp null1: |a|=" << a1.norm() << " [AU/(yr/2pi)^2] expected 0.0\n";
        std::cout << "incl_damp null2: |a|=" << a2.norm() << " [AU/(yr/2pi)^2] expected 0.0\n";
    }
    return {null1, null2};
}

std::vector<bool> test_gas_drag() {
    // Gas-drag acceleration is proportional to the gas density rho and to the relative
    // velocity v_rel, so either being zero removes the drag entirely.
    // Natural units: rho [Ms/AU^3], speed [AU/(yr/2pi)], mass [Ms], length [AU];
    // I_factor (dynamical-friction form factor) and Mach are dimensionless.

    double cs    = 15_kms;                // disk sound speed (~15 km/s for an AGN disk)
    double m_i   = 10_Ms;                 // satellite (stellar-mass BH)
    double r_eff = 1e-4_AU;               // effective (accretion) radius
    double I_factor = 0.1;                // dynamical-friction form factor (dimensionless)

    // null-effect 1: zero local gas density (rho = 0). The other inputs are physical and
    // self-consistent (vmag, v2, Mach derived from v_rel and cs); with no gas there is
    // nothing to drag against.
    Vector v_rel{100_kms, 0.0, 50_kms};   // body's velocity relative to the gas (~112 km/s)
    double vmag = v_rel.norm();
    double v2   = vmag * vmag;
    double Mach = vmag / cs;
    auto a1 = AGNDisk::accel_gas_drag(/*rho=*/0.0, cs, m_i, r_eff, vmag, v2, I_factor, Mach, v_rel);
    bool null1 = (a1.norm() == 0.0);

    // null-effect 2: zero relative velocity (the body moves with the gas). rho is nonzero
    // here (1e-8 Ms/AU^3) to show it is the lack of relative motion, not the gas, that
    // removes the drag.
    Vector v_rel_zero{0.0, 0.0, 0.0};
    auto a2 = AGNDisk::accel_gas_drag(/*rho=*/1e-8, cs, m_i, r_eff,
                                      /*vmag=*/0.0, /*v2=*/0.0, I_factor, /*Mach=*/0.0, v_rel_zero);
    bool null2 = (a2.norm() == 0.0);

    if (verbose) {
        std::cout << "gas_drag null1: |a|=" << a1.norm() << " [AU/(yr/2pi)^2] expected 0.0\n";
        std::cout << "gas_drag null2: |a|=" << a2.norm() << " [AU/(yr/2pi)^2] expected 0.0\n";
    }
    return {null1, null2};
}

std::vector<bool> test_regime_switching() {
    /*
            Case Description            |       Expected result
    I: Fully retrograde orbit (i=180)   |   No migration/damping, gas drag forces active
    II: i=0, e=0, a ~ 200Rg             |   No gas drag, migration active, i/e damping does not cause any issues
    III: i~H/R(~0.1deg), e=0, a ~ 200Rg |   inclination damping occurs, migration occurs, no gas drag, ecc damping no issues caused
    IV: i ~ 70, e=0, a ~ 200Rg          |   only gas drag
    V: i=0, e>0.3, a~200Rg              |   Only gas drag

    The regime decision lives in the local bool in_typeI_regime inside add_acc_to(); it is the
    public predicate AGNDisk::is_typeI_regime() that add_acc_to() actually calls. Part A checks the
    boolean is SET correctly; Part B checks it SELECTS the right force set end-to-end.
    */
    std::vector<bool> out;
    auto push = [&](bool b) { out.push_back(b); };

    // ============================================================================================
    // Part A: in_typeI_regime is SET correctly -- exercise the predicate directly.
    // aspect_ratio (H/R) = 0.05  ->  default i_max = 1.5*h = 0.075 rad, e_max = 4*h = 0.2.
    // ============================================================================================
    const double h     = 0.05;
    const double i_max = 1.5 * h;           // 0.075 rad (~4.3 deg)
    const double e_max = 4.0 * h;           // 0.2
    const double deg   = consts::pi / 180.0;

    push(AGNDisk::is_typeI_regime(consts::pi, 0.0, h) == false);  // I   retrograde -> drag
    push(AGNDisk::is_typeI_regime(0.0,        0.0, h) == true);   // II  embedded   -> type I
    push(AGNDisk::is_typeI_regime(0.02,       0.0, h) == true);   // III low incl   -> type I
    push(AGNDisk::is_typeI_regime(70.0*deg,   0.0, h) == false);  // IV  i=70       -> drag
    push(AGNDisk::is_typeI_regime(0.0,        0.3, h) == false);  // V   e>e_max    -> drag

    // a retrograde orbit is NEVER in the type I regime, at any inclination > 90 deg
    push(AGNDisk::is_typeI_regime(91.0*deg,  0.0, h) == false);
    push(AGNDisk::is_typeI_regime(135.0*deg, 0.0, h) == false);

    // boundaries straddling i_max and e_max
    push(AGNDisk::is_typeI_regime(0.99*i_max, 0.0,        h) == true);
    push(AGNDisk::is_typeI_regime(1.01*i_max, 0.0,        h) == false);
    push(AGNDisk::is_typeI_regime(0.0,        0.99*e_max, h) == true);
    push(AGNDisk::is_typeI_regime(0.0,        1.01*e_max, h) == false);

    // static overrides take precedence over the default h-based thresholds (reset after each)
    AGNDisk::type_i_inclination_max = 10.0_deg;                         // 0.175 rad > default 0.075
    bool inc_override = (AGNDisk::is_typeI_regime( 8.0*deg, 0.0, h) == true)    // inside override
                     && (AGNDisk::is_typeI_regime(12.0*deg, 0.0, h) == false);  // outside override
    AGNDisk::type_i_inclination_max = 0.0_deg;
    push(inc_override);

    AGNDisk::type_i_eccentricity_max = 0.5;                             // > default 0.2
    bool ecc_override = (AGNDisk::is_typeI_regime(0.0, 0.3, h) == true)        // inside override
                     && (AGNDisk::is_typeI_regime(0.0, 0.6, h) == false);      // outside override
    AGNDisk::type_i_eccentricity_max = 0.0;
    push(ecc_override);

    // ============================================================================================
    // Part B: the boolean SELECTS the right force set -- checked end-to-end through add_acc_to.
    // The mock disk has no v_disk column -> v_disk = 0 -> v_rel = v[1] - v[0]. Pure gas drag is
    // exactly anti-parallel to v_rel, whereas the type I force-set is not collinear with v_rel.
    // ============================================================================================
    AGNDisk::mute_diagnostics = true;
    write_mock_csv("/tmp/agn_regime.csv");                   // H/R = 0.05 at R = 1, 2, 3 AU
    AGNDisk::init_from_file("/tmp/agn_regime.csv");
    setup_particles();

    const double a       = 2.0;     // semi-major axis [AU], inside [Rmin, Rmax] = [1, 3]
    const double rel_tol = 1e-9;    // "collinear / perpendicular / zero" relative tolerance
    const double sig_tol = 1e-6;    // "clearly non-collinear" threshold

    // build one orbit and return {acc on the secondary, v_rel, dr}
    auto run_case = [&](double e, double incl, double nu) {
        Particle c = central_smbh, s = sbh;
        auto orb = orbit::Elliptic(M, m_sbh, a, e, incl, 0.0, 0.0, nu);
        orbit::move_particles(orb, s);
        orbit::move_to_COM_frame(c, s);
        Particles ptcs{0.0, orbit::group(c, s)};
        VectorArray acc(ptcs.number(), Vector{0, 0, 0});
        AGNDisk::add_acc_to(ptcs, acc);
        Vector vrel = ptcs.vel()[1] - ptcs.vel()[0];
        Vector dr   = ptcs.pos()[1] - ptcs.pos()[0];
        return std::make_tuple(acc[1], vrel, dr);
    };

    // pure gas drag:  a = -coef * v_rel  ->  collinear with, and opposing, v_rel
    auto is_drag = [&](const Vector &a1, const Vector &vrel) {
        bool collinear = cross(a1, vrel).norm() <= rel_tol * a1.norm() * vrel.norm();
        bool opposing  = dot(a1, vrel) < 0.0;
        return a1.norm() > 0.0 && collinear && opposing;
    };
    // type I force-set: not collinear with v_rel
    auto not_drag = [&](const Vector &a1, const Vector &vrel) {
        return a1.norm() > 0.0 && cross(a1, vrel).norm() > sig_tol * a1.norm() * vrel.norm();
    };

    // Case I: retrograde (i=180), circular -> gas drag only.
    auto [aI, vI, drI] = run_case(0.0, consts::pi, 0.0);
    push(is_drag(aI, vI));

    // Case II: i=0, e=0 -> type I; e/i damping vanish (v.r = 0, v_z = 0) so only the migration
    // torque survives. A torque is in-plane and perpendicular to the radius. Geometry alone cannot
    // separate migration from drag here (both tangential) -> Part A pins this to the type I branch.
    auto [aII, vII, drII] = run_case(0.0, 0.0, 0.0);
    bool II_in_plane   = std::abs(aII.z) <= rel_tol * aII.norm();
    bool II_tangential = std::abs(dot(aII, drII)) <= rel_tol * aII.norm() * drII.norm();
    push(aII.norm() > 0.0 && II_in_plane && II_tangential);

    // Case III: small inclination, e=0, at the node (v_z != 0) -> type I with active inclination
    // damping. Migration alone (circular) would be collinear with v_rel; the non-collinearity is
    // precisely the vertical inclination-damping contribution.
    auto [aIII, vIII, drIII] = run_case(0.0, 1.0*deg, 0.0);
    push(not_drag(aIII, vIII));

    // Case IV: i=70, e=0 -> gas drag only.
    auto [aIV, vIV, drIV] = run_case(0.0, 70.0*deg, 0.0);
    push(is_drag(aIV, vIV));

    // Case V: i=0, e=0.3 (> e_max) -> gas drag only. Evaluated off the apsides (nu=90) so v_rel
    // has both radial and tangential parts (relied on by the flip test below).
    auto [aV, vV, drV] = run_case(0.3, 0.0, consts::pi/2);
    push(is_drag(aV, vV));

    // Regime-flip differential (the strongest bool->force proof): ONE fixed orbit (Case V); flip
    // only the regime via the eccentricity override. The boolean -- nothing else about the orbit --
    // switches the acceleration from drag (collinear with v_rel) to type I (not collinear).
    auto [aFlipOff, vFlipOff, drFlipOff] = run_case(0.3, 0.0, consts::pi/2);   // e=0.3 > e_max -> drag
    AGNDisk::type_i_eccentricity_max = 0.5;                                    // now e=0.3 < 0.5 -> type I
    auto [aFlipOn, vFlipOn, drFlipOn] = run_case(0.3, 0.0, consts::pi/2);
    AGNDisk::type_i_eccentricity_max = 0.0;                                    // reset
    push(is_drag(aFlipOff, vFlipOff) && not_drag(aFlipOn, vFlipOn));

    if (verbose) {
        std::cout << "regime_switching Part B accelerations [AU/(yr/2pi)^2]:\n";
        std::cout << "  I   (i=180, e=0)   |a|=" << aI.norm()   << "  |a x v_rel|=" << cross(aI, vI).norm()   << "  a.v_rel=" << dot(aI, vI)   << " (drag)\n";
        std::cout << "  II  (i=0,   e=0)   |a|=" << aII.norm()  << "  a.z="         << aII.z                  << "  a.r="     << dot(aII, drII) << " (migration)\n";
        std::cout << "  III (i~1d,  e=0)   |a|=" << aIII.norm() << "  |a x v_rel|=" << cross(aIII, vIII).norm() << " (type I, incl damp)\n";
        std::cout << "  IV  (i=70,  e=0)   |a|=" << aIV.norm()  << "  |a x v_rel|=" << cross(aIV, vIV).norm()  << "  a.v_rel=" << dot(aIV, vIV) << " (drag)\n";
        std::cout << "  V   (i=0, e=0.3)   |a|=" << aV.norm()   << "  |a x v_rel|=" << cross(aV, vV).norm()    << "  a.v_rel=" << dot(aV, vV)   << " (drag)\n";
        std::cout << "  flip e=0.3: drag |axv|=" << cross(aFlipOff, vFlipOff).norm()
                  << " -> typeI |axv|="          << cross(aFlipOn, vFlipOn).norm() << "\n";
    }
    return out;
}

std::vector<bool> test_acceleration_vector() {
    /*
    Cases:

    */
    return {};
}

// Compile and run:
//   g++ -std=c++17 -O0 -pthread SpaceHub/test/agn-disk-test.cpp -o SpaceHub/test/agn-disk-test
//   ./SpaceHub/test/agn-disk-test
