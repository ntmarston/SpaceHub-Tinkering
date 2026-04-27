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
#include <fstream>
#include <string>

#include "../../src/type-class.hpp"
#include "../../src/particles/finite-size.hpp"
#include "../../src/orbits/orbits.hpp"
#include "../../src/orbits/particle-manip.hpp"
#include "../../src/interaction/disk-migration-switching.hpp"
#include "../catch.hpp"
#include "utest.hpp"

using namespace hub;
using namespace hub::unit;
using namespace hub::force;

using Type        = hub::Types<double>;
using Particles   = hub::particles::SizeParticles<Type>;
using Particle    = Particles::Particle;
using Vector      = Particles::Vector;
using VectorArray = Type::VectorArray;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static void write_mock_csv(const std::string &path, bool include_optional_cols) {
    std::ofstream f(path);
    if (include_optional_cols)
        f << "R,R_Rg,Tc,rho,P,cs,H,visc,Sigma,Q,grad_T,grad_Sigma,grad_P,gamma,f_thermal\n";
    else
        f << "R,R_Rg,Tc,rho,P,cs,H,visc,Sigma,Q,grad_T,grad_Sigma,grad_P\n";

    // H/R = 0.05 throughout (H = 0.05*R)
    auto row = [&](double R, double H, double rho, double Tc, double cs, double Sigma) {
        f << R << "," << R*10 << "," << Tc << "," << rho << ","
          << rho*cs*cs << "," << cs << "," << H << ",0.01,"
          << Sigma << ",1.0,-0.5,-1.5,-2.5";
        if (include_optional_cols) f << ",1.66,1.0";
        f << "\n";
    };
    row(1.0, 0.05,  1e-8,  300.0, 1.0,  100.0);
    row(2.0, 0.10,  5e-9,  150.0, 0.5,   50.0);
    row(3.0, 0.15,  1e-9,   75.0, 0.25,  25.0);
}

// Reset all DISABLE_ overrides to false (all forces active, governed by switching logic).
static void reset_flags() {
    DiskMigration::DISABLE_MIGRATION          = false;
    DiskMigration::DISABLE_E_DAMPING          = false;
    DiskMigration::DISABLE_I_DAMPING          = false;
    DiskMigration::DISABLE_DYNAMICAL_FRICTION = false;
    DiskMigration::DISABLE_AERODYNAMIC_DRAG   = false;
    DiskMigration::DISABLE_BONDI_HOYLE        = false;
    DiskMigration::ignore_dynamical_friction  = false;
}

// Build a 2-body system (SMBH + stellar-mass BH) on a Keplerian orbit.
// sma in AU; default sma=2 places body at R=2 (disk midpoint in mock CSV).
// TA=0 at construction, so the initial position is at periapsis.
static Particles make_particles(double ecc, double incl_deg, double sma = 2.0) {
    double M = 1e8, m = 10.0;  // SMBH + stellar-mass BH (M☉)
    double r_star = 1e-6, r_body = 1e-6;
    Particle star{M, r_star};
    Particle body{m, r_body};

    auto orb = orbit::Elliptic(M, m, sma, ecc, incl_deg * consts::pi / 180.0, 0.0, 0.0, 0.0);
    orbit::move_particles(orb, body);
    orbit::move_to_COM_frame(star, body);
    return Particles(0.0, orbit::group(star, body));
}

// Call add_acc_to and return the resulting acceleration array.
static VectorArray compute_acc(const Particles &ptcs) {
    VectorArray acc(ptcs.number(), Vector{0, 0, 0});
    DiskMigration::add_acc_to(ptcs, acc);
    return acc;
}

// Check whether body[1] received a nonzero net acceleration.
static bool acc_nonzero(const VectorArray &acc) {
    auto &a = acc[1];  // acc[0] is the central star; disk forces act on body at index 1
    return (a.x != 0.0 || a.y != 0.0 || a.z != 0.0);
}

// Compute OrbitalRegime for body[1] from its current particle state and the
// loaded disk table.  Mirrors the orbital-element extraction in add_acc_to.
static DiskMigration::OrbitalRegime get_regime(const Particles &ptcs) {
    auto const &p = ptcs.pos();
    auto const &v = ptcs.vel();
    auto const &m = ptcs.mass();
    auto dr      = p[1] - p[0];
    auto dv      = v[1] - v[0];
    double R_cyl = sqrt(dr.x*dr.x + dr.y*dr.y);
    double u     = consts::G * (m[0] + m[1]);
    auto [a_orb, ecc] = orbit::calc_a_e(u, dr, dv);
    auto h_vec   = cross(dr, dv);
    double h_mag = norm(h_vec);
    double incl  = acos(std::clamp(h_vec.z / h_mag, -1.0, 1.0));
    double sin_i = sqrt(h_vec.x*h_vec.x + h_vec.y*h_vec.y) / h_mag;
    auto props   = DiskMigration::interp_all(R_cyl);
    double aspect_ratio = props.H / R_cyl;
    return DiskMigration::classify(ecc, incl, sin_i, aspect_ratio);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("DiskMigration", "[DiskMigration]") {

    const std::string mock_csv      = "/tmp/test_disk_profile.csv";
    const std::string mock_csv_bare = "/tmp/test_disk_profile_bare.csv";

    // ------------------------------------------------------------------
    SECTION("load_disk_data") {
        write_mock_csv(mock_csv, true);
        DiskMigration::init_from_file(mock_csv);

        REQUIRE(DiskMigration::disk_table.size() == 3);
        REQUIRE(DiskMigration::disk_table[1].Sigma == APPROX(50.0));
        REQUIRE(DiskMigration::Rmin == APPROX(1.0));
        REQUIRE(DiskMigration::Rmax == APPROX(3.0));
    }

    // ------------------------------------------------------------------
    SECTION("load_disk_data defaults") {
        write_mock_csv(mock_csv_bare, false);
        DiskMigration::init_from_file(mock_csv_bare);

        REQUIRE(DiskMigration::disk_table[0].gamma     == APPROX(5.0 / 3.0));
        REQUIRE(DiskMigration::disk_table[0].f_thermal == APPROX(1.0));
    }

    // ------------------------------------------------------------------
    SECTION("interp_all at node") {
        write_mock_csv(mock_csv, true);
        DiskMigration::init_from_file(mock_csv);

        auto props = DiskMigration::interp_all(2.0);
        REQUIRE(props.Sigma == APPROX(50.0));
        REQUIRE(props.H     == APPROX(0.10));
        REQUIRE(props.cs    == APPROX(0.5));
    }

    // ------------------------------------------------------------------
    SECTION("interp_all midpoint") {
        write_mock_csv(mock_csv, true);
        DiskMigration::init_from_file(mock_csv);

        auto props = DiskMigration::interp_all(1.5);
        REQUIRE(props.Sigma > 50.0);
        REQUIRE(props.Sigma < 100.0);
    }

    // ------------------------------------------------------------------
    SECTION("interp_all out of bounds") {
        write_mock_csv(mock_csv, true);
        DiskMigration::init_from_file(mock_csv);

        auto lo = DiskMigration::interp_all(0.5);
        REQUIRE(lo.Sigma     == APPROX(0.0));
        REQUIRE(lo.gamma     == APPROX(5.0 / 3.0));
        REQUIRE(lo.f_thermal == APPROX(1.0));

        auto hi = DiskMigration::interp_all(3.5);
        REQUIRE(hi.Sigma == APPROX(0.0));
    }

    // ------------------------------------------------------------------
    // T3 switching tests.
    // All use the 3-row mock CSV (H/R = 0.05 throughout).
    // Default body position: sma=2.0 AU → R_cyl≈2.0 AU, H/R≈0.05.
    //
    // Each section verifies two things:
    //   1. classify() returns the expected OrbitalRegime flags
    //   2. add_acc_to modifies (or leaves zero) the acceleration vector correctly
    // ------------------------------------------------------------------

    SECTION("switching: over-inclined") {
        // i=10°, e=0 → sin(10°)≈0.174 > H/R=0.05 → not embedded → dyn_fric ON
        write_mock_csv(mock_csv, true);
        DiskMigration::init_from_file(mock_csv);
        auto ptcs = make_particles(0.0, 10.0);

        auto regime = get_regime(ptcs);
        REQUIRE(regime.retrograde         == false);
        REQUIRE(regime.embedded           == false);
        REQUIRE(regime.dynamical_friction == true);

        // Drag fires when migration/damping overridden off
        reset_flags();
        DiskMigration::DISABLE_MIGRATION = true;
        DiskMigration::DISABLE_E_DAMPING = true;
        DiskMigration::DISABLE_I_DAMPING = true;
        REQUIRE(acc_nonzero(compute_acc(ptcs)));

        // Migration/damping produce no acc: switching suppresses them (not embedded)
        reset_flags();
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        DiskMigration::DISABLE_E_DAMPING          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        REQUIRE_FALSE(acc_nonzero(compute_acc(ptcs)));
    }

    SECTION("switching: over-eccentric") {
        // e > DiskMigration::e_max, i=0 → embedded but over-eccentric → dyn_fric ON
        // sma=7: periapsis = sma*(1-e) = 1.4 AU, inside [1,3]; sin(0)=0 < H/R
        write_mock_csv(mock_csv, true);
        DiskMigration::init_from_file(mock_csv);
        auto ptcs = make_particles(0.8, 0.0, 7.0);

        auto regime = get_regime(ptcs);
        REQUIRE(regime.retrograde         == false);
        REQUIRE(regime.embedded           == true);
        REQUIRE(regime.dynamical_friction == true);

        // Drag fires
        reset_flags();
        DiskMigration::DISABLE_MIGRATION = true;
        DiskMigration::DISABLE_E_DAMPING = true;
        DiskMigration::DISABLE_I_DAMPING = true;
        REQUIRE(acc_nonzero(compute_acc(ptcs)));

        // Migration/damping produce no acc: switching suppresses them (over-eccentric)
        reset_flags();
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        DiskMigration::DISABLE_E_DAMPING          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        REQUIRE_FALSE(acc_nonzero(compute_acc(ptcs)));
    }

    SECTION("switching: retrograde") {
        // i=180°, e=0 → retrograde → dyn_fric ON immediately
        write_mock_csv(mock_csv, true);
        DiskMigration::init_from_file(mock_csv);
        auto ptcs = make_particles(0.0, 180.0);

        auto regime = get_regime(ptcs);
        REQUIRE(regime.retrograde         == true);
        REQUIRE(regime.dynamical_friction == true);

        // Drag fires
        reset_flags();
        DiskMigration::DISABLE_MIGRATION = true;
        DiskMigration::DISABLE_E_DAMPING = true;
        DiskMigration::DISABLE_I_DAMPING = true;
        REQUIRE(acc_nonzero(compute_acc(ptcs)));

        // Migration/damping produce no acc: switching suppresses them (retrograde)
        reset_flags();
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        DiskMigration::DISABLE_E_DAMPING          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        REQUIRE_FALSE(acc_nonzero(compute_acc(ptcs)));
    }

    SECTION("switching: embedded circular coplanar") {
        // i=0, e=0 → regime: embedded, in_plane, dyn_fric OFF
        // Migration ON; e-damping acc is zero (T_bar ∝ e² → 0); i-damping blocked (in_plane)
        write_mock_csv(mock_csv, true);
        DiskMigration::init_from_file(mock_csv);
        auto ptcs = make_particles(0.0, 0.0);

        auto regime = get_regime(ptcs);
        REQUIRE(regime.retrograde         == false);
        REQUIRE(regime.embedded           == true);
        REQUIRE(regime.in_plane           == true);
        REQUIRE(regime.dynamical_friction == false);

        // Migration produces nonzero acc
        reset_flags();
        DiskMigration::DISABLE_E_DAMPING          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        REQUIRE(acc_nonzero(compute_acc(ptcs)));

        // e-damping alone → zero acc: formula gives T_bar ∝ e² = 0 for circular orbit
        reset_flags();
        DiskMigration::DISABLE_MIGRATION          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        REQUIRE_FALSE(acc_nonzero(compute_acc(ptcs)));

        // All forces off → zero acc
        reset_flags();
        DiskMigration::DISABLE_MIGRATION          = true;
        DiskMigration::DISABLE_E_DAMPING          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        REQUIRE_FALSE(acc_nonzero(compute_acc(ptcs)));
    }

    SECTION("switching: embedded eccentric coplanar") {
        // i=0, e=0.02 → embedded, in_plane, dyn_fric OFF → migration + e-damping ON
        write_mock_csv(mock_csv, true);
        DiskMigration::init_from_file(mock_csv);
        auto ptcs = make_particles(0.02, 0.0);

        auto regime = get_regime(ptcs);
        REQUIRE(regime.retrograde         == false);
        REQUIRE(regime.embedded           == true);
        REQUIRE(regime.in_plane           == true);
        REQUIRE(regime.dynamical_friction == false);

        // e-damping alone → nonzero acc (e > 0)
        reset_flags();
        DiskMigration::DISABLE_MIGRATION          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        REQUIRE(acc_nonzero(compute_acc(ptcs)));

        // Migration alone → nonzero acc
        reset_flags();
        DiskMigration::DISABLE_E_DAMPING          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        REQUIRE(acc_nonzero(compute_acc(ptcs)));
    }

    SECTION("switching: embedded circular inclined") {
        // i=1°, e=0 → embedded, not in_plane, dyn_fric OFF → migration + i-damping ON
        write_mock_csv(mock_csv, true);
        DiskMigration::init_from_file(mock_csv);
        auto ptcs = make_particles(0.0, 1.0);

        auto regime = get_regime(ptcs);
        REQUIRE(regime.retrograde         == false);
        REQUIRE(regime.embedded           == true);
        REQUIRE(regime.in_plane           == false);
        REQUIRE(regime.dynamical_friction == false);

        // i-damping alone → nonzero acc
        reset_flags();
        DiskMigration::DISABLE_MIGRATION          = true;
        DiskMigration::DISABLE_E_DAMPING          = true;
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        REQUIRE(acc_nonzero(compute_acc(ptcs)));

        // Migration alone → nonzero acc
        reset_flags();
        DiskMigration::DISABLE_E_DAMPING          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        REQUIRE(acc_nonzero(compute_acc(ptcs)));
    }

    SECTION("switching: embedded eccentric inclined") {
        // i=1°, e=0.02 → embedded, not in_plane, dyn_fric OFF → migration + e + i-damping ON
        write_mock_csv(mock_csv, true);
        DiskMigration::init_from_file(mock_csv);
        auto ptcs = make_particles(0.02, 1.0);

        auto regime = get_regime(ptcs);
        REQUIRE(regime.retrograde         == false);
        REQUIRE(regime.embedded           == true);
        REQUIRE(regime.in_plane           == false);
        REQUIRE(regime.dynamical_friction == false);

        // e-damping alone → nonzero acc
        reset_flags();
        DiskMigration::DISABLE_MIGRATION          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        REQUIRE(acc_nonzero(compute_acc(ptcs)));

        // i-damping alone → nonzero acc
        reset_flags();
        DiskMigration::DISABLE_MIGRATION          = true;
        DiskMigration::DISABLE_E_DAMPING          = true;
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        REQUIRE(acc_nonzero(compute_acc(ptcs)));

        // Migration alone → nonzero acc
        reset_flags();
        DiskMigration::DISABLE_E_DAMPING          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        REQUIRE(acc_nonzero(compute_acc(ptcs)));
    }

    SECTION("switching: stalled migration") {
        // i=0, e=0.06 → embedded, in_plane, dyn_fric OFF
        // e ≥ 1.0999 * H/R (0.055) → PL00 guard issues continue, stalling migration
        // Expected: e-damping ON, migration acc zero (stalled)
        write_mock_csv(mock_csv, true);
        DiskMigration::init_from_file(mock_csv);
        auto ptcs = make_particles(0.06, 0.0);

        auto regime = get_regime(ptcs);
        REQUIRE(regime.retrograde         == false);
        REQUIRE(regime.embedded           == true);
        REQUIRE(regime.in_plane           == true);
        REQUIRE(regime.dynamical_friction == false);

        // e-damping alone → nonzero acc (embedded, e > 0)
        reset_flags();
        DiskMigration::DISABLE_MIGRATION          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        REQUIRE(acc_nonzero(compute_acc(ptcs)));

        // Migration alone → zero acc (PL00 guard issues continue, skipping acc update)
        reset_flags();
        DiskMigration::DISABLE_E_DAMPING          = true;
        DiskMigration::DISABLE_I_DAMPING          = true;
        DiskMigration::DISABLE_DYNAMICAL_FRICTION = true;
        DiskMigration::DISABLE_AERODYNAMIC_DRAG   = true;
        DiskMigration::DISABLE_BONDI_HOYLE        = true;
        REQUIRE_FALSE(acc_nonzero(compute_acc(ptcs)));
    }
}
