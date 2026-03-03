// Standalone interpolation method comparison test
//
// Compares the current disk-model.hpp interpolation (Catmull-Rom + Fritsch-Carlson)
// against three well-known alternatives:
//   A. Log-log linear interpolation
//   B. Clean PCHIP (Fritsch-Carlson with standard 3-point slopes)
//   C. Steffen (1990) monotone interpolation
//
// Tests: positivity, grid-point accuracy, boundedness, max/mean relative difference,
//        zone boundary behavior, near-boundary behavior, non-monotonic handling.
//
// Compile: g++ -std=c++17 -O3 -pthread test_interp_comparison.cpp -o test_interp
// Run:     ./test_interp [disk_csv_path]
//   Default: ../../SpaceHub/src/interaction/disk_tab/disk_default.csv

#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include <limits>
#include <functional>
#include <numeric>
#include <cassert>
#include <chrono>

// ============================================================================
// Minimal disk table loader (matches disk-model.hpp DiskRow)
// ============================================================================
struct DiskRow {
    double R, R_Rg, Tc, rho, P, cs, H, visc, Sigma, Q, grad_T, grad_Sigma, grad_P;
};

static std::vector<DiskRow> load_disk(const std::string& filename) {
    std::vector<DiskRow> table;
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "ERROR: Cannot open " << filename << "\n";
        std::exit(1);
    }
    std::string line;
    std::getline(file, line); // skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        DiskRow row;
        char comma;
        ss >> row.R >> comma >> row.R_Rg >> comma >> row.Tc >> comma
           >> row.rho >> comma >> row.P >> comma >> row.cs >> comma
           >> row.H >> comma >> row.visc >> comma >> row.Sigma >> comma
           >> row.Q >> comma >> row.grad_T >> comma >> row.grad_Sigma >> comma
           >> row.grad_P;
        table.push_back(row);
    }
    return table;
}

// ============================================================================
// Field accessors: pairs of (pointer-to-member, name, is_positive_definite)
// ============================================================================
struct FieldInfo {
    double DiskRow::*ptr;
    const char* name;
    bool expect_positive; // if true, we test for positivity violations
};

static const FieldInfo fields[] = {
    {&DiskRow::rho,    "rho",    true},
    {&DiskRow::cs,     "cs",     true},
    {&DiskRow::H,      "H",     true},
    {&DiskRow::grad_P, "grad_P", false},
};
static constexpr size_t N_FIELDS = 4;

// ============================================================================
// METHOD 0: Current (Catmull-Rom slopes + Fritsch-Carlson monotonicity)
// Exact copy from disk-model.hpp interp()
// ============================================================================
namespace current {

static double interp(const std::vector<DiskRow>& table, double R, double DiskRow::*field) {
    auto it = std::lower_bound(table.begin(), table.end(), R,
        [](const DiskRow& row, double r) { return row.R < r; });

    size_t i = std::clamp<size_t>(it - table.begin(), 1, table.size() - 2);
    size_t i0 = (i > 1) ? i - 1 : 0;
    size_t i1 = i;
    size_t i2 = i + 1;
    size_t i3 = std::min(i + 2, table.size() - 1);

    double x0 = table[i0].R, x1 = table[i1].R, x2 = table[i2].R, x3 = table[i3].R;
    double y0 = table[i0].*field, y1 = table[i1].*field,
           y2 = table[i2].*field, y3 = table[i3].*field;

    double m1 = (y2 - y0) / (x2 - x0);
    double m2 = (y3 - y1) / (x3 - x1);
    double h = x2 - x1;

    // Fritsch-Carlson monotonicity constraint
    double delta = (y2 - y1) / h;
    if (std::abs(delta) < 1e-30) {
        m1 = m2 = 0.0;
    } else {
        double alpha = m1 / delta;
        double beta  = m2 / delta;
        if (alpha <= 0.0) m1 = 0.0;
        if (beta  <= 0.0) m2 = 0.0;
        double r2 = alpha * alpha + beta * beta;
        if (r2 > 9.0) {
            double tau = 3.0 / std::sqrt(r2);
            m1 = tau * alpha * delta;
            m2 = tau * beta  * delta;
        }
    }

    double t = (R - x1) / h;
    double t2 = t * t, t3 = t2 * t;

    double h00 = 2*t3 - 3*t2 + 1;
    double h10 = t3 - 2*t2 + t;
    double h01 = -2*t3 + 3*t2;
    double h11 = t3 - t2;

    return h00*y1 + h10*h*m1 + h01*y2 + h11*h*m2;
}

} // namespace current

// ============================================================================
// Shared helper: clean binary search (matches fairbairn-datacube.hpp pattern)
// ============================================================================
static size_t find_interval(const std::vector<DiskRow>& table, double R) {
    auto it = std::upper_bound(table.begin(), table.end(), R,
        [](double r, const DiskRow& row) { return r < row.R; });
    if (it == table.begin()) return 0;
    size_t idx = static_cast<size_t>(it - table.begin()) - 1;
    return std::min(idx, table.size() - 2);
}

// ============================================================================
// METHOD A: Log-log linear interpolation
// For positive values: y = y0 * (y1/y0)^t  where t = log(R/R0)/log(R1/R0)
// For signed/zero values: falls back to linear in R-space
// ============================================================================
namespace loglog_linear {

static double interp(const std::vector<DiskRow>& table, double R, double DiskRow::*field,
                     bool positive_field) {
    size_t i = find_interval(table, R);
    double R0 = table[i].R, R1 = table[i+1].R;
    double y0 = table[i].*field, y1 = table[i+1].*field;

    if (positive_field && y0 > 0.0 && y1 > 0.0) {
        // Log-log: interpolate log(y) linearly in log(R)
        double t = std::log(R / R0) / std::log(R1 / R0);
        return y0 * std::pow(y1 / y0, t);
    } else {
        // Fallback to linear for signed or zero-valued fields
        double t = (R - R0) / (R1 - R0);
        return y0 + t * (y1 - y0);
    }
}

} // namespace loglog_linear

// ============================================================================
// METHOD B: Clean PCHIP (Fritsch-Carlson with standard 3-point slopes)
// Uses weighted harmonic mean of adjacent secants (standard PCHIP initialization),
// then applies Fritsch-Carlson monotonicity constraints.
// ============================================================================
namespace pchip {

// Precompute slopes for the entire table for a given field
static std::vector<double> compute_slopes(const std::vector<DiskRow>& table, double DiskRow::*field) {
    size_t n = table.size();
    std::vector<double> slopes(n, 0.0);

    if (n < 2) return slopes;

    // Compute secants
    std::vector<double> h(n-1), delta(n-1);
    for (size_t i = 0; i < n-1; ++i) {
        h[i] = table[i+1].R - table[i].R;
        delta[i] = (table[i+1].*field - table[i].*field) / h[i];
    }

    // Interior slopes: standard 3-point weighted formula (de Boor, Fritsch-Carlson)
    // d_i = weighted harmonic-style mean that accounts for non-uniform spacing
    for (size_t i = 1; i < n-1; ++i) {
        double w1 = 2.0 * h[i] + h[i-1];
        double w2 = h[i] + 2.0 * h[i-1];
        if (delta[i-1] * delta[i] > 0.0) {
            // Both secants same sign: weighted harmonic mean
            slopes[i] = (w1 + w2) / (w1 / delta[i-1] + w2 / delta[i]);
        } else {
            slopes[i] = 0.0;
        }
    }

    // Endpoint slopes: one-sided
    slopes[0] = delta[0];
    slopes[n-1] = delta[n-2];

    // Fritsch-Carlson monotonicity constraint on each interval
    for (size_t i = 0; i < n-1; ++i) {
        if (std::abs(delta[i]) < 1e-30) {
            slopes[i] = 0.0;
            slopes[i+1] = 0.0;
        } else {
            double alpha = slopes[i] / delta[i];
            double beta  = slopes[i+1] / delta[i];
            if (alpha <= 0.0) slopes[i] = 0.0;
            if (beta  <= 0.0) slopes[i+1] = 0.0;
            double r2 = alpha * alpha + beta * beta;
            if (r2 > 9.0) {
                double tau = 3.0 / std::sqrt(r2);
                slopes[i]   = tau * alpha * delta[i];
                slopes[i+1] = tau * beta  * delta[i];
            }
        }
    }

    return slopes;
}

static double interp(const std::vector<DiskRow>& table, const std::vector<double>& slopes,
                     double R, double DiskRow::*field) {
    size_t i = find_interval(table, R);
    double x1 = table[i].R, x2 = table[i+1].R;
    double y1 = table[i].*field, y2 = table[i+1].*field;
    double m1 = slopes[i], m2 = slopes[i+1];
    double h = x2 - x1;
    double t = (R - x1) / h;
    double t2 = t * t, t3 = t2 * t;

    double h00 = 2*t3 - 3*t2 + 1;
    double h10 = t3 - 2*t2 + t;
    double h01 = -2*t3 + 3*t2;
    double h11 = t3 - t2;

    return h00*y1 + h10*h*m1 + h01*y2 + h11*h*m2;
}

} // namespace pchip

// ============================================================================
// METHOD C: Steffen (1990) monotone interpolation
// Slopes determined by the Steffen formula — guaranteed monotone by construction
// without a separate correction pass.
// Reference: Steffen, M. (1990), A&A 239, 443-450
// ============================================================================
namespace steffen {

static std::vector<double> compute_slopes(const std::vector<DiskRow>& table, double DiskRow::*field) {
    size_t n = table.size();
    std::vector<double> slopes(n, 0.0);

    if (n < 2) return slopes;

    std::vector<double> h(n-1), s(n-1);
    for (size_t i = 0; i < n-1; ++i) {
        h[i] = table[i+1].R - table[i].R;
        s[i] = (table[i+1].*field - table[i].*field) / h[i];
    }

    // Interior: Steffen formula
    for (size_t i = 1; i < n-1; ++i) {
        double p = (s[i-1] * h[i] + s[i] * h[i-1]) / (h[i-1] + h[i]);
        // Monotonicity constraint built into slope
        double sign_val = (std::signbit(s[i-1]) == std::signbit(s[i])) ?
            (s[i-1] >= 0 ? 1.0 : -1.0) : 0.0;
        if (sign_val == 0.0) {
            slopes[i] = 0.0;
        } else {
            slopes[i] = sign_val * std::min({std::abs(s[i-1]), std::abs(s[i]), 0.5 * std::abs(p)});
        }
    }

    // Endpoints: Steffen's one-sided formula with monotonicity check
    slopes[0] = s[0]; // simple one-sided
    if (n >= 3) {
        double p0 = s[0] + (s[0] - s[1]) * h[0] / (h[0] + h[1]);
        if (p0 * s[0] <= 0.0) {
            slopes[0] = 0.0;
        } else if (std::abs(p0) > 2.0 * std::abs(s[0])) {
            slopes[0] = 2.0 * s[0];
        } else {
            slopes[0] = p0;
        }
    }

    slopes[n-1] = s[n-2]; // simple one-sided
    if (n >= 3) {
        double pn = s[n-2] + (s[n-2] - s[n-3]) * h[n-2] / (h[n-2] + h[n-3]);
        if (pn * s[n-2] <= 0.0) {
            slopes[n-1] = 0.0;
        } else if (std::abs(pn) > 2.0 * std::abs(s[n-2])) {
            slopes[n-1] = 2.0 * s[n-2];
        } else {
            slopes[n-1] = pn;
        }
    }

    return slopes;
}

// Hermite evaluation is the same as PCHIP
static double interp(const std::vector<DiskRow>& table, const std::vector<double>& slopes,
                     double R, double DiskRow::*field) {
    return pchip::interp(table, slopes, R, field);
}

} // namespace steffen

// ============================================================================
// Test infrastructure
// ============================================================================
struct TestStats {
    const char* method_name;
    const char* field_name;
    double max_rel_diff;       // max |new - old| / |old| vs current method
    double mean_rel_diff;      // mean of same
    double max_abs_diff;       // max |new - old|
    size_t n_negative;         // count of negative values for positive-definite fields
    size_t n_unbounded;        // count of values outside [min(y0,y1), max(y0,y1)] bracket
    double max_gridpoint_err;  // max error when evaluating exactly at grid points
    size_t n_samples;
};

static void print_stats(const TestStats& s) {
    std::cout << std::left << std::setw(16) << s.method_name
              << std::setw(10) << s.field_name
              << std::scientific << std::setprecision(3)
              << "  max_rel=" << std::setw(12) << s.max_rel_diff
              << "  mean_rel=" << std::setw(12) << s.mean_rel_diff
              << "  max_abs=" << std::setw(12) << s.max_abs_diff
              << std::fixed
              << "  neg=" << std::setw(6) << s.n_negative
              << "  unbnd=" << std::setw(6) << s.n_unbounded
              << "  gp_err=" << std::scientific << std::setw(12) << s.max_gridpoint_err
              << "  (N=" << s.n_samples << ")\n";
}

// ============================================================================
// Main test driver
// ============================================================================
int main(int argc, char** argv) {
    std::string csv_path = "../../SpaceHub/src/interaction/disk_tab/disk_default.csv";
    if (argc > 1) csv_path = argv[1];

    std::cout << "Loading disk table: " << csv_path << "\n";
    auto table = load_disk(csv_path);
    std::cout << "Loaded " << table.size() << " rows, R range: ["
              << table.front().R << ", " << table.back().R << "]\n\n";

    if (table.size() < 4) {
        std::cerr << "ERROR: Table too small for cubic interpolation (need >= 4 rows)\n";
        return 1;
    }

    // Precompute slopes for PCHIP and Steffen (per field)
    struct PrecomputedSlopes {
        std::vector<double> pchip_slopes;
        std::vector<double> steffen_slopes;
    };
    PrecomputedSlopes precomp[N_FIELDS];
    for (size_t f = 0; f < N_FIELDS; ++f) {
        precomp[f].pchip_slopes   = pchip::compute_slopes(table, fields[f].ptr);
        precomp[f].steffen_slopes = steffen::compute_slopes(table, fields[f].ptr);
    }

    // Number of sub-samples per interval
    constexpr int N_SUB = 20;

    // Collect stats for each method x field
    // Methods: 0=loglog_linear, 1=pchip, 2=steffen
    const char* method_names[] = {"LogLog-Linear", "PCHIP", "Steffen"};
    constexpr int N_METHODS = 3;

    TestStats stats[N_METHODS][N_FIELDS];
    for (int m = 0; m < N_METHODS; ++m)
        for (size_t f = 0; f < N_FIELDS; ++f) {
            stats[m][f] = {method_names[m], fields[f].name, 0, 0, 0, 0, 0, 0, 0};
        }

    // ---- Test 1: Sub-interval sampling ----
    std::cout << "=== Test 1: Sub-interval comparison (N_SUB=" << N_SUB << " per interval) ===\n";

    for (size_t i = 0; i < table.size() - 1; ++i) {
        double R0 = table[i].R, R1 = table[i+1].R;

        for (int s = 1; s < N_SUB; ++s) {  // skip endpoints (tested separately)
            double frac = static_cast<double>(s) / N_SUB;
            // Use geometric mean for log-spaced grids
            double R = R0 * std::pow(R1 / R0, frac);

            for (size_t f = 0; f < N_FIELDS; ++f) {
                double y_current = current::interp(table, R, fields[f].ptr);

                double y_candidates[N_METHODS];
                y_candidates[0] = loglog_linear::interp(table, R, fields[f].ptr, fields[f].expect_positive);
                y_candidates[1] = pchip::interp(table, precomp[f].pchip_slopes, R, fields[f].ptr);
                y_candidates[2] = steffen::interp(table, precomp[f].steffen_slopes, R, fields[f].ptr);

                double y0 = table[i].*(fields[f].ptr);
                double y1 = table[i+1].*(fields[f].ptr);
                double bracket_min = std::min(y0, y1);
                double bracket_max = std::max(y0, y1);

                for (int m = 0; m < N_METHODS; ++m) {
                    double y_new = y_candidates[m];
                    stats[m][f].n_samples++;

                    // Relative difference from current method
                    double denom = std::abs(y_current);
                    if (denom > 1e-300) {
                        double rel = std::abs(y_new - y_current) / denom;
                        stats[m][f].max_rel_diff = std::max(stats[m][f].max_rel_diff, rel);
                        stats[m][f].mean_rel_diff += rel;
                    }
                    double abs_diff = std::abs(y_new - y_current);
                    stats[m][f].max_abs_diff = std::max(stats[m][f].max_abs_diff, abs_diff);

                    // Positivity check
                    if (fields[f].expect_positive && y_new < 0.0) {
                        stats[m][f].n_negative++;
                    }

                    // Boundedness check (is it within the bracket of endpoints?)
                    // Allow 5% overshoot tolerance for cubic methods
                    double range = bracket_max - bracket_min;
                    double tol = 0.05 * range;
                    if (y_new < bracket_min - tol || y_new > bracket_max + tol) {
                        stats[m][f].n_unbounded++;
                    }
                }
            }
        }
    }

    // Finalize mean
    for (int m = 0; m < N_METHODS; ++m)
        for (size_t f = 0; f < N_FIELDS; ++f)
            if (stats[m][f].n_samples > 0)
                stats[m][f].mean_rel_diff /= stats[m][f].n_samples;

    // ---- Test 2: Grid-point accuracy ----
    std::cout << "=== Test 2: Grid-point accuracy ===\n";

    for (size_t i = 0; i < table.size(); ++i) {
        double R = table[i].R;
        for (size_t f = 0; f < N_FIELDS; ++f) {
            double y_exact = table[i].*(fields[f].ptr);

            double y_candidates[N_METHODS];
            y_candidates[0] = loglog_linear::interp(table, R, fields[f].ptr, fields[f].expect_positive);
            y_candidates[1] = pchip::interp(table, precomp[f].pchip_slopes, R, fields[f].ptr);
            y_candidates[2] = steffen::interp(table, precomp[f].steffen_slopes, R, fields[f].ptr);

            for (int m = 0; m < N_METHODS; ++m) {
                double err = std::abs(y_candidates[m] - y_exact);
                if (std::abs(y_exact) > 1e-300)
                    err /= std::abs(y_exact); // relative
                stats[m][f].max_gridpoint_err = std::max(stats[m][f].max_gridpoint_err, err);
            }
        }
    }

    // ---- Print results ----
    std::cout << "\n=== RESULTS ===\n\n";
    std::cout << std::left << std::setw(16) << "Method"
              << std::setw(10) << "Field"
              << "  " << std::setw(16) << "max_rel_diff"
              << "  " << std::setw(16) << "mean_rel_diff"
              << "  " << std::setw(16) << "max_abs_diff"
              << "  " << std::setw(10) << "neg"
              << "  " << std::setw(10) << "unbnd"
              << "  " << std::setw(16) << "gp_err"
              << "\n";
    std::cout << std::string(120, '-') << "\n";

    bool any_fail = false;
    for (int m = 0; m < N_METHODS; ++m) {
        for (size_t f = 0; f < N_FIELDS; ++f) {
            print_stats(stats[m][f]);
            if (fields[f].expect_positive && stats[m][f].n_negative > 0) {
                std::cout << "  *** FAIL: " << stats[m][f].n_negative
                          << " negative values in positive-definite field!\n";
                any_fail = true;
            }
            if (stats[m][f].max_gridpoint_err > 1e-10) {
                std::cout << "  *** WARNING: Grid-point error > 1e-10\n";
            }
        }
        std::cout << "\n";
    }

    // ---- Test 3: Boundary behavior ----
    std::cout << "=== Test 3: Near-boundary queries ===\n";
    {
        double Rmin = table.front().R;
        double Rmax = table.back().R;

        // Query very close to boundaries
        double R_near_min = Rmin * 1.00001;
        double R_near_max = Rmax * 0.99999;

        for (double R_test : {R_near_min, R_near_max}) {
            std::cout << "  R = " << R_test << (R_test < Rmin * 1.001 ? " (near Rmin)" : " (near Rmax)") << ":\n";
            for (size_t f = 0; f < N_FIELDS; ++f) {
                double y_cur = current::interp(table, R_test, fields[f].ptr);
                double y_ll  = loglog_linear::interp(table, R_test, fields[f].ptr, fields[f].expect_positive);
                double y_pc  = pchip::interp(table, precomp[f].pchip_slopes, R_test, fields[f].ptr);
                double y_st  = steffen::interp(table, precomp[f].steffen_slopes, R_test, fields[f].ptr);

                std::cout << "    " << std::setw(8) << fields[f].name
                          << "  current=" << std::scientific << std::setprecision(6) << y_cur
                          << "  loglog=" << y_ll
                          << "  pchip=" << y_pc
                          << "  steffen=" << y_st << "\n";

                // Check for negatives
                if (fields[f].expect_positive) {
                    if (y_ll < 0) { std::cout << "    *** FAIL: loglog negative near boundary\n"; any_fail = true; }
                    if (y_pc < 0) { std::cout << "    *** FAIL: pchip negative near boundary\n"; any_fail = true; }
                    if (y_st < 0) { std::cout << "    *** FAIL: steffen negative near boundary\n"; any_fail = true; }
                }
            }
        }
    }

    // ---- Test 4: Zone boundary region (look for largest gradient jump) ----
    std::cout << "\n=== Test 4: Zone boundary / largest gradient jump region ===\n";
    {
        // Find the interval with the largest relative change in rho
        size_t worst_idx = 0;
        double worst_ratio = 0;
        for (size_t i = 0; i < table.size() - 1; ++i) {
            double r = std::abs(table[i+1].rho / table[i].rho - 1.0);
            if (r > worst_ratio) { worst_ratio = r; worst_idx = i; }
        }

        std::cout << "  Largest rho transition at interval " << worst_idx
                  << " (R/Rg=" << table[worst_idx].R_Rg << " to " << table[worst_idx+1].R_Rg << ")"
                  << ", ratio=" << std::fixed << std::setprecision(4) << worst_ratio << "\n";

        // Sample densely around this region
        size_t start = (worst_idx > 5) ? worst_idx - 5 : 0;
        size_t end = std::min(worst_idx + 6, table.size() - 1);

        std::cout << "  Dense sampling around transition (R/Rg range "
                  << table[start].R_Rg << " to " << table[end].R_Rg << "):\n";

        int n_neg_ll = 0, n_neg_pc = 0, n_neg_st = 0;
        double max_overshoot_pc = 0, max_overshoot_st = 0;

        for (size_t i = start; i < end; ++i) {
            double R0 = table[i].R, R1 = table[i+1].R;
            for (int s = 1; s < N_SUB; ++s) {
                double frac = static_cast<double>(s) / N_SUB;
                double R = R0 * std::pow(R1 / R0, frac);

                double y_ll = loglog_linear::interp(table, R, &DiskRow::rho, true);
                double y_pc = pchip::interp(table, precomp[0].pchip_slopes, R, &DiskRow::rho);
                double y_st = steffen::interp(table, precomp[0].steffen_slopes, R, &DiskRow::rho);

                if (y_ll < 0) n_neg_ll++;
                if (y_pc < 0) n_neg_pc++;
                if (y_st < 0) n_neg_st++;

                // Check overshoot beyond bracket
                double bracket_min = std::min(table[i].rho, table[i+1].rho);
                double bracket_max = std::max(table[i].rho, table[i+1].rho);
                if (bracket_max > bracket_min) {
                    double range = bracket_max - bracket_min;
                    double os_pc = std::max(0.0, std::max(bracket_min - y_pc, y_pc - bracket_max)) / range;
                    double os_st = std::max(0.0, std::max(bracket_min - y_st, y_st - bracket_max)) / range;
                    max_overshoot_pc = std::max(max_overshoot_pc, os_pc);
                    max_overshoot_st = std::max(max_overshoot_st, os_st);
                }
            }
        }

        std::cout << "  Negatives in transition zone:  loglog=" << n_neg_ll
                  << "  pchip=" << n_neg_pc << "  steffen=" << n_neg_st << "\n";
        std::cout << "  Max overshoot (fraction of bracket range):  pchip="
                  << std::scientific << max_overshoot_pc << "  steffen=" << max_overshoot_st << "\n";
    }

    // ---- Test 5: Smoothness proxy (estimate derivative discontinuity) ----
    std::cout << "\n=== Test 5: Smoothness (derivative continuity at grid points) ===\n";
    {
        // For each method, estimate d/dR at grid point i from left and right
        // and measure the jump |d_left - d_right| / |d_avg|
        for (int m = 0; m < N_METHODS; ++m) {
            double max_jump[N_FIELDS] = {};
            double mean_jump[N_FIELDS] = {};

            for (size_t i = 1; i < table.size() - 1; ++i) {
                double R = table[i].R;
                double eps_left  = (table[i].R - table[i-1].R) * 1e-4;
                double eps_right = (table[i+1].R - table[i].R) * 1e-4;

                for (size_t f = 0; f < N_FIELDS; ++f) {
                    double y_at = 0, y_left = 0, y_right = 0;

                    if (m == 0) {
                        y_at    = loglog_linear::interp(table, R, fields[f].ptr, fields[f].expect_positive);
                        y_left  = loglog_linear::interp(table, R - eps_left, fields[f].ptr, fields[f].expect_positive);
                        y_right = loglog_linear::interp(table, R + eps_right, fields[f].ptr, fields[f].expect_positive);
                    } else if (m == 1) {
                        y_at    = pchip::interp(table, precomp[f].pchip_slopes, R, fields[f].ptr);
                        y_left  = pchip::interp(table, precomp[f].pchip_slopes, R - eps_left, fields[f].ptr);
                        y_right = pchip::interp(table, precomp[f].pchip_slopes, R + eps_right, fields[f].ptr);
                    } else {
                        y_at    = steffen::interp(table, precomp[f].steffen_slopes, R, fields[f].ptr);
                        y_left  = steffen::interp(table, precomp[f].steffen_slopes, R - eps_left, fields[f].ptr);
                        y_right = steffen::interp(table, precomp[f].steffen_slopes, R + eps_right, fields[f].ptr);
                    }

                    double d_left  = (y_at - y_left) / eps_left;
                    double d_right = (y_right - y_at) / eps_right;
                    double d_avg = 0.5 * (std::abs(d_left) + std::abs(d_right));

                    if (d_avg > 1e-300) {
                        double jump = std::abs(d_right - d_left) / d_avg;
                        max_jump[f] = std::max(max_jump[f], jump);
                        mean_jump[f] += jump;
                    }
                }
            }

            std::cout << "  " << method_names[m] << ":\n";
            for (size_t f = 0; f < N_FIELDS; ++f) {
                size_t n = table.size() - 2;
                std::cout << "    " << std::setw(8) << fields[f].name
                          << "  max_deriv_jump=" << std::scientific << std::setprecision(3) << max_jump[f]
                          << "  mean_deriv_jump=" << (n > 0 ? mean_jump[f] / n : 0.0) << "\n";
            }
        }
    }

    // ---- Test 6: Performance comparison ----
    std::cout << "\n=== Test 6: Performance (timing 1M random lookups) ===\n";
    {
        constexpr size_t N_PERF = 1000000;
        double Rmin = table.front().R;
        double Rmax = table.back().R;

        // Generate random R values (log-uniform)
        std::vector<double> R_vals(N_PERF);
        for (size_t i = 0; i < N_PERF; ++i) {
            double frac = static_cast<double>(i) / N_PERF;
            R_vals[i] = Rmin * std::pow(Rmax / Rmin, frac);
        }

        // Time current method
        volatile double sink = 0;
        auto t0 = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < N_PERF; ++i) {
            sink += current::interp(table, R_vals[i], &DiskRow::rho);
            sink += current::interp(table, R_vals[i], &DiskRow::cs);
            sink += current::interp(table, R_vals[i], &DiskRow::H);
            sink += current::interp(table, R_vals[i], &DiskRow::grad_P);
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double dt_current = std::chrono::duration<double, std::milli>(t1 - t0).count();

        // Time loglog linear
        t0 = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < N_PERF; ++i) {
            sink += loglog_linear::interp(table, R_vals[i], &DiskRow::rho, true);
            sink += loglog_linear::interp(table, R_vals[i], &DiskRow::cs, true);
            sink += loglog_linear::interp(table, R_vals[i], &DiskRow::H, true);
            sink += loglog_linear::interp(table, R_vals[i], &DiskRow::grad_P, false);
        }
        t1 = std::chrono::high_resolution_clock::now();
        double dt_loglog = std::chrono::duration<double, std::milli>(t1 - t0).count();

        // Time PCHIP (includes precomputed slopes)
        t0 = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < N_PERF; ++i) {
            sink += pchip::interp(table, precomp[0].pchip_slopes, R_vals[i], &DiskRow::rho);
            sink += pchip::interp(table, precomp[1].pchip_slopes, R_vals[i], &DiskRow::cs);
            sink += pchip::interp(table, precomp[2].pchip_slopes, R_vals[i], &DiskRow::H);
            sink += pchip::interp(table, precomp[3].pchip_slopes, R_vals[i], &DiskRow::grad_P);
        }
        t1 = std::chrono::high_resolution_clock::now();
        double dt_pchip = std::chrono::duration<double, std::milli>(t1 - t0).count();

        // Time Steffen
        t0 = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < N_PERF; ++i) {
            sink += steffen::interp(table, precomp[0].steffen_slopes, R_vals[i], &DiskRow::rho);
            sink += steffen::interp(table, precomp[1].steffen_slopes, R_vals[i], &DiskRow::cs);
            sink += steffen::interp(table, precomp[2].steffen_slopes, R_vals[i], &DiskRow::H);
            sink += steffen::interp(table, precomp[3].steffen_slopes, R_vals[i], &DiskRow::grad_P);
        }
        t1 = std::chrono::high_resolution_clock::now();
        double dt_steffen = std::chrono::duration<double, std::milli>(t1 - t0).count();

        std::cout << "  Current:       " << std::fixed << std::setprecision(1) << dt_current << " ms\n";
        std::cout << "  LogLog-Linear: " << dt_loglog << " ms (" << std::setprecision(2) << dt_loglog/dt_current << "x)\n";
        std::cout << "  PCHIP:         " << dt_pchip << " ms (" << dt_pchip/dt_current << "x)\n";
        std::cout << "  Steffen:       " << dt_steffen << " ms (" << dt_steffen/dt_current << "x)\n";
        (void)sink; // prevent optimizer from removing the loop
    }

    // ---- Summary ----
    std::cout << "\n=== SUMMARY ===\n";
    if (any_fail) {
        std::cout << "*** FAILURES DETECTED - see above for details ***\n";
    } else {
        std::cout << "All positivity and boundary checks PASSED.\n";
    }

    // Quick comparison table
    std::cout << "\nMethod comparison (rho field):\n";
    std::cout << std::left << std::setw(18) << "Method"
              << std::setw(14) << "max_rel_diff"
              << std::setw(14) << "mean_rel_diff"
              << std::setw(8) << "neg"
              << std::setw(8) << "unbnd"
              << "smoothness\n";
    std::cout << std::string(74, '-') << "\n";
    for (int m = 0; m < N_METHODS; ++m) {
        std::cout << std::setw(18) << method_names[m]
                  << std::scientific << std::setprecision(3)
                  << std::setw(14) << stats[m][0].max_rel_diff
                  << std::setw(14) << stats[m][0].mean_rel_diff
                  << std::fixed
                  << std::setw(8) << stats[m][0].n_negative
                  << std::setw(8) << stats[m][0].n_unbounded
                  << (m == 0 ? "C0" : "C1") << "\n";
    }

    return any_fail ? 1 : 0;
}
