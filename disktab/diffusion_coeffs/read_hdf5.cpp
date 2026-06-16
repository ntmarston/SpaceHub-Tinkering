#include "../../SpaceHub/src/spaceHub.hpp"
#include <highfive/H5File.hpp>

// Return the one slab whose attributes match all four conditions.
// Attributes are stored as int64/float64, so read them in those native types;
// compare floats with a tolerance (exact == on floats is unsafe).
HighFive::Group access_slab(HighFive::File& file, int M_exp, float m_star = 0.3,
                            float gamma = 1.75, float r_inf_scaled = 1.0) {
    auto approx = [](double a, double b) { return std::abs(a - b) < 1e-6; };

    for (const auto& name : file.listObjectNames()) {
        HighFive::Group g = file.getGroup(name);
        if (g.getAttribute("M_exp").read<long long>()       == M_exp        &&
            approx(g.getAttribute("gamma").read<double>(),        gamma)     &&
            approx(g.getAttribute("m_star_Msun").read<double>(),  m_star)    &&
            approx(g.getAttribute("r_inf_scaled").read<double>(), r_inf_scaled))
            return g;
    }
    throw std::runtime_error("no slab matches the requested attributes");
}

// A tabulation axis: the sample grid plus the constants for an O(1) bracketing-
// cell lookup, computed once at construction (the hoist). Geometric grids (e.g.
// r) are indexed in log space; evenly spaced grids (e.g. v) in linear space.
template <typename Scalar>
struct Axis {
    std::vector<Scalar> grid;
    bool log_spaced = false;
    Scalar t0 = 0;        // transformed first node:  log(grid.front()) or grid.front()
    Scalar inv_step = 0;  // (n - 1) / (transformed last - transformed first)

    Axis() = default;
    Axis(std::vector<Scalar> g, bool is_log) : grid(std::move(g)), log_spaced(is_log) {
        t0 = tx(grid.front());
        inv_step = (grid.size() - 1) / (tx(grid.back()) - t0);
    }
    Scalar tx(Scalar x) const { return log_spaced ? std::log(x) : x; }
    // Lower index of the cell bracketing x, clamped to a valid [i, i+1] pair.
    size_t cell(Scalar x) const {
        size_t i = static_cast<size_t>((tx(x) - t0) * inv_step);
        return i < grid.size() - 1 ? i : grid.size() - 2;
    }
};

// All standard-integral tables for one slab, plus the axes they're sampled on.
// Tables are indexed [i_r][i_v]; the axes carry their precomputed lookup constants.
template <typename Scalar = double>
struct StdIntegrals {
    Axis<Scalar> r_axis, v_axis;
    std::vector<std::vector<Scalar>> E1, F2, F4;
};

// Read the E1, F2, F4 tables (and the r, v axes) for a slab in one shot. The axes
// are wrapped as Axis objects here, so the index constants are computed just once.
template <typename Scalar = double>
StdIntegrals<Scalar> extract_Std_Integrals(HighFive::Group& slab) {
    StdIntegrals<Scalar> s;
    std::vector<Scalar> r, v;
    slab.getDataSet("r").read(r);
    slab.getDataSet("v_norm").read(v);
    s.r_axis = Axis<Scalar>(std::move(r), /*log_spaced=*/true);
    s.v_axis = Axis<Scalar>(std::move(v), /*log_spaced=*/false);
    slab.getDataSet("E1").read(s.E1);
    slab.getDataSet("F2").read(s.F2);
    slab.getDataSet("F4").read(s.F4);
    return s;
}

/**
 * Bilinear interpolation of a diffusion-coefficient table at (r, v).
 *
 * Blends along v on the two bracketing r-rows, then along r. Cell lookup is O(1)
 * via the Axis's hoisted constants -- no binary search.
 *
 * @tparam Scalar  Float type. double in practice; bump it if SpaceHub wants more.
 *
 * @param table   One standard-integral slab (E1, F2, or F4), indexed [i_r][i_v]
 *                (outer = radius, inner = normalized speed). Shape must match
 *                r_axis.grid.size() x v_axis.grid.size().
 * @param r_axis  Radius axis, AU (SpaceHub natural length). Log-spaced grid.
 * @param v_axis  Normalized-speed axis, x = v/v_esc, dimensionless. Linear grid.
 * @param r       Query radius, AU.
 * @param v       Query speed -- ALREADY normalized: x = v/v_esc with
 *                v_esc = sqrt(2 G M / r), G=1. Yes, dimensionless. No, not m/s.
 *                Normalize it yourself; you need M to do it, and we don't have M here.
 *
 * @return Interpolated value, in SpaceHub-natural units (yr/(2*pi))^2 / AU^5 --
 *         the tables are stored that way, so it's returned as-is, no conversion.
 *         Returns exactly 0 if r or v lands on or outside the tabulated range
 *         (silent, like disk-model.hpp::interp_all -- so don't trust a bare 0).
 */
template <typename Scalar>
Scalar interp_table(const std::vector<std::vector<Scalar>>& table,
                    const Axis<Scalar>& r_axis, const Axis<Scalar>& v_axis,
                    Scalar r, Scalar v) {
    const auto& rg = r_axis.grid;
    const auto& vg = v_axis.grid;
    if (r <= rg.front() || r >= rg.back()) return 0;
    if (v <= vg.front() || v >= vg.back()) return 0;

    const size_t i = r_axis.cell(r);  // O(1), uses the hoisted constants
    const size_t j = v_axis.cell(v);

    // Blend along v on each bracketing r-row, then along r.
    const Scalar fr = (r - rg[i]) / (rg[i + 1] - rg[i]);
    const Scalar fv = (v - vg[j]) / (vg[j + 1] - vg[j]);
    auto lerp = [](Scalar f0, Scalar f1, Scalar t) { return f0 + (f1 - f0) * t; };

    const Scalar f_lo = lerp(table[i][j],     table[i][j + 1],     fv);
    const Scalar f_hi = lerp(table[i + 1][j], table[i + 1][j + 1], fv);
    return lerp(f_lo, f_hi, fr);
}

int main() {
    using namespace HighFive;

    // Open the file read-only.
    File file("dc_tables.h5", File::ReadOnly);

    Group slab = access_slab(file, 6);
    std::cout << "matched: " << slab.getPath() << "\n";

    StdIntegrals si = extract_Std_Integrals(slab);


    //===================PRINTING FOR DEBUG=======================
    // Print the E1 table as an aligned grid: top row is the v-axis, each subsequent
    // row is r then E1(r, v). Fixed-width columns keep it readable in a terminal.
    std::cout << "E1 table (" << si.r_axis.grid.size() << " r x " << si.v_axis.grid.size() << " v):\n";

    std::cout << std::setw(10) << "r \\ v" << std::fixed << std::setprecision(3);
    for (double v : si.v_axis.grid) std::cout << std::setw(10) << v;
    std::cout << '\n';

    for (size_t i = 0; i < si.r_axis.grid.size(); ++i) {
        std::cout << std::fixed << std::setprecision(3) << std::setw(10) << si.r_axis.grid[i];
        std::cout << std::scientific << std::setprecision(2);
        for (double e : si.E1[i]) std::cout << std::setw(10) << e;
        std::cout << '\n';
    }

    //===================INTERPOLATION CHECK======================
    std::cout << std::scientific << std::setprecision(6);
    // Exact hit on a grid node should return the stored value (sanity check).
    std::cout << "E1 at node (r[2], v[3]) = "
              << interp_table(si.E1, si.r_axis, si.v_axis, si.r_axis.grid[2], si.v_axis.grid[3])
              << "  (table = " << si.E1[2][3] << ")\n";
    // Off-grid query.
    std::cout << "E1(r=0.5, v=0.30)       = "
              << interp_table(si.E1, si.r_axis, si.v_axis, 0.5, 0.30) << "\n";

    return 0;
}

// Build (HighFive headers vendored at project root; HDF5 C lib is system-installed):
// g++ -std=c++17 -O3 -pthread -I extern disktab/diffusion_coeffs/read_hdf5.cpp -lhdf5 -o disktab/diffusion_coeffs/read_hdf5 (run from repo root)
// Run from disktab/diffusion_coeffs/ so "dc_tables.h5" resolves:
//   cd disktab/diffusion_coeffs && ./read_hdf5
