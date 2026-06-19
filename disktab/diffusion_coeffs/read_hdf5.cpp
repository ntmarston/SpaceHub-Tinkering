#include "../../SpaceHub/src/spaceHub.hpp"
#include <highfive/H5File.hpp>
#include <random>
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
using namespace HighFive;

using f = Interactions<NewtonianGrav, PN1>;
using Solver = methods::Sym6<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;
auto G = consts::G;
auto pi = consts::pi;

/**
 * Return the HDF5 group (slab) whose attributes exactly match the four physical parameters.
 *
 * Iterates all top-level groups in the file and returns the first one that satisfies
 * all four attribute conditions. Integer attributes are compared exactly; float attributes
 * use a tolerance of 1e-6 because HDF5 stores them as float64.
 *
 * @param file          Open HDF5 file to search.
 * @param M_exp         Black-hole mass exponent (log10 M_BH / M_sun), matched as int64.
 * @param m_star        Stellar mass in solar masses, default 0.3.
 * @param gamma         Stellar density power-law index, default 1.75.
 * @param r_inf_scaled  Influence radius scaling factor, default 1.0.
 *
 * @return The matching HDF5 group.
 * @throws std::runtime_error if no group matches all four conditions.
 */
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

// All standard-integral tables for one slab, plus the axes they're sampled on.
// Tables are indexed [i_r][i_v]. r is log-spaced; v is linear-spaced.
// Lookup constants (t0, inv_step) are hoisted at construction for O(1) cell().
template <typename Scalar = double>
struct StdIntegrals {
    std::vector<Scalar> r_grid, v_grid;
    Scalar r_t0 = 0, r_inv_step = 0;  // log-space hoist for r
    Scalar v_t0 = 0, v_inv_step = 0;  // linear-space hoist for v
    std::vector<std::vector<Scalar>> E1, F2, F4;

    /**
     * Store the radius grid and precompute log-space lookup constants.
     *
     * @param g  Radius grid in AU (SpaceHub natural length). Must be log-uniformly spaced
     *           and strictly increasing.
     */
    void set_r(std::vector<Scalar> g) {
        r_grid = std::move(g);
        r_t0 = std::log(r_grid.front());
        r_inv_step = (r_grid.size() - 1) / (std::log(r_grid.back()) - r_t0);
    }

    /**
     * Store the normalized-speed grid and precompute linear-space lookup constants.
     *
     * @param g  Normalized-speed grid, x = v/v_esc, dimensionless. Must be uniformly
     *           spaced and strictly increasing.
     */
    void set_v(std::vector<Scalar> g) {
        v_grid = std::move(g);
        v_t0 = v_grid.front();
        v_inv_step = (v_grid.size() - 1) / (v_grid.back() - v_t0);
    }

    /**
     * Return the lower index of the log-space cell bracketing r, clamped to [0, n-2].
     *
     * @param r  Query radius in AU. Assumed in-range; clamping handles boundary rounding.
     * @return   Index i such that r_grid[i] <= r < r_grid[i+1] (clamped at the last cell).
     */
    size_t r_cell(Scalar r) const {
        size_t i = static_cast<size_t>((std::log(r) - r_t0) * r_inv_step);
        return i < r_grid.size() - 1 ? i : r_grid.size() - 2;
    }

    /**
     * Return the lower index of the linear-space cell bracketing v, clamped to [0, n-2].
     *
     * @param v  Query normalized speed x = v/v_esc. Assumed in-range; clamping handles
     *           boundary rounding.
     * @return   Index j such that v_grid[j] <= v < v_grid[j+1] (clamped at the last cell).
     */
    size_t v_cell(Scalar v) const {
        size_t i = static_cast<size_t>((v - v_t0) * v_inv_step);
        return i < v_grid.size() - 1 ? i : v_grid.size() - 2;
    }
};

/**
 * Read the E1, F2, F4 standard-integral tables and their axes from an HDF5 slab group.
 *
 * Datasets "r" and "v_norm" are read into the StdIntegrals axes via set_r/set_v,
 * which also precompute the O(1) cell-lookup constants. The three table datasets
 * are read directly into the 2-D vectors.
 *
 * @tparam Scalar  Float type for all grid and table data. Defaults to double.
 *
 * @param slab  Open HDF5 group containing datasets: "r" (radius, AU), "v_norm"
 *              (normalized speed, dimensionless), "E1", "F2", "F4" (standard
 *              integrals in SpaceHub natural units), all indexed [i_r][i_v].
 *
 * @return Fully populated StdIntegrals struct, ready for interp_table calls.
 */
template <typename Scalar = double>
StdIntegrals<Scalar> extract_Std_Integrals(HighFive::Group& slab) {
    StdIntegrals<Scalar> s;
    std::vector<Scalar> r, v;
    slab.getDataSet("r").read(r);
    slab.getDataSet("v_norm").read(v);
    s.set_r(std::move(r));
    s.set_v(std::move(v));
    slab.getDataSet("E1").read(s.E1);
    slab.getDataSet("F2").read(s.F2);
    slab.getDataSet("F4").read(s.F4);
    return s;
}

/**
 * Bilinear interpolation of a diffusion-coefficient table at (r, v).
 * /!\NOTE v is particle velocity SCALED BY escape velocity.
 * Blends along v on the two bracketing r-rows, then along r. Cell lookup is O(1)
 * via the hoisted constants stored in si -- no binary search.
 *
 * @tparam Scalar  Float type. double in practice; bump it if SpaceHub wants more.
 *
 * @param si     StdIntegrals holding the axes and their precomputed lookup constants.
 *               The r and v grid sizes must match the outer/inner dimensions of table.
 * @param table  One standard-integral slab (E1, F2, or F4), indexed [i_r][i_v]
 *               (outer = radius, inner = normalized speed).
 * @param r      Query radius in AU.
 * @param v      Query speed -- ALREADY normalized: x = v/v_esc with
 *               v_esc = sqrt(2 G M / r), G=1. Yes, dimensionless. No, not m/s.
 *               Normalize it yourself; you need M to do it, and we don't have M here.
 *
 * @return Interpolated value, in SpaceHub-natural units (yr/(2*pi))^2 / AU^5 --
 *         the tables are stored that way, so it's returned as-is, no conversion.
 *         Returns exactly 0 if r or v lands on or outside the tabulated range
 *         (silent, like disk-model.hpp::interp_all -- so don't trust a bare 0).
 */
template <typename Scalar>
Scalar interp_table(const StdIntegrals<Scalar>& si,
                    const std::vector<std::vector<Scalar>>& table,
                    Scalar r, Scalar v) {
    const auto& rg = si.r_grid;
    const auto& vg = si.v_grid;
    if (r <= rg.front() || r >= rg.back()) return 0;
    if (v <= vg.front() || v >= vg.back()) return 0;

    const size_t i = si.r_cell(r);  // O(1), uses the hoisted constants
    const size_t j = si.v_cell(v);

    // Blend along v on each bracketing r-row, then along r.
    const Scalar fr = (r - rg[i]) / (rg[i + 1] - rg[i]);
    const Scalar fv = (v - vg[j]) / (vg[j + 1] - vg[j]);
    auto lerp = [](Scalar f0, Scalar f1, Scalar t) { return f0 + (f1 - f0) * t; };

    const Scalar f_lo = lerp(table[i][j],     table[i][j + 1],     fv);
    const Scalar f_hi = lerp(table[i + 1][j], table[i + 1][j + 1], fv);
    return lerp(f_lo, f_hi, fr);
}


Vec3<Scalar> delta_vel(auto r_vec, Scalar v_vec, auto dt){
    static inline Scalar M_exp = 7; //use a global in the real sim
    Scalar m1 = pow(10, M_exp) * unit::Ms; //access from ptc in real sim
    Scalar m2 = 10_Ms; //access from ptc in real sim
     // diff time with previous step in real sim to get dt
    static inline float mf = 0.3_Ms;
    static inline float gamma = 7/4;
    auto v = norm(v_vec);
    auto v_hat = v_vec / norm(v_vec);
    
    auto r = norm(r_vec);
    auto r_hat = r_vec / norm(r_vec);

    auto n_hat = cross(r_vec, v_vec) / norm(cross(r_vec, v_vec)); // vhat, rhat, nhat should form orthonormal basis
    
    Scalar v_esc = sqrt(2 * G * m1 / r);
    auto v_norm = v/v_esc;

    StdIntegrals std_ints = extract_Std_Integrals(access_slab(dc_table, M_exp, m_star = mf,
                            gamma = gamma, r_inf_scaled = 1.0));

    //---- coulomb log -----
    auto b_maxmin_ratio = M / gamma / (1+gamma) / (m+mf);
    auto lnLambda = log(b_maxmin_ratio);

    Scalar E1 = interp_table(si = std_ints, table = si.E1, r=r, v=v_norm);
    Scalar F2 = interp_table(si = std_ints, table = si.E1, r=r, v=v_norm);
    Scalar F4 = interp_table(si = std_ints, table = si.E1, r=r, v=v_norm);
    Scalar pi2G2lnLambda = pi * pi * G * G * lnLambda;

    static std::mt19937 rng{std::random_device{}()};
    Scalar X1     = std::uniform_real_distribution<Scalar>{-1,   1  }(rng);
    Scalar X2     = std::uniform_real_distribution<Scalar>{-1,   1  }(rng);
    Scalar theta2 = std::uniform_real_distribution<Scalar>{ 0, 2*pi }(rng);

    auto DC_1par = -16 * pi2G2lnLambda * (mf + m) * mf * F2;
    auto DC_2par = (32/3) * pi2G2lnLambda * (F4 + E1) * v;
    auto DC_2perp = (32/3) * pi2G2lnLambda * (3*F2 - F4 + 2*E1) * v;

    auto dv_i = (DC_1par * dt) + (sqrt(DC_2par * dt) * X1); // v_hat component
    auto dv_j = -sqrt(DC_2perp * dt) * X2 * cos(theta2); // r_hat - negative to point towards central mass 
    auto dv_k = sqrt(DC_2perp * dt) * X2 * sin(theta2); // n_hat
    
    auto delta_v = dv_i * v_hat + dv_j * r_hat * dv_k * n_hat;
    

    return delta_v;
}

int main() {
    

    // Open the file read-only.
    File file("dc_tables.h5", File::ReadOnly);

    Group slab = access_slab(file, 6);
    std::cout << "matched: " << slab.getPath() << "\n";

    StdIntegrals si = extract_Std_Integrals(slab);


    //===================PRINTING FOR DEBUG=======================
    // Print the E1 table as an aligned grid: top row is the v-axis, each subsequent
    // row is r then E1(r, v). Fixed-width columns keep it readable in a terminal.
    std::cout << "E1 table (" << si.r_grid.size() << " r x " << si.v_grid.size() << " v):\n";

    std::cout << std::setw(10) << "r \\ v" << std::fixed << std::setprecision(3);
    for (double v : si.v_grid) std::cout << std::setw(10) << v;
    std::cout << '\n';

    for (size_t i = 0; i < si.r_grid.size(); ++i) {
        std::cout << std::fixed << std::setprecision(3) << std::setw(10) << si.r_grid[i];
        std::cout << std::scientific << std::setprecision(2);
        for (double e : si.E1[i]) std::cout << std::setw(10) << e;
        std::cout << '\n';
    }

    //===================INTERPOLATION CHECK======================
    std::cout << std::scientific << std::setprecision(6);
    // Exact hit on a grid node should return the stored value (sanity check).
    std::cout << "E1 at node (r[2], v[3]) = "
              << interp_table(si, si.E1, si.r_grid[2], si.v_grid[3])
              << "  (table = " << si.E1[2][3] << ")\n";
    // Off-grid query.
    std::cout << "E1(r=0.5, v=0.30)       = "
              << interp_table(si, si.E1, 0.5, 0.30) << "\n";

    return 0;
}

// Build (HighFive headers vendored at project root; HDF5 C lib is system-installed):
// g++ -std=c++17 -O3 -pthread -I extern disktab/diffusion_coeffs/read_hdf5.cpp -lhdf5 -o disktab/diffusion_coeffs/read_hdf5 (run from repo root)
// Run from disktab/diffusion_coeffs/ so "dc_tables.h5" resolves:
//   cd disktab/diffusion_coeffs && ./read_hdf5
