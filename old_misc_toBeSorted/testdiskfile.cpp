#include "src/spaceHub.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace hub;
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

struct DiskRow {
    double R, R_Rg, Tc, rho, P, cs, H, visc, Sigma, Q;
    std::string zone;
};


//------Helper functions-----------
double interp(const std::vector<DiskRow>& table, double R, double DiskRow::*field) {
    // Find interval, works because R is monotonic
    auto it = std::lower_bound(table.begin(), table.end(), R,
        [](const DiskRow& row, double r) { return row.R < r; });

    size_t i = std::clamp<size_t>(it - table.begin(), 1, table.size() - 2); //Handles interpolation around extreme Rs

    // Get 4 neighbors for cubic Hermite (Catmull-Rom)
    size_t i0 = (i > 1) ? i - 1 : 0; // prevent negative index
    size_t i1 = i; 
    size_t i2 = i + 1;
    size_t i3 = std::min(i + 2, table.size() - 1); //handle running into upper bound

    double x0 = table[i0].R, x1 = table[i1].R, x2 = table[i2].R, x3 = table[i3].R; //get 4 r values
    double y0 = table[i0].*field, y1 = table[i1].*field, y2 = table[i2].*field, y3 = table[i3].*field; //get 4 dependent variable values

    // slope 1 and slope 2 for Catmull-Rom spline
    double m1 = (y2 - y0) / (x2 - x0);
    double m2 = (y3 - y1) / (x3 - x1);

    // calc Hermite basis
    double h = x2 - x1; //interval width
    double t = (R - x1) / h; //normalized position in interval (t=0 at x1, t=1 at x2)
    double t2 = t * t, t3 = t2 * t;

    // Hermite basis functions
    double h00 = 2*t3 - 3*t2 + 1;
    double h10 = t3 - 2*t2 + t;
    double h01 = -2*t3 + 3*t2;
    double h11 = t3 - t2;

    return h00*y1 + h10*h*m1 + h01*y2 + h11*h*m2;
}


//----------Called on init to load pre-tabulated disk file to memory--------
std::vector<DiskRow> load_disk_table(const std::string& filename) {
    std::vector<DiskRow> table;
    std::ifstream file(filename);
    std::string line;
    std::getline(file, line); // skip header

    //Note this is dependent on the ordering of the columns in the csv file, so if that ever changes this will need to be updated
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        DiskRow row;
        char comma;
        ss >> row.R >> comma >> row.R_Rg >> comma >> row.Tc >> comma
           >> row.rho >> comma >> row.P >> comma >> row.cs >> comma
           >> row.H >> comma >> row.visc >> comma >> row.Sigma >> comma
           >> row.Q >> comma;
        std::getline(ss, row.zone);
        table.push_back(row);
    }
    return table;
}

int main(int argc, char **argv) {
    print(std::cout, "disk tab generation test (log)\n");

    auto disk = load_disk_table("src/interaction/disk_tab/disk_test_log.csv");
    print(std::cout, "Loaded ", disk.size(), " rows\n");

    Vec3<Scalar> pos = Vec3<Scalar>(1000, 1000, 25);
    auto r = sqrt(pos.x*pos.x + pos.y*pos.y);
    auto z = pos.z;

    

    double rho_interp = interp(disk, r, &DiskRow::rho);
    double cs_interp = interp(disk, r, &DiskRow::cs);
    double H_interp = interp(disk, r, &DiskRow::H);

    auto rho_z = rho_interp * exp(-1 * (z*z) / (2*H_interp * H_interp));

    print(std::cout, "At R=", r, ", z=", z, ", H=", H_interp, "\n");
    print(std::cout, "rho_c=", rho_interp, ", rho_z=", rho_z, ", cs=", cs_interp, "\n");

    return 0;
}
//g++-std=c++17 -O3 -pthread testdiskfile.cpp -o testdiskfile
// ./testdiskfile