#include "AMRStructure.hpp"
#include <fstream>
#include <iomanip>
#include <algorithm>

// Reconnected flux on the midplane y = 0.
//
// With lap psi = j and B = grad^perp psi, the code's kernel convention gives
//     b1 = -d_y psi,   b2 = +d_x psi,
// so along the line y = 0
//     psi(x,0) - psi(x_min,0) = int_{x_min}^{x} b2(x',0) dx'.
// X- and O-points on the midplane are where b2 = 0, i.e. the extrema of psi,
// so Psi_rec = max psi - min psi over the line. psi is single-valued on the
// periodic line only if int over the full period vanishes; that integral is
// returned as Psi_res and should stay at round-off.
//
// y = 0 is a node line at every refinement level (level-0 midpoint, halved by
// each refine), so the row is picked out by |y| < tol_y. Panel-shared nodes are
// duplicated in xs/ys, so equal x are merged before the trapezoid sum.
static void midplane_flux(const std::vector<double>& xs,
                          const std::vector<double>& ys,
                          const std::vector<double>& b2s,
                          double Lx, double Ly, double& Psi_rec, double& Psi_res)
{
    Psi_rec = 0.0;
    Psi_res = 0.0;

    const double tol_y = 1e-9 * Ly;   // picks out the row; y scale is Ly
    const double tol_x = 1e-9 * Lx;   // merges shared nodes; x scale is Lx

    std::vector<std::pair<double,double>> row;   // (x, b2) on y = 0
    for (size_t i = 0; i < xs.size(); ++i) {
        if (std::abs(ys[i]) < tol_y) {
            row.push_back(std::make_pair(xs[i], b2s[i]));
        }
    }
    if (row.size() < 3) {
        static bool warned = false;
        if (!warned) {
            std::cerr << "[diag] no midplane row found at y = 0; Psi_rec disabled"
                      << std::endl;
            warned = true;
        }
        return;
    }

    std::sort(row.begin(), row.end());

    // merge duplicated nodes at panel boundaries
    std::vector<double> x, b2;
    size_t i = 0;
    while (i < row.size()) {
        size_t k = i;
        double sum = 0.0;
        while (k < row.size() && row[k].first - row[i].first < tol_x) {
            sum += row[k].second;
            ++k;
        }
        x.push_back(row[i].first);
        b2.push_back(sum / (double)(k - i));
        i = k;
    }

    // cumulative trapezoid: psi_k = int_{x_0}^{x_k} b2 dx
    const size_t M = x.size();
    std::vector<double> psi(M, 0.0);
    for (size_t k = 1; k < M; ++k) {
        psi[k] = psi[k-1] + 0.5 * (b2[k] + b2[k-1]) * (x[k] - x[k-1]);
    }

    Psi_res = psi[M-1];
    Psi_rec = *std::max_element(psi.begin(), psi.end())
            - *std::min_element(psi.begin(), psi.end());
}

MHDDiagnostics AMRStructure::compute_diagnostics() {
    MHDDiagnostics d{};
    d.iter = iter_num;
    d.t    = t;

    const size_t N = weights.size();
    double E_kin = 0.0, E_mag = 0.0, H_C = 0.0;
    double I_j = 0.0, I_w = 0.0;

    #pragma omp parallel for reduction(+:E_kin, E_mag, H_C, I_j, I_w)
    for (size_t i = 0; i < N; ++i) {
        const double wi = weights[i];
        E_kin += 0.5 * wi * (u1s[i]*u1s[i] + u2s[i]*u2s[i]);
        E_mag += 0.5 * wi * (b1s[i]*b1s[i] + b2s[i]*b2s[i]);
        H_C   +=       wi * (u1s[i]*b1s[i] + u2s[i]*b2s[i]);
        I_j   +=       wi * j0s[i];
        I_w   +=       wi * w0s[i];
    }

    d.E_kin = E_kin;
    d.E_mag = E_mag;
    d.E_tot = E_kin + E_mag;
    d.H_C   = H_C;
    d.I_j   = I_j;
    d.I_w   = I_w;

    midplane_flux(xs, ys, b2s, Lx, Ly, d.Psi_rec, d.Psi_res);
    return d;
}

int AMRStructure::write_diagnostics(const MHDDiagnostics& d) {
    static bool header_written = false;
    const std::string path = sim_dir + "simulation_output/diagnostics.csv";

    std::ofstream f;
    if (!header_written) {
        f.open(path, std::ios::out | std::ios::trunc);
        f << "iter,t,E_kin,E_mag,E_tot,H_C,I_j,I_w,Psi_rec,Psi_res\n";
        header_written = true;
    } else {
        f.open(path, std::ios::out | std::ios::app);
    }

    if (!f.is_open()) {
        std::cerr << "[diag] failed to open " << path << std::endl;
        return 1;
    }

    f << std::setprecision(16);
    f << d.iter << "," << d.t << ","
      << d.E_kin << "," << d.E_mag << "," << d.E_tot << ","
      << d.H_C << "," << d.I_j << "," << d.I_w << ","
      << d.Psi_rec << "," << d.Psi_res << "\n";
    return 0;
}