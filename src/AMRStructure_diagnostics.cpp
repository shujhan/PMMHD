#include "AMRStructure.hpp"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

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

// Half width at half maximum of f along one mesh line (coordinate s), e.g. the
// column x = x_center or the row y = 0. The line is picked out by |t - t0| < tol
// (t = the other coordinate); shared panel nodes are merged as in midplane_flux.
// Walks out from the maximum to the first crossings of f_max/2, linear in between.
// Returns 0 if the line has fewer than 3 nodes.
static double line_hwhm(const std::vector<double>& s_all,
                        const std::vector<double>& t_all,
                        const std::vector<double>& f_all,
                        double t0, double tol_t, double tol_s)
{
    std::vector<std::pair<double,double>> line;
    for (size_t i = 0; i < s_all.size(); ++i) {
        if (std::abs(t_all[i] - t0) < tol_t) {
            line.push_back(std::make_pair(s_all[i], f_all[i]));
        }
    }
    if (line.size() < 3) { return 0.0; }
    std::sort(line.begin(), line.end());

    std::vector<double> s, f;
    size_t i = 0;
    while (i < line.size()) {
        size_t k = i;
        double sum = 0.0;
        while (k < line.size() && line[k].first - line[i].first < tol_s) {
            sum += line[k].second;
            ++k;
        }
        s.push_back(line[i].first);
        f.push_back(sum / (double)(k - i));
        i = k;
    }

    const size_t M = s.size();
    size_t im = 0;
    for (size_t k = 1; k < M; ++k) { if (f[k] > f[im]) { im = k; } }
    const double half = 0.5 * f[im];
    if (half <= 0.0) { return 0.0; }

    double s_lo = s[0], s_hi = s[M-1];
    for (size_t k = im; k > 0; --k) {
        if (f[k-1] < half) {
            s_lo = s[k-1] + (half - f[k-1]) * (s[k] - s[k-1]) / (f[k] - f[k-1]);
            break;
        }
    }
    for (size_t k = im; k + 1 < M; ++k) {
        if (f[k+1] < half) {
            s_hi = s[k] + (f[k] - half) * (s[k+1] - s[k]) / (f[k] - f[k+1]);
            break;
        }
    }
    return 0.5 * (s_hi - s_lo);
}

MHDDiagnostics AMRStructure::compute_diagnostics() {
    MHDDiagnostics d{};
    d.iter = iter_num;
    d.t    = t;

    const size_t N = weights.size();
    double E_kin = 0.0, E_mag = 0.0, H_C = 0.0;
    double I_j = 0.0, I_w = 0.0;
    double w_max = 0.0, j_max = 0.0;

    #pragma omp parallel for reduction(+:E_kin, E_mag, H_C, I_j, I_w) \
                             reduction(max:w_max, j_max)
    for (size_t i = 0; i < N; ++i) {
        const double wi = weights[i];
        // self-consistent flow only: the external stagnation flow (if any) has
        // box-size-dependent energy and is not part of the dynamics' budget
        const double v1 = u1s[i] - u_ext_x(xs[i]);
        const double v2 = u2s[i] - u_ext_y(ys[i]);
        E_kin += 0.5 * wi * (v1*v1 + v2*v2);
        E_mag += 0.5 * wi * (b1s[i]*b1s[i] + b2s[i]*b2s[i]);
        H_C   +=       wi * (v1*b1s[i] + v2*b2s[i]);
        I_j   +=       wi * j0s[i];
        I_w   +=       wi * w0s[i];
        // peak amplitudes: unweighted, so these are nodal maxima rather than
        // quadrature sums and stay comparable between uniform and AMR meshes
        w_max = std::max(w_max, std::abs(w0s[i]));
        j_max = std::max(j_max, std::abs(j0s[i]));
    }

    d.E_kin = E_kin;
    d.E_mag = E_mag;
    d.E_tot = E_kin + E_mag;
    d.H_C   = H_C;
    d.I_j   = I_j;
    d.I_w   = I_w;
    d.w_max = w_max;
    d.j_max = j_max;

    if (bcs == free_bcs) {
        // midplane_flux assumes periodic x; not meaningful on a free line
        d.Psi_rec = std::nan("");
        d.Psi_res = std::nan("");
    } else {
        midplane_flux(xs, ys, b2s, Lx, Ly, d.Psi_rec, d.Psi_res);
    }

    // sheet size from j: half thickness along the center column, half length along y = 0
    const double x_c = 0.5 * (x_min + x_max);
    d.sheet_a = 2 * line_hwhm(ys, xs, j0s, x_c, 1e-9 * Lx, 1e-9 * Ly) / std::acosh(std::sqrt(2.0));
    d.sheet_b = 2 * line_hwhm(xs, ys, j0s, 0.0, 1e-9 * Ly, 1e-9 * Lx);
    return d;
}

int AMRStructure::write_diagnostics(const MHDDiagnostics& d) {
    static bool header_written = false;
    const std::string path = sim_dir + "simulation_output/diagnostics.csv";

    std::ofstream f;
    if (!header_written) {
        f.open(path, std::ios::out | std::ios::trunc);
        f << "iter,t,E_kin,E_mag,E_tot,H_C,I_j,I_w,Psi_rec,Psi_res,w_max,j_max,sheet_a,sheet_b\n";
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
      << d.Psi_rec << "," << d.Psi_res << ","
      << d.w_max << "," << d.j_max << ","
      << d.sheet_a << "," << d.sheet_b << "\n";
    return 0;
}