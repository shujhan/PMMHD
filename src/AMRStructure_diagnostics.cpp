#include "AMRStructure.hpp"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

// Half width at half maximum of f along one mesh line (coordinate s), e.g. the
// column x = x_center or the row y = 0. The line is picked out by |t - t0| < tol
// (t = the other coordinate); shared panel nodes are merged as in midplane_flux.
// Works on |f|, so a sheet of either sign is measured (the Loureiro sheet has
// j < 0 with j > 0 return currents). Walks out from the maximum of |f| to the
// first crossings of |f|_max/2, linear in between.
// Returns 0 if the line has fewer than 3 nodes.
static double line_hwhm(const std::vector<double>& s_all,
                        const std::vector<double>& t_all,
                        const std::vector<double>& f_all,
                        double t0, double tol_t, double tol_s)
{
    std::vector<std::pair<double,double>> line;
    for (size_t i = 0; i < s_all.size(); ++i) {
        if (std::abs(t_all[i] - t0) < tol_t) {
            line.push_back(std::make_pair(s_all[i], std::abs(f_all[i])));
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


// Collect the nodes on one mesh line (|t - t0| < tol_t), sorted in s, with
// panel-shared duplicates (|ds| < tol_s) merged by averaging f.
static void mesh_line(const std::vector<double>& s_all,
                      const std::vector<double>& t_all,
                      const std::vector<double>& f_all,
                      double t0, double tol_t, double tol_s,
                      std::vector<double>& s, std::vector<double>& f)
{
    s.clear(); f.clear();
    std::vector<std::pair<double,double>> line;
    for (size_t i = 0; i < s_all.size(); ++i) {
        if (std::abs(t_all[i] - t0) < tol_t) {
            line.push_back(std::make_pair(s_all[i], f_all[i]));
        }
    }
    std::sort(line.begin(), line.end());
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
}

// Trapezoid x-mean of f over the merged nodes of one mesh line (s, f from mesh_line).
static double line_mean(const std::vector<double>& s, const std::vector<double>& f)
{
    const size_t M = s.size();
    if (M < 2) { return std::nan(""); }
    double area = 0.0;
    for (size_t k = 1; k < M; ++k) { area += 0.5 * (f[k] + f[k-1]) * (s[k] - s[k-1]); }
    return area / (s[M-1] - s[0]);
}

// Reconnected flux, Loureiro et al., PRL 95, 235003 (2005):
//     Psi(t) = psi_0 - psi(t, X),   W = 4 sqrt(Psi / |psi_0''(0)|),   gamma_eff = d ln Psi / dt
// psi_0 = psi_eq(X), the unperturbed flux at the X-point.
//
// psi is psis from evaluate_potential (lap psi = j, zero quadrature-weighted mean).
// The domain mean of psi is conserved by the dynamics (the mean of u.grad psi and
// of j vanish on the periodic box), so the zero-mean gauge is fixed in time, as
// the k = 0 mode of a spectral psi is. No equilibrium-diffusion subtraction, as in
// the paper: Psi includes the eta |j(X)| t drift of the diffusing sheet.
//
// Fixed on the first call (t = 0), on the midplane y = 0:
//   psi_0      = x-mean of psi(x, 0, 0): the seed ~ cos(kx) averages out, leaving psi_eq(0)
//   |psi_0''|  = |x-mean of j(x, 0, 0)| = |j_eq(0)|
//   X-point    = midplane node where psi is extremal opposite to the O-point:
//                the minimum of sgn*psi, sgn = -sign(j_eq(0)) (sheet with j < 0 has psi max at y = 0)
// The X-point node is kept for the whole run (symmetric runs: the primary X-point does not move).
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

    d.Psi = d.W = d.gamma_eff = std::nan("");
    const double tol_x = 1e-9 * Lx, tol_y = 1e-9 * Ly;

    if (bcs != free_bcs && psis.size() == xs.size()) {
        if (!rec_ref_set) {
            std::vector<double> x, psi_line, j_line;
            mesh_line(xs, ys, psis, 0.0, tol_y, tol_x, x, psi_line);
            mesh_line(xs, ys, j0s,  0.0, tol_y, tol_x, x, j_line);
            const double j_eq0 = line_mean(x, j_line);
            if (x.size() >= 3 && std::isfinite(j_eq0) && j_eq0 != 0.0) {
                rec_sgn  = (j_eq0 < 0.0) ? 1.0 : -1.0;
                rec_psi0 = line_mean(x, psi_line);
                rec_d2   = std::abs(j_eq0);
                size_t kX = 0;
                for (size_t k = 1; k < x.size(); ++k) {
                    if (rec_sgn * psi_line[k] < rec_sgn * psi_line[kX]) { kX = k; }
                }
                rec_xX = x[kX];
                rec_ref_set = true;
            }
        }
        if (rec_ref_set) {
            // psi at the X-point node (shared panel copies averaged)
            double s = 0.0; int c = 0;
            for (size_t i = 0; i < xs.size(); ++i) {
                if (std::abs(xs[i] - rec_xX) < tol_x && std::abs(ys[i]) < tol_y) { s += psis[i]; ++c; }
            }
            if (c > 0) {
                d.Psi = rec_sgn * (rec_psi0 - s / c);
                d.W   = 4.0 * std::sqrt(std::max(d.Psi, 0.0) / rec_d2);
                // backward difference between consecutive diagnostic dumps,
                // so gamma_eff on this row is centred at t - diag_interval/2
                if (std::isfinite(rec_Psi_prev) && rec_Psi_prev > 0.0 && d.Psi > 0.0 && t > rec_t_prev) {
                    d.gamma_eff = std::log(d.Psi / rec_Psi_prev) / (t - rec_t_prev);
                }
                rec_Psi_prev = d.Psi;
                rec_t_prev   = t;
            }
        }
    }

    // sheet width: FWHM of |j| along the column through the X-point (placeholder, to be revised)
    const double x_col = rec_ref_set ? rec_xX : 0.5 * (x_min + x_max);
    d.delta_cs = 2 * line_hwhm(ys, xs, j0s, x_col, tol_x, tol_y);
    return d;
}

// diagnostics.csv: '#' preamble describing every column, then one header row and
// fixed-width, comma-separated rows grouped as  time | energy | invariants | peaks | reconnection.
// Read with  pandas.read_csv(path, comment='#', skipinitialspace=True).
namespace {
struct DiagColumn { const char* name; const char* meaning; };
const DiagColumn diag_columns[] = {
    {"step",            "time step index"},
    {"time",            "simulation time t (Alfven units)"},
    {"E_kin",           "kinetic energy   0.5 * int |u|^2 dA  (self-consistent flow, external stagnation flow excluded)"},
    {"E_mag",           "magnetic energy  0.5 * int |B|^2 dA"},
    {"E_tot",           "total energy     E_kin + E_mag"},
    {"cross_helicity",  "int u . B dA"},
    {"total_current",   "int j dA      (conserved; nonzero drift = boundary flux or error)"},
    {"total_vorticity", "int omega dA  (conserved)"},
    {"w_max",           "max |omega| over mesh nodes"},
    {"j_max",           "max |j| over mesh nodes"},
    {"Psi",             "reconnected flux, Loureiro et al. 2005: psi_0 - psi(t, X); NaN for free bcs"},
    {"W",               "island width 4 sqrt(Psi / |psi_0''(0)|), Loureiro Eq. (3)"},
    {"gamma_eff",       "effective growth rate d ln Psi / dt, backward difference over one diag interval"},
    {"delta_cs",        "sheet width: FWHM of |j| along the column through the X-point (placeholder)"},
};
const int diag_ncol = sizeof(diag_columns) / sizeof(diag_columns[0]);
const int diag_w = 24;   // column width: fits % .16e
}

int AMRStructure::write_diagnostics(const MHDDiagnostics& d) {
    static bool header_written = false;
    const std::string path = sim_dir + "simulation_output/diagnostics.csv";

    std::ofstream f;
    if (!header_written) {
        f.open(path, std::ios::out | std::ios::trunc);
    } else {
        f.open(path, std::ios::out | std::ios::app);
    }
    if (!f.is_open()) {
        std::cerr << "[diag] failed to open " << path << std::endl;
        return 1;
    }

    if (!header_written) {
        f << "# PMMHD diagnostics, written every " << n_steps_diag << " steps (dt = " << dt << ")\n";
        f << "# box [" << x_min << ", " << x_max << "] x [" << y_min << ", " << y_max << "], bcs = " << bcs
          << ", nu = " << nu << ", eta = " << mu << "\n";
        f << "# columns:\n";
        for (int c = 0; c < diag_ncol; ++c) {
            f << "#   " << std::left << std::setw(16) << diag_columns[c].name << " " << diag_columns[c].meaning << "\n";
        }
        f << std::right;
        for (int c = 0; c < diag_ncol; ++c) {
            f << (c ? "," : "") << std::setw(c ? diag_w : 8) << diag_columns[c].name;
        }
        f << "\n";
        header_written = true;
    }

    const double vals[] = {d.t, d.E_kin, d.E_mag, d.E_tot, d.H_C, d.I_j, d.I_w, d.w_max, d.j_max,
                           d.Psi, d.W, d.gamma_eff, d.delta_cs};
    f << std::setw(8) << d.iter;
    f << std::scientific << std::setprecision(16);
    for (double v : vals) { f << "," << std::setw(diag_w) << v; }
    f << "\n";
    return 0;
}