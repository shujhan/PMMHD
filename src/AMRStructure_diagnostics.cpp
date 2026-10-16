#include "AMRStructure.hpp"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

// Midplane flux variation max psi - min psi on y = 0. Used only at t = 0 to fix
// the reference of the X-point reconnected flux (see compute_diagnostics).
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

// Flux function at the X-point (x_c, 0) relative to the x-mean of psi on the
// boundary row y = y_min:
//     D = psi(x_c, 0) - <psi(., y_min)>_x
//       = [psi(x_c,0) - psi(x_c,y_min)] - [<psi(.,y_min)>_x - psi(x_c,y_min)]
// from line integrals of B (d_y psi = -b1 up the column, d_x psi = b2 along
// the row), so it needs no absolute psi. The seed cos(kx) averages out of the
// row mean, so the reference is x-independent like the gauge of Loureiro's
// spectral psi. Returns NaN if either line is missing.
static double xpoint_flux(const std::vector<double>& xs,
                          const std::vector<double>& ys,
                          const std::vector<double>& b1s,
                          const std::vector<double>& b2s,
                          double x_c, double y_min, double Lx, double Ly)
{
    const double tol_x = 1e-9 * Lx, tol_y = 1e-9 * Ly;

    // column x = x_c, from y_min up to y = 0
    std::vector<double> y, b1;
    mesh_line(ys, xs, b1s, x_c, tol_x, tol_y, y, b1);
    if (y.size() < 3 || std::abs(y.front() - y_min) > tol_y) { return std::nan(""); }
    double col = 0.0;
    bool hit0 = false;
    for (size_t k = 1; k < y.size(); ++k) {
        if (y[k] > tol_y) { break; }
        col -= 0.5 * (b1[k] + b1[k-1]) * (y[k] - y[k-1]);
        if (std::abs(y[k]) < tol_y) { hit0 = true; break; }
    }
    if (!hit0) { return std::nan(""); }

    // row y = y_min: psi(x) - psi(x_min) by cumulative trapezoid of b2
    std::vector<double> x, b2;
    mesh_line(xs, ys, b2s, y_min, tol_y, tol_x, x, b2);
    const size_t M = x.size();
    if (M < 3) { return std::nan(""); }
    std::vector<double> p(M, 0.0);
    for (size_t k = 1; k < M; ++k) {
        p[k] = p[k-1] + 0.5 * (b2[k] + b2[k-1]) * (x[k] - x[k-1]);
    }
    double area = 0.0, p_c = std::nan("");
    for (size_t k = 0; k < M; ++k) {
        if (k > 0) { area += 0.5 * (p[k] + p[k-1]) * (x[k] - x[k-1]); }
        if (std::abs(x[k] - x_c) < tol_x) { p_c = p[k]; }
    }
    const double p_mean = area / (x[M-1] - x[0]);

    return col - (p_mean - p_c);
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

    const double x_c = 0.5 * (x_min + x_max);

    if (bcs == free_bcs) {
        // midplane_flux assumes periodic x; not meaningful on a free line
        d.Psi_rec = std::nan("");
        d.Psi_res = std::nan("");
        d.W       = std::nan("");
    } else {
        double Psi_mid = 0.0;
        midplane_flux(xs, ys, b2s, Lx, Ly, Psi_mid, d.Psi_res);

        // Reconnected flux as in Loureiro et al. 2005: Psi(t) = psi_eq(X) - psi(X, t)
        // at the X-point (x_c, 0). With the gauge of xpoint_flux, psi_eq(X) is
        // fixed once at t = 0 from
        //     Psi(0) = |seed flux at X| = (max - min of psi on y = 0) / 2,
        // exact for a separable seed j_0(y) (1 + amp cos kx): psi(x, 0) on the
        // midplane is psi_0(0) + A cos kx, so max - min = 2|A|.
        // No equilibrium-diffusion subtraction, as in the paper: Psi also
        // grows at ~ eta |j(X)| from diffusion of the sheet itself.
        const double D = xpoint_flux(xs, ys, b1s, b2s, x_c, y_min, Lx, Ly);
        if (!psiL_ref_set && std::isfinite(D)) {
            double jX = std::nan("");
            for (size_t i = 0; i < N; ++i) {
                if (std::abs(xs[i] - x_c) < 1e-9 * Lx && std::abs(ys[i]) < 1e-9 * Ly) {
                    jX = j0s[i];
                    break;
                }
            }
            if (std::isfinite(jX) && jX != 0.0) {
                psiL_sgn = (jX < 0.0) ? 1.0 : -1.0;
                psiL_jX0 = std::abs(jX);
                psiL_ref = D + psiL_sgn * 0.5 * Psi_mid;
                psiL_ref_set = true;
            }
        }
        if (psiL_ref_set && std::isfinite(D)) {
            d.Psi_rec = psiL_sgn * (psiL_ref - D);
            d.W       = 4.0 * std::sqrt(std::max(d.Psi_rec, 0.0) / psiL_jX0);
        } else {
            d.Psi_rec = std::nan("");
            d.W       = std::nan("");
        }
    }

    // sheet size from |j| along the center column (through the X-point) and y = 0
    const double fwhm_x = 2 * line_hwhm(ys, xs, j0s, x_c, 1e-9 * Lx, 1e-9 * Ly);
    // d.sheet_a  = fwhm_x / std::acosh(std::sqrt(2.0));
    d.delta_cs = fwhm_x;
    // d.sheet_b  = 2 * line_hwhm(xs, ys, j0s, 0.0, 1e-9 * Ly, 1e-9 * Lx);
    return d;
}

int AMRStructure::write_diagnostics(const MHDDiagnostics& d) {
    static bool header_written = false;
    const std::string path = sim_dir + "simulation_output/diagnostics.csv";

    std::ofstream f;
    if (!header_written) {
        f.open(path, std::ios::out | std::ios::trunc);
        f << "iter,t,E_kin,E_mag,E_tot,H_C,I_j,I_w,Psi_rec,Psi_res,w_max,j_max,sheet_width,sheet_length,W,delta_cs\n";
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
    //   << d.sheet_a << "," << d.sheet_b << ","
      << d.W << "," << d.delta_cs << "\n";
    return 0;
}