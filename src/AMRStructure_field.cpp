#include "AMRStructure.hpp"

// int AMRStructure::evaluate_u_field(std::vector<double>& u1s_local, std::vector<double>& u2s_local, 
//                      std::vector<double>& xs_local,std::vector<double>& ys_local,
//                      std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = original;
//     calculate_e->set_mode(m);
//     (*calculate_e)(u1s_local.data(), u2s_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }



// int AMRStructure::evaluate_b_field(std::vector<double>& b1s_local, std::vector<double>& b2s_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = original;
//     calculate_e->set_mode(m);
//     (*calculate_e)(b1s_local.data(), b2s_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }

int AMRStructure::evaluate_u_field(std::vector<double>& u1s_local, std::vector<double>& u2s_local,
                                   std::vector<double>& xs_local, std::vector<double>& ys_local,
                                   std::vector<double>& ws_local, double t)
{
    const int n_local = (int)xs_local.size();

    if (bcs == periodic_bcs) {
        // ---------- Step A: build 9x image list, 3x3 near sum at PARTICLES ----------
        const int n9 = 9 * n_local;
        std::vector<double> xtmp(n9), ytmp(n9), wtmp(n9);
        const double sx[9] = {0.0, -Lx,  Lx, 0.0, 0.0, -Lx, -Lx,  Lx,  Lx};
        const double sy[9] = {0.0,  0.0, 0.0, -Ly, Ly, -Ly,  Ly, -Ly,  Ly};
        for (int img = 0; img < 9; ++img) {
            for (int i = 0; i < n_local; ++i) {
                xtmp[img * n_local + i] = xs_local[i] + sx[img];
                ytmp[img * n_local + i] = ys_local[i] + sy[img];
                wtmp[img * n_local + i] = ws_local[i];
            }
        }

        std::vector<double> u1tmp(n9, 0.0), u2tmp(n9, 0.0);
        KernelMode km = periodic_xy;   // free-space 2D Biot-Savart kernel branch
        calculate_e->set_mode(km);
        // sources: all n9 image points; targets: the first n_local (center block)
        (*calculate_e)(u1tmp.data(), u2tmp.data(),
                       xtmp.data(), n_local,
                       ytmp.data(), wtmp.data(), n9);

        for (int i = 0; i < n_local; ++i) {
            u1s_local[i] = u1tmp[i];
            u2s_local[i] = u2tmp[i];
        }

        // ---------- Step B: same near sum at the 4m wall collocation targets ----------
        const std::vector<double>& wxs = periodizer->get_wall_xs();
        const std::vector<double>& wys = periodizer->get_wall_ys();
        const int mm = periodizer->get_m();
        const int nw = 4 * mm;

        std::vector<double> near_walls_u1(nw, 0.0);
        std::vector<double> near_walls_u2(nw, 0.0);

        const double pi = 3.14159265358979323846;
        const double eps = greens_epsilon;

        #pragma omp parallel for
        for (int i = 0; i < nw; ++i) {
            double tx = wxs[i];
            double ty = wys[i];
            double s1 = 0.0, s2 = 0.0;
            for (int k = 0; k < n9; ++k) {
                double dx = tx - xtmp[k];
                double dy = ty - ytmp[k];
                // was: if (x_diff == 0.0 && y_diff == 0.0) continue;
                double r2 = dx * dx + dy * dy + eps * eps;
                if (r2 < 1e-16) continue;   // skip self + coincident periodic images
                s1 -= (1.0 / (2.0 * pi)) * dy / r2 * wtmp[k];
                s2 += (1.0 / (2.0 * pi)) * dx / r2 * wtmp[k];
            }
            near_walls_u1[i] = s1;
            near_walls_u2[i] = s2;
        }

        // ---------- Step C: pack 8m vector for Periodizer ----------
        std::vector<double> near_u_on_walls(8 * mm);
        for (int i = 0; i < mm; ++i) {
            near_u_on_walls[          i] = near_walls_u1[         i];  // L u1
            near_u_on_walls[   mm  + i] = near_walls_u2[         i];  // L u2
            near_u_on_walls[ 2*mm  + i] = near_walls_u1[  mm   + i];  // D u1
            near_u_on_walls[ 3*mm  + i] = near_walls_u2[  mm   + i];  // D u2
            near_u_on_walls[ 4*mm  + i] = near_walls_u1[2*mm   + i];  // R u1
            near_u_on_walls[ 5*mm  + i] = near_walls_u2[2*mm   + i];  // R u2
            near_u_on_walls[ 6*mm  + i] = near_walls_u1[3*mm   + i];  // U u1
            near_u_on_walls[ 7*mm  + i] = near_walls_u2[3*mm   + i];  // U u2
        }

        // ---------- Step D: solve Q xi = g for proxy coefficients ----------
        periodizer->compute_xi(near_u_on_walls);

        // ---------- Step E: add proxy correction at the particles ----------
        periodizer->add_correction(xs_local, ys_local, u1s_local, u2s_local);

        // ---------- Step F: project onto the zero-mean gauge ----------
        // The doubly-periodic Biot-Savart velocity of a zero-mean vorticity is
        // defined only up to a uniform constant. The 3x3-image near sum carries a
        // spurious mean (truncation of the infinite periodic sum) that the
        // jump-only wall matching does not remove, leaving a phantom mean flow
        // that inflates E_kin and uniformly advects the field. Subtract the
        // quadrature-weighted domain mean to restore the zero-mean gauge.
        {
            const bool use_w = ((int)weights.size() == n_local);
            double m1 = 0.0, m2 = 0.0, wsum = 0.0;
            for (int i = 0; i < n_local; ++i) {
                const double wq = use_w ? weights[i] : 1.0;
                m1 += wq * u1s_local[i];
                m2 += wq * u2s_local[i];
                wsum += wq;
            }
            if (wsum != 0.0) {
                m1 /= wsum; m2 /= wsum;
                for (int i = 0; i < n_local; ++i) {
                    u1s_local[i] -= m1;
                    u2s_local[i] -= m2;
                }
            }
        }
    }
    else if (bcs == free_bcs) {
        // free case (open x and y): free-space kernel, sources = targets = particles.
        // No images, no Periodizer, no zero-mean projection (the free-space field
        // is unique: it decays at infinity).
        std::vector<double> xtmp(xs_local), ytmp(ys_local), wtmp(ws_local);
        std::vector<double> u1tmp(n_local, 0.0), u2tmp(n_local, 0.0);
        KernelMode km = free_xy;
        calculate_e->set_mode(km);
        (*calculate_e)(u1tmp.data(), u2tmp.data(),
                       xtmp.data(), n_local,
                       ytmp.data(), wtmp.data(), n_local);
        u1s_local = u1tmp;
        u2s_local = u2tmp;
    }
    else {
        // channel case (periodic x, open y): unchanged, use analytic channel kernel
        std::vector<double> u1tmp(n_local, 0.0), u2tmp(n_local, 0.0);
        KernelMode km = original;
        calculate_e->set_mode(km);
        (*calculate_e)(u1tmp.data(), u2tmp.data(),
                       xs_local.data(), n_local,
                       ys_local.data(), ws_local.data(), n_local);
        u1s_local = u1tmp;
        u2s_local = u2tmp;
    }
    return 0;
}


int AMRStructure::evaluate_b_field(std::vector<double>& b1s_local, std::vector<double>& b2s_local,
                                   std::vector<double>& xs_local, std::vector<double>& ys_local,
                                   std::vector<double>& ws_local, double t)
{
    const int n_local = (int)xs_local.size();

    if (bcs == periodic_bcs) {
        // ---------- Step A: build 9x image list, 3x3 near sum at PARTICLES ----------
        const int n9 = 9 * n_local;
        std::vector<double> xtmp(n9), ytmp(n9), wtmp(n9);
        const double sx[9] = {0.0, -Lx,  Lx, 0.0, 0.0, -Lx, -Lx,  Lx,  Lx};
        const double sy[9] = {0.0,  0.0, 0.0, -Ly, Ly, -Ly,  Ly, -Ly,  Ly};
        for (int img = 0; img < 9; ++img) {
            for (int i = 0; i < n_local; ++i) {
                xtmp[img * n_local + i] = xs_local[i] + sx[img];
                ytmp[img * n_local + i] = ys_local[i] + sy[img];
                wtmp[img * n_local + i] = ws_local[i];
            }
        }

        std::vector<double> b1tmp(n9, 0.0), b2tmp(n9, 0.0);
        KernelMode km = periodic_xy;
        calculate_e->set_mode(km);
        // sources: all n9 image points; targets: the first n_local (center block)
        (*calculate_e)(b1tmp.data(), b2tmp.data(),
                       xtmp.data(), n_local,
                       ytmp.data(), wtmp.data(), n9);

        for (int i = 0; i < n_local; ++i) {
            b1s_local[i] = b1tmp[i];
            b2s_local[i] = b2tmp[i];
        }

        // ---------- Step B: same near sum at the 4m wall collocation targets ----------
        const std::vector<double>& wxs = periodizer->get_wall_xs();
        const std::vector<double>& wys = periodizer->get_wall_ys();
        const int mm = periodizer->get_m();
        const int nw = 4 * mm;

        std::vector<double> near_walls_b1(nw, 0.0);
        std::vector<double> near_walls_b2(nw, 0.0);

        const double pi = 3.14159265358979323846;
        const double eps = greens_epsilon;

        #pragma omp parallel for
        for (int i = 0; i < nw; ++i) {
            double tx = wxs[i];
            double ty = wys[i];
            double s1 = 0.0, s2 = 0.0;
            for (int k = 0; k < n9; ++k) {
                double dx = tx - xtmp[k];
                double dy = ty - ytmp[k];
                // was: if (x_diff == 0.0 && y_diff == 0.0) continue;
                double r2 = dx * dx + dy * dy + eps * eps;
                if (r2 < 1e-16) continue;
                s1 -= (1.0 / (2.0 * pi)) * dy / r2 * wtmp[k];
                s2 += (1.0 / (2.0 * pi)) * dx / r2 * wtmp[k];
            }
            near_walls_b1[i] = s1;
            near_walls_b2[i] = s2;
        }

        // ---------- Step C: pack 8m vector for Periodizer ----------
        std::vector<double> near_b_on_walls(8 * mm);
        for (int i = 0; i < mm; ++i) {
            near_b_on_walls[          i] = near_walls_b1[         i];  // L b1
            near_b_on_walls[   mm  + i] = near_walls_b2[         i];  // L b2
            near_b_on_walls[ 2*mm  + i] = near_walls_b1[  mm   + i];  // D b1
            near_b_on_walls[ 3*mm  + i] = near_walls_b2[  mm   + i];  // D b2
            near_b_on_walls[ 4*mm  + i] = near_walls_b1[2*mm   + i];  // R b1
            near_b_on_walls[ 5*mm  + i] = near_walls_b2[2*mm   + i];  // R b2
            near_b_on_walls[ 6*mm  + i] = near_walls_b1[3*mm   + i];  // U b1
            near_b_on_walls[ 7*mm  + i] = near_walls_b2[3*mm   + i];  // U b2
        }

        // ---------- Step D: solve Q xi = g for proxy coefficients (for B) ----------
        periodizer->compute_xi(near_b_on_walls);

        // ---------- Step E: add proxy correction at the particles ----------
        periodizer->add_correction(xs_local, ys_local, b1s_local, b2s_local);

        // ---------- Step F: project onto the zero-mean gauge ----------
        // Same fix as in evaluate_u_field: remove the spurious mean left by the
        // near-sum truncation. This zeroes the mean of the Biot-Savart part of B;
        // the physical background field B0x/B0y is added by the caller afterward,
        // so the total mean B is exactly (B0x, B0y) as intended.
        {
            const bool use_w = ((int)weights.size() == n_local);
            double m1 = 0.0, m2 = 0.0, wsum = 0.0;
            for (int i = 0; i < n_local; ++i) {
                const double wq = use_w ? weights[i] : 1.0;
                m1 += wq * b1s_local[i];
                m2 += wq * b2s_local[i];
                wsum += wq;
            }
            if (wsum != 0.0) {
                m1 /= wsum; m2 /= wsum;
                for (int i = 0; i < n_local; ++i) {
                    b1s_local[i] -= m1;
                    b2s_local[i] -= m2;
                }
            }
        }
    }
    else if (bcs == free_bcs) {
        // free case (open x and y): free-space kernel, no images / Periodizer / gauge
        std::vector<double> xtmp(xs_local), ytmp(ys_local), wtmp(ws_local);
        std::vector<double> b1tmp(n_local, 0.0), b2tmp(n_local, 0.0);
        KernelMode km = free_xy;
        calculate_e->set_mode(km);
        (*calculate_e)(b1tmp.data(), b2tmp.data(),
                       xtmp.data(), n_local,
                       ytmp.data(), wtmp.data(), n_local);
        b1s_local = b1tmp;
        b2s_local = b2tmp;
    }
    else {
        // channel case (periodic x, open y): unchanged
        std::vector<double> b1tmp(n_local, 0.0), b2tmp(n_local, 0.0);
        KernelMode km = original;
        calculate_e->set_mode(km);
        (*calculate_e)(b1tmp.data(), b2tmp.data(),
                       xs_local.data(), n_local,
                       ys_local.data(), ws_local.data(), n_local);
        b1s_local = b1tmp;
        b2s_local = b2tmp;
    }
    return 0;
}



int AMRStructure::evaluate_potential(std::vector<double>& potentials_local,
                                   std::vector<double>& xs_local, std::vector<double>& ys_local,
                                   std::vector<double>& ws_local, double t)
{
    const int n_local = (int)xs_local.size();
    potentials_local.assign(n_local, 0.0);
    if (n_local == 0) { return 0; }

    if (bcs == periodic_bcs) {
        // ---------- Step A: build 9x image list, 3x3 near sum at PARTICLES ----------
        const int n9 = 9 * n_local;
        std::vector<double> xtmp(n9), ytmp(n9), wtmp(n9);
        const double sx[9] = {0.0, -Lx,  Lx, 0.0, 0.0, -Lx, -Lx,  Lx,  Lx};
        const double sy[9] = {0.0,  0.0, 0.0, -Ly, Ly, -Ly,  Ly, -Ly,  Ly};
        for (int img = 0; img < 9; ++img) {
            for (int i = 0; i < n_local; ++i) {
                xtmp[img * n_local + i] = xs_local[i] + sx[img];
                ytmp[img * n_local + i] = ys_local[i] + sy[img];
                wtmp[img * n_local + i] = ws_local[i];
            }
        }

        std::vector<double> potentialtmp(n9, 0.0), b2tmp(n9, 0.0);
        KernelMode km = periodic_xy_potentials;
        calculate_e->set_mode(km);
        // sources: all n9 image points; targets: the first n_local (center block)
        (*calculate_e)(potentialtmp.data(), b2tmp.data(),
                       xtmp.data(), n_local,
                       ytmp.data(), wtmp.data(), n9);

        for (int i = 0; i < n_local; ++i) {
            potentials_local[i] = potentialtmp[i];
        }

        // ---------- Step B: near sum at the 4m wall collocation targets ----------
        // Field (for the proxy fit, as in evaluate_b_field) and potential (for
        // the linear-gauge fix in Step F), from the same 3x3 image list.
        const std::vector<double>& wxs = periodizer->get_wall_xs();
        const std::vector<double>& wys = periodizer->get_wall_ys();
        const int mm = periodizer->get_m();
        const int nw = 4 * mm;

        std::vector<double> near_walls_f1(nw, 0.0);
        std::vector<double> near_walls_f2(nw, 0.0);
        std::vector<double> walls_pot(nw, 0.0);

        const double pi = 3.14159265358979323846;
        const double eps = greens_epsilon;

        #pragma omp parallel for
        for (int i = 0; i < nw; ++i) {
            double tx = wxs[i];
            double ty = wys[i];
            double s1 = 0.0, s2 = 0.0, sp = 0.0;
            for (int k = 0; k < n9; ++k) {
                double dx = tx - xtmp[k];
                double dy = ty - ytmp[k];
                double r2 = dx * dx + dy * dy + eps * eps;
                if (r2 < 1e-16) continue;   // skip coincident points, as in the field solve
                s1 -= (1.0 / (2.0 * pi)) * dy / r2 * wtmp[k];
                s2 += (1.0 / (2.0 * pi)) * dx / r2 * wtmp[k];
                sp += (0.25 / pi) * std::log(r2) * wtmp[k];
            }
            near_walls_f1[i] = s1;
            near_walls_f2[i] = s2;
            walls_pot[i] = sp;
        }

        // ---------- Step C: pack 8m vector for Periodizer ----------
        std::vector<double> near_f_on_walls(8 * mm);
        for (int i = 0; i < mm; ++i) {
            near_f_on_walls[          i] = near_walls_f1[         i];  // L f1
            near_f_on_walls[   mm  + i] = near_walls_f2[         i];  // L f2
            near_f_on_walls[ 2*mm  + i] = near_walls_f1[  mm   + i];  // D f1
            near_f_on_walls[ 3*mm  + i] = near_walls_f2[  mm   + i];  // D f2
            near_f_on_walls[ 4*mm  + i] = near_walls_f1[2*mm   + i];  // R f1
            near_f_on_walls[ 5*mm  + i] = near_walls_f2[2*mm   + i];  // R f2
            near_f_on_walls[ 6*mm  + i] = near_walls_f1[3*mm   + i];  // U f1
            near_f_on_walls[ 7*mm  + i] = near_walls_f2[3*mm   + i];  // U f2
        }

        // ---------- Step D: proxy coefficients from the field jumps ----------
        // Same fit as evaluate_u_field / evaluate_b_field for these weights.
        // Overwrites the periodizer's xi; every field solve recomputes its own.
        periodizer->compute_xi(near_f_on_walls);

        // ---------- Step E: proxy potential at the particles and the walls ----------
        periodizer->add_correction_potential(xs_local, ys_local, potentials_local);
        periodizer->add_correction_potential(wxs, wys, walls_pot);

        // ---------- Step F: remove the linear part, then the mean ----------
        // With the field jumps cancelled, grad(potential) is periodic, so the
        // potential is periodic up to c_x (x-xc)/Lx + c_y (y-yc)/Ly with
        // c_x = <pot(R) - pot(L)>, c_y = <pot(U) - pot(D)> over the wall pairs.
        // Removing it is the zero-mean-field gauge of the field solve (Step F
        // there). The constant is then fixed by a zero quadrature-weighted mean,
        // as in a spectral solve (k = 0 mode dropped).
        {
            double c_x = 0.0, c_y = 0.0;
            for (int i = 0; i < mm; ++i) {
                c_x += walls_pot[2*mm + i] - walls_pot[       i];   // R - L
                c_y += walls_pot[3*mm + i] - walls_pot[  mm + i];   // U - D
            }
            c_x /= mm; c_y /= mm;
            const double xc = 0.5 * (x_min + x_max), yc = 0.5 * (y_min + y_max);

            const bool use_w = ((int)weights.size() == n_local);
            double mean = 0.0, wsum = 0.0;
            for (int i = 0; i < n_local; ++i) {
                potentials_local[i] -= c_x * (xs_local[i] - xc) / Lx
                                     + c_y * (ys_local[i] - yc) / Ly;
                const double wq = use_w ? weights[i] : 1.0;
                mean += wq * potentials_local[i];
                wsum += wq;
            }
            if (wsum != 0.0) {
                mean /= wsum;
                for (int i = 0; i < n_local; ++i) { potentials_local[i] -= mean; }
            }
        }
    }
    else if (bcs == free_bcs) {
        // free case (open x and y): free-space log kernel, no images
        std::vector<double> xtmp(xs_local), ytmp(ys_local), wtmp(ws_local);
        std::vector<double> potentialtmp(n_local, 0.0), b2tmp(n_local, 0.0);
        KernelMode km = free_xy_potentials;
        calculate_e->set_mode(km);
        (*calculate_e)(potentialtmp.data(), b2tmp.data(),
                       xtmp.data(), n_local,
                       ytmp.data(), wtmp.data(), n_local);
        potentials_local = potentialtmp;
    }
    else {
        // channel case (periodic x, open y): unchanged
        std::vector<double> potentialtmp(n_local, 0.0), b2tmp(n_local, 0.0);
        KernelMode km = original_potentials;
        calculate_e->set_mode(km);
        (*calculate_e)(potentialtmp.data(), b2tmp.data(),
                       xs_local.data(), n_local,
                       ys_local.data(), ws_local.data(), n_local);
        potentials_local = potentialtmp;
    }
    return 0;
}























// int AMRStructure::evaluate_u1s_grad(std::vector<double>& u1s_grad_x_local, std::vector<double>& u1s_grad_y_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = u1_grad;
//     calculate_e->set_mode(m);
//     (*calculate_e)(u1s_grad_x_local.data(), u1s_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }

// int AMRStructure::evaluate_u2s_grad(std::vector<double>& u2s_grad_x_local, std::vector<double>& u2s_grad_y_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = u2_grad;
//     calculate_e->set_mode(m);
//     (*calculate_e)(u2s_grad_x_local.data(), u2s_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }

// int AMRStructure::evaluate_b1s_grad(std::vector<double>& b1s_grad_x_local, std::vector<double>& b1s_grad_y_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = u1_grad;
//     calculate_e->set_mode(m);
//     (*calculate_e)(b1s_grad_x_local.data(), b1s_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }

// int AMRStructure::evaluate_b2s_grad(std::vector<double>& b2s_grad_x_local, std::vector<double>& b2s_grad_y_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = u2_grad;
//     calculate_e->set_mode(m);
//     (*calculate_e)(b2s_grad_x_local.data(), b2s_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }


// int AMRStructure::evaluate_vorticity_grad(std::vector<double>& vorticity_grad_x_local, std::vector<double>& vorticity_grad_y_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = vorticity_grad;
//     calculate_e->set_mode(m);
//     (*calculate_e)(vorticity_grad_x_local.data(), vorticity_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }


// int AMRStructure::evaluate_j_grad(std::vector<double>& j_grad_x_local, std::vector<double>& j_grad_y_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = vorticity_grad;
//     calculate_e->set_mode(m);
//     (*calculate_e)(j_grad_x_local.data(), j_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }


// int AMRStructure::evaluate_vorticity_laplacian(std::vector<double>& vorticity_laplacian_local, std::vector<double>& vorticity_none_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = laplacian;
//     calculate_e->set_mode(m);
//     (*calculate_e)(vorticity_laplacian_local.data(), vorticity_none_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }

// int AMRStructure::evaluate_j_laplacian(std::vector<double>& j_laplacian_local, std::vector<double>& j_none_local, 
//                     std::vector<double>& xs_local,std::vector<double>& ys_local,
//                         std::vector<double>& ws_local, double t) {
//     std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
//     std::vector<double> ws_temp_cpy(ws_local);
//     KernelMode m = laplacian;
//     calculate_e->set_mode(m);
//     (*calculate_e)(j_laplacian_local.data(),j_none_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
//                     ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
//     return 0;
// }