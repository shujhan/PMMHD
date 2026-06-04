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
        (*calculate_e)(u1tmp.data(), u2tmp.data(),
                       xtmp.data(), n9,
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
                double r2 = dx * dx + dy * dy + eps * eps;
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
        (*calculate_e)(b1tmp.data(), b2tmp.data(),
                       xtmp.data(), n9,
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
                double r2 = dx * dx + dy * dy + eps * eps;
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

