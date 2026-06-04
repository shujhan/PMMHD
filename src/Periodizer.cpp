#include "Periodizer.hpp"
#include <cmath>
#include <cassert>
#include <iostream>

// -----------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------
Periodizer::Periodizer(double x_min, double x_max,
                       double y_min, double y_max,
                       double eps,
                       Field* free_space_kernel,
                       int M_, int m_, double Rp_factor)
    : Lx(x_max - x_min), Ly(y_max - y_min), eps(eps),
      cx(0.5 * (x_min + x_max)), cy(0.5 * (y_min + y_max)),
      M(M_), m(m_), kernel(free_space_kernel)
{
    // Proxy circle radius: must sit outside the unit cell but inside the
    // nearest shell-2 image.  r_B = half-diagonal of the cell.
    const double rB = 0.5 * std::sqrt(Lx * Lx + Ly * Ly);
    Rp = Rp_factor * rB;

    // ---- proxy points on the circle ----
    proxy_x.resize(M);
    proxy_y.resize(M);
    for (int j = 0; j < M; ++j) {
        double th = 2.0 * M_PI * j / M;
        proxy_x[j] = cx + Rp * std::cos(th);
        proxy_y[j] = cy + Rp * std::sin(th);
    }

    // ---- wall collocation nodes ----
    // Layout of wall_xs_4m / wall_ys_4m (each of size 4m):
    //   [0   .. m)    L wall:  x = x_min, y at m nodes on [y_min, y_max]
    //   [m   .. 2m)   D wall:  y = y_min, x at m nodes on [x_min, x_max]
    //   [2m  .. 3m)   R wall:  x = x_max, same y nodes as L  (= L + e1)
    //   [3m  .. 4m)   U wall:  y = y_max, same x nodes as D  (= D + e2)
    wall_xs_4m.resize(4 * m);
    wall_ys_4m.resize(4 * m);

    for (int i = 0; i < m; ++i) {
        // Gauss-Legendre nodes would be ideal; uniform midpoints used here.
        // Replace s_y / s_x with your GL nodes for slightly higher accuracy.
        double s_y = (i + 0.5) / m;   // parameter in [0,1] along L/R walls
        double s_x = (i + 0.5) / m;   // parameter in [0,1] along D/U walls

        // L wall
        wall_xs_4m[i]         = x_min;
        wall_ys_4m[i]         = y_min + s_y * Ly;
        // D wall
        wall_xs_4m[m + i]     = x_min + s_x * Lx;
        wall_ys_4m[m + i]     = y_min;
        // R wall  (= L + e1)
        wall_xs_4m[2 * m + i] = x_max;
        wall_ys_4m[2 * m + i] = y_min + s_y * Ly;
        // U wall  (= D + e2)
        wall_xs_4m[3 * m + i] = x_min + s_x * Lx;
        wall_ys_4m[3 * m + i] = y_max;
    }

    xi.assign(M, 0.0);
    // Q_matrix and Q_solver are default-constructed; filled in precompute_Q()
}

// -----------------------------------------------------------------------
// precompute_Q  — fill Q and factorize it.  Call once after construction.
// -----------------------------------------------------------------------
void Periodizer::precompute_Q()
{
    // Q has 4m rows and M columns.
    // Row layout:
    //   [0   .. m)   u1 jump at L/R pair:  u1(R_i) - u1(L_i)
    //   [m   .. 2m)  u2 jump at L/R pair:  u2(R_i) - u2(L_i)
    //   [2m  .. 3m)  u1 jump at D/U pair:  u1(U_i) - u1(D_i)
    //   [3m  .. 4m)  u2 jump at D/U pair:  u2(U_i) - u2(D_i)

    Q_matrix.resize(4 * m, M);

    const double pi  = 3.14159265358979323846;
    // const double eps = 1e-12;   // proxy is far from all walls, no singularity
    // const double eps = 0.1;

    for (int j = 0; j < M; ++j) {
        double yx = proxy_x[j];
        double yy = proxy_y[j];

        // velocity at target (tx, ty) due to unit-strength source at (yx, yy)
        auto biot = [&](double tx, double ty, double& u1, double& u2) {
            double dx = tx - yx;
            double dy = ty - yy;
            double r2 = dx * dx + dy * dy + eps * eps;
            u1 = -(1.0 / (2.0 * pi)) * dy / r2;
            u2 =  (1.0 / (2.0 * pi)) * dx / r2;
        };

        for (int i = 0; i < m; ++i) {
            double u1L, u2L, u1R, u2R;
            biot(wall_xs_4m[i],         wall_ys_4m[i],         u1L, u2L);  // L
            biot(wall_xs_4m[2 * m + i], wall_ys_4m[2 * m + i], u1R, u2R);  // R

            Q_matrix(i,       j) = u1R - u1L;   // u1 jump L->R
            Q_matrix(m + i,   j) = u2R - u2L;   // u2 jump L->R

            double u1D, u2D, u1U, u2U;
            biot(wall_xs_4m[m + i],     wall_ys_4m[m + i],     u1D, u2D);  // D
            biot(wall_xs_4m[3 * m + i], wall_ys_4m[3 * m + i], u1U, u2U);  // U

            Q_matrix(2 * m + i, j) = u1U - u1D;  // u1 jump D->U
            Q_matrix(3 * m + i, j) = u2U - u2D;  // u2 jump D->U
        }
    }

    // Pivoted QR — handles the expected ill-conditioning in Q correctly.
    // Eigen's ColPivHouseholderQR is the direct analogue of LAPACK's dgeqp3
    // and gives a min-norm LS solution when the system is rank-deficient.
    Q_solver.compute(Q_matrix);

    std::cout << "[Periodizer] Q matrix rank: "
              << Q_solver.rank() << " / " << M
              << "  (4m=" << 4*m << " rows, M=" << M << " cols)\n";
}

// -----------------------------------------------------------------------
// compute_xi
// -----------------------------------------------------------------------
// near_u_on_walls has size 8m, packed as:
//   [0   .. m)   u1 at L wall
//   [m   .. 2m)  u2 at L wall
//   [2m  .. 3m)  u1 at D wall
//   [3m  .. 4m)  u2 at D wall
//   [4m  .. 5m)  u1 at R wall
//   [5m  .. 6m)  u2 at R wall
//   [6m  .. 7m)  u1 at U wall
//   [7m  .. 8m)  u2 at U wall
void Periodizer::compute_xi(const std::vector<double>& near_u_on_walls)
{
    assert((int)near_u_on_walls.size() == 8 * m);

    // Build rhs g = -(near_jump), so that Q*xi + near_jump = 0  =>  periodic
    Eigen::VectorXd g(4 * m);
    for (int i = 0; i < m; ++i) {
        double near_u1R = near_u_on_walls[4 * m + i];
        double near_u2R = near_u_on_walls[5 * m + i];
        double near_u1L = near_u_on_walls[      i];
        double near_u2L = near_u_on_walls[  m + i];

        double near_u1U = near_u_on_walls[6 * m + i];
        double near_u2U = near_u_on_walls[7 * m + i];
        double near_u1D = near_u_on_walls[2 * m + i];
        double near_u2D = near_u_on_walls[3 * m + i];

        g(i)         = -(near_u1R - near_u1L);   // u1: cancel R-L jump
        g(m + i)     = -(near_u2R - near_u2L);   // u2: cancel R-L jump
        g(2 * m + i) = -(near_u1U - near_u1D);   // u1: cancel U-D jump
        g(3 * m + i) = -(near_u2U - near_u2D);   // u2: cancel U-D jump
    }

    // Solve Q * xi = g in least-squares sense
    Eigen::VectorXd xi_eigen = Q_solver.solve(g);
    
    double g_norm = g.norm();
    double xi_norm = xi_eigen.norm();
    double residual = (Q_matrix * xi_eigen - g).norm();
    std::cout << "[Periodizer] |g|=" << g_norm 
            << " |xi|=" << xi_norm 
            << " residual=" << residual
            << " rel_residual=" << (g_norm > 0 ? residual/g_norm : 0) << std::endl;

    for (int j = 0; j < M; ++j) xi[j] = xi_eigen(j);
}

// -----------------------------------------------------------------------
// add_correction
// -----------------------------------------------------------------------
void Periodizer::add_correction(const std::vector<double>& tx,
                                const std::vector<double>& ty,
                                std::vector<double>& u1,
                                std::vector<double>& u2) const
{
    const double pi  = 3.14159265358979323846;
    // const double eps = 1e-12;
    // const double eps = 0.1;
    const int N = (int)tx.size();

    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        double s1 = 0.0, s2 = 0.0;
        for (int j = 0; j < M; ++j) {
            double dx = tx[i] - proxy_x[j];
            double dy = ty[i] - proxy_y[j];
            double r2 = dx * dx + dy * dy + eps * eps;
            s1 += xi[j] * (-(1.0 / (2.0 * pi)) * dy / r2);
            s2 += xi[j] * ( (1.0 / (2.0 * pi)) * dx / r2);
        }
        u1[i] += s1;
        u2[i] += s2;
    }
}

// void Periodizer::add_correction_grad_u1(const std::vector<double>& tx,
//                                         const std::vector<double>& ty,
//                                         std::vector<double>& du1_dx,
//                                         std::vector<double>& du1_dy) const
// {
//     const double pi = 3.14159265358979323846;
//     const int N = (int)tx.size();

//     #pragma omp parallel for
//     for (int i = 0; i < N; ++i) {
//         double s1 = 0.0, s2 = 0.0;
//         for (int j = 0; j < M; ++j) {
//             double dx = tx[i] - proxy_x[j];
//             double dy = ty[i] - proxy_y[j];
//             double R  = dx*dx + dy*dy + eps*eps;
//             double R2 = R * R;
//             // d/dx of (-1/(2pi) * dy / R) summed over j:
//             s1 += xi[j] * ( (1.0 / pi) * (dx * dy) / R2 );
//             // d/dy of (-1/(2pi) * dy / R):
//             s2 += xi[j] * ( -(1.0 / (2.0 * pi)) * (dx*dx - dy*dy + eps*eps) / R2 );
//         }
//         du1_dx[i] += s1;
//         du1_dy[i] += s2;
//     }
// }

// void Periodizer::add_correction_grad_u2(const std::vector<double>& tx,
//                                         const std::vector<double>& ty,
//                                         std::vector<double>& du2_dx,
//                                         std::vector<double>& du2_dy) const
// {
//     const double pi = 3.14159265358979323846;
//     const int N = (int)tx.size();

//     #pragma omp parallel for
//     for (int i = 0; i < N; ++i) {
//         double s1 = 0.0, s2 = 0.0;
//         for (int j = 0; j < M; ++j) {
//             double dx = tx[i] - proxy_x[j];
//             double dy = ty[i] - proxy_y[j];
//             double R  = dx*dx + dy*dy + eps*eps;
//             double R2 = R * R;
//             // d/dx of (1/(2pi) * dx / R):
//             s1 += xi[j] * ( (1.0 / (2.0 * pi)) * (dy*dy - dx*dx + eps*eps) / R2 );
//             // d/dy of (1/(2pi) * dx / R):
//             s2 += xi[j] * ( -(1.0 / pi) * (dx * dy) / R2 );
//         }
//         du2_dx[i] += s1;
//         du2_dy[i] += s2;
//     }
// }