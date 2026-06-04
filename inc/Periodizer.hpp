#ifndef PERIODIZER_HPP
#define PERIODIZER_HPP

#include <vector>
#include <Eigen/Dense>
#include "FieldStructure.hpp"

class Periodizer {
public:
    Periodizer(double x_min, double x_max,
               double y_min, double y_max,
               double eps,
               Field* free_space_kernel,
               int M = 80,
               int m = 22,
               double Rp_factor = 1.4);

    // Call once at startup: fills Q and factorises it.
    void precompute_Q();

    // Given the 3x3 near-sum velocity already evaluated at the 4 wall pairs
    // (8m values, see comment in .cpp), solve Q xi = g.
    void compute_xi(const std::vector<double>& near_u_on_walls_8m);

    // Add xi * K_proxy to velocity at arbitrary target points.
    void add_correction(const std::vector<double>& tx,
                        const std::vector<double>& ty,
                        std::vector<double>& u1,
                        std::vector<double>& u2) const;

    // // Add the gradient of the proxy correction at target points.
    // // Outputs the (x,y) components of grad(u1) and grad(u2).
    // void add_correction_grad_u1(const std::vector<double>& tx,
    //                             const std::vector<double>& ty,
    //                             std::vector<double>& du1_dx,
    //                             std::vector<double>& du1_dy) const;

    // void add_correction_grad_u2(const std::vector<double>& tx,
    //                             const std::vector<double>& ty,
    //                             std::vector<double>& du2_dx,
    //                             std::vector<double>& du2_dy) const;                  

    // Accessors so the caller can build the near-sum wall targets
    const std::vector<double>& get_wall_xs() const { return wall_xs_4m; }
    const std::vector<double>& get_wall_ys() const { return wall_ys_4m; }
    int get_m() const { return m; }
    int get_M() const { return M; }

private:
    double Lx, Ly, cx, cy, Rp;
    double eps;
    int M, m;
    Field* kernel;   // not owned — just used for type; actual evaluation is inline

    std::vector<double> proxy_x, proxy_y;   // size M
    std::vector<double> wall_xs_4m, wall_ys_4m;  // size 4m

    Eigen::MatrixXd Q_matrix;                                   // 4m x M
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Q_solver;       // factorization

    std::vector<double> xi;   // size M, updated each compute_xi call
};

#endif // PERIODIZER_HPP