#include "AMRStructure.hpp"
 
int AMRStructure::run() {
    while (iter_num < num_steps) {
        step();
    }
    return 0;
}

int AMRStructure::step() {
    iter_num += 1;
    std::cout << "step " << iter_num << std::endl;

    if (method > 0) {
        euler();
    }
    else {
        rk4();
    }

    t += dt;

    // if remesh:
    if (iter_num % n_steps_remesh == 0) {
        remesh();
    }


    // if dump : write to file
    if (iter_num % n_steps_diag == 0) {
        write_to_file();
    }

    return 0;
}


int AMRStructure::euler() {
    cout << "enter euler" << endl;
    u1s.assign(xs.size(), 0.0);
    u2s.assign(xs.size(), 0.0);
    evaluate_u_field(u1s, u2s, xs, ys, u_weights, t);

    #ifdef DEBUG
        cout << "u1s/u2s first 5:" << endl;
        for (int i = 0; i < std::min<int>(5, (int)u1s.size()); ++i) {
            cout << i << " u1=" << u1s[i] << " u2=" << u2s[i] << endl;
        }
    #endif

    // B evaluation
    b1s.assign(xs.size(), 0.0);
    b2s.assign(xs.size(), 0.0);
    evaluate_b_field(b1s, b2s, xs, ys, b_weights, t);

    // start interpolation
    u1s_grad_x.assign(xs.size(), 0.0);
    u1s_grad_y.assign(xs.size(), 0.0);
    u2s_grad_x.assign(xs.size(), 0.0);
    u2s_grad_y.assign(xs.size(), 0.0);
    b1s_grad_x.assign(xs.size(), 0.0);
    b1s_grad_y.assign(xs.size(), 0.0);
    b2s_grad_x.assign(xs.size(), 0.0);
    b2s_grad_y.assign(xs.size(), 0.0);

    vorticity_grad_x.assign(xs.size(), 0.0);
    vorticity_grad_y.assign(xs.size(), 0.0);
    j_grad_x.assign(xs.size(), 0.0);
    j_grad_y.assign(xs.size(), 0.0);


    // to average: 
     std::vector<int> grad_count(xs.size(), 0);


    for (int panel_ind = 0; panel_ind < panels.size(); panel_ind++) {
        // interpolate_from_panel_to_points
        Panel* panel = &(panels[panel_ind]);
        // only use leaf panels 
        if (panel->child_inds_start > -1) {
            continue;
        }
        const int* panel_point_inds = panel->point_inds;
        double panel_xs[9], panel_ys[9];
        double panel_w0s[9], panel_j0s[9];
        double panel_u1s[9], panel_u2s[9];
        double panel_b1s[9], panel_b2s[9];

        for (int ii = 0; ii < 9; ++ii) {
            int pind = panel_point_inds[ii];
            panel_xs[ii] = xs[pind];
            panel_ys[ii] = ys[pind];
            panel_w0s[ii] = w0s[pind];
            panel_j0s[ii] = j0s[pind];
            panel_u1s[ii] = u1s[pind];
            panel_u2s[ii] = u2s[pind];
            panel_b1s[ii] = b1s[pind];
            panel_b2s[ii] = b2s[pind];
        }

        double panel_dx[9], panel_dy[9];
        for (int ii = 0; ii < 9; ++ii) {
            panel_dx[ii] = panel_xs[ii] - panel_xs[4];
            panel_dy[ii] = panel_ys[ii] - panel_ys[4];
        }

        // Build A for biquadratic basis:
        // [1, x, xy, y, x^2, x^2 y, x^2 y^2, x y^2, y^2]
        Eigen::Matrix<double,9,9> A;
        for (int ii = 0; ii < 9; ++ii) {
            A(ii,0) = 1; A(ii,1) = panel_dx[ii];
            A(ii,2) = panel_dx[ii] * panel_dy[ii];
            A(ii,3) = panel_dy[ii];
            A(ii,4) = panel_dx[ii] * panel_dx[ii];
            A(ii,5) = panel_dx[ii] * panel_dx[ii] * panel_dy[ii];
            A(ii,6) = panel_dx[ii] * panel_dx[ii] * panel_dy[ii] * panel_dy[ii];
            A(ii,7) = panel_dx[ii] * panel_dy[ii] * panel_dy[ii];
            A(ii,8) = panel_dy[ii] * panel_dy[ii];
        }

        // create RHS vector for both w0 and j0, u1s, u2s, b1s, b2s
        Eigen::Map<Eigen::Matrix<double,9,1>> f_w0(panel_w0s);
        Eigen::Map<Eigen::Matrix<double,9,1>> f_j0(panel_j0s);
        Eigen::Map<Eigen::Matrix<double,9,1>> f_u1s(panel_u1s);
        Eigen::Map<Eigen::Matrix<double,9,1>> f_u2s(panel_u2s);
        Eigen::Map<Eigen::Matrix<double,9,1>> f_b1s(panel_b1s);
        Eigen::Map<Eigen::Matrix<double,9,1>> f_b2s(panel_b2s);


        // solve for both w0 and j0 
        Eigen::PartialPivLU<Eigen::Matrix<double, 9, 9>> lu(A);
        // Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 9, 9>> lu(A);
        Eigen::Matrix<double,9,1> c_w0 = lu.solve(f_w0);
        Eigen::Matrix<double,9,1> c_j0 = lu.solve(f_j0);
        Eigen::Matrix<double,9,1> c_u1s = lu.solve(f_u1s);
        Eigen::Matrix<double,9,1> c_u2s = lu.solve(f_u2s);
        Eigen::Matrix<double,9,1> c_b1s = lu.solve(f_b1s);
        Eigen::Matrix<double,9,1> c_b2s = lu.solve(f_b2s);


        // d/dx of basis:
        // d(1)=0, d(x)=1, d(xy)=y, d(y)=0, d(x^2)=2x, d(x^2 y)=2x y,
        // d(x^2 y^2)=2x y^2, d(x y^2)=y^2, d(y^2)=0
        Eigen::Matrix<double,9,9> Dx;
        for (int ii = 0; ii < 9; ++ii) {
            Dx(ii,0) = 0; Dx(ii,1) = 1;
            Dx(ii,2) = panel_dy[ii];
            Dx(ii,3) = 0;
            Dx(ii,4) = 2 * panel_dx[ii];
            Dx(ii,5) = 2 * panel_dx[ii] * panel_dy[ii];
            Dx(ii,6) = 2 * panel_dx[ii] * panel_dy[ii] * panel_dy[ii];
            Dx(ii,7) = panel_dy[ii] * panel_dy[ii];
            Dx(ii,8) = 0;
        }


        Eigen::Matrix<double, 9,1> dx_w0 = Dx * c_w0;
        Eigen::Matrix<double, 9,1> dx_j0 = Dx * c_j0;
        Eigen::Matrix<double, 9,1> dx_u1s = Dx * c_u1s;
        Eigen::Matrix<double, 9,1> dx_u2s = Dx * c_u2s;
        Eigen::Matrix<double, 9,1> dx_b1s = Dx * c_b1s;
        Eigen::Matrix<double, 9,1> dx_b2s = Dx * c_b2s;

        cout << "Dx matrix:\n";
        for (int ii = 0; ii < 9; ++ii) {
            for (int jj = 0; jj < 9; ++jj) {
                cout << Dx(ii,jj) << " ";
            }
            cout << endl;
        }

        cout << "constant coeffcients j0 vector:\n";
        for (int ii = 0; ii < 9; ++ii) {
            cout << c_j0(ii) << endl;
        }



        // d/dy of basis:
        // d(1)=0, d(x)=0, d(xy)=x, d(y)=1, d(x^2)=0, d(x^2 y)=x^2,
        // d(x^2 y^2)=2 x^2 y, d(x y^2)=2 x y, d(y^2)=2y
        Eigen::Matrix<double,9,9> Dy;
        for (int ii = 0; ii < 9; ++ii) {
            Dy(ii,0) = 0; Dy(ii,1) = 0;
            Dy(ii,2) = panel_dx[ii];
            Dy(ii,3) = 1;
            Dy(ii,4) = 0;
            Dy(ii,5) = panel_dx[ii] * panel_dx[ii];
            Dy(ii,6) = 2 * panel_dx[ii] * panel_dx[ii] * panel_dy[ii];
            Dy(ii,7) = 2 * panel_dx[ii] * panel_dy[ii];
            Dy(ii,8) = 2 * panel_dy[ii];
        }

        Eigen::Matrix<double, 9,1> dy_w0 = Dy * c_w0;
        Eigen::Matrix<double, 9,1> dy_j0 = Dy * c_j0;
        Eigen::Matrix<double, 9,1> dy_u1s = Dy * c_u1s;
        Eigen::Matrix<double, 9,1> dy_u2s = Dy * c_u2s;
        Eigen::Matrix<double, 9,1> dy_b1s = Dy * c_b1s;
        Eigen::Matrix<double, 9,1> dy_b2s = Dy * c_b2s;


        for (int ii = 0; ii < 9; ++ii) {
            int pind = panel_point_inds[ii];
            vorticity_grad_x[pind] += dx_w0(ii,0);
            j_grad_x[pind]         += dx_j0(ii,0);
            u1s_grad_x[pind]       += dx_u1s(ii,0);
            u2s_grad_x[pind]       += dx_u2s(ii,0);
            b1s_grad_x[pind]       += dx_b1s(ii,0);
            b2s_grad_x[pind]       += dx_b2s(ii,0);

            vorticity_grad_y[pind] += dy_w0(ii,0);
            j_grad_y[pind]         += dy_j0(ii,0);
            u1s_grad_y[pind]       += dy_u1s(ii,0);
            u2s_grad_y[pind]       += dy_u2s(ii,0); 
            b1s_grad_y[pind]       += dy_b1s(ii,0);
            b2s_grad_y[pind]       += dy_b2s(ii,0);

            grad_count[pind] += 1;
        }


        // print out the gradient for each panel 
        cout << "leaf panel: " << panel->panel_ind <<endl;
        for (int ii = 0; ii < 9; ++ii) {
            int pind = panel_point_inds[ii];
            cout << "i=" << ii
            << " x=" << xs[pind]
            << " y=" << ys[pind]
            << " b1=" << b1s[pind]
            << " b2=" << b2s[pind]
            << " j=" << j0s[pind]
            << " w=" << w0s[pind]
            << " w_dx=" << dx_w0[ii]
            << " w_dy=" << dy_w0[ii]
            << " j_dx=" << dx_j0[ii]
            << " j_dy=" << dy_j0[ii]
            << "\n";
        }

    }

    // average gradients
    for (int i = 0; i < xs.size(); ++i) {
        if (grad_count[i] > 0) {
            const double inv = 1.0 / static_cast<double>(grad_count[i]);
            vorticity_grad_x[i] *= inv; vorticity_grad_y[i] *= inv;
            j_grad_x[i] *= inv;         j_grad_y[i] *= inv;
            u1s_grad_x[i] *= inv;       u1s_grad_y[i] *= inv;
            u2s_grad_x[i] *= inv;       u2s_grad_y[i] *= inv;
            b1s_grad_x[i] *= inv;       b1s_grad_y[i] *= inv;
            b2s_grad_x[i] *= inv;       b2s_grad_y[i] *= inv;
        }
    }

    // treat boundary
    for (int panel_ind = 0; panel_ind < panels.size(); panel_ind++) {
        Panel* panel = &(panels[panel_ind]);
        // only use leaf panels 
        if (panel->child_inds_start > -1) {
            continue;
        }
        if (panel->is_left_bdry) {
            cout << " left bdry panel: " << panel->panel_ind << endl;
            const int* panel_point_inds = panel->point_inds;
            Panel* left_panel = &(panels[panel->left_nbr_ind]);
            const int* left_panel_point_inds = left_panel->point_inds;
            // three left boundary points of the panel, point index in the panel is 0,1,2
            // three right boundary points of the left panel, point index in the panel is 6,7,8
            int pind_0 = panel_point_inds[0];
            int pind_1 = panel_point_inds[1];
            int pind_2 = panel_point_inds[2];

            int left_pind_6 = left_panel_point_inds[6];
            int left_pind_7 = left_panel_point_inds[7];
            int left_pind_8 = left_panel_point_inds[8];

            vorticity_grad_x[pind_0] = (vorticity_grad_x[pind_0] + vorticity_grad_x[left_pind_6])/2;
            vorticity_grad_y[pind_0] = (vorticity_grad_y[pind_0] + vorticity_grad_y[left_pind_6])/2;
            j_grad_x[pind_0] = (j_grad_x[pind_0] + j_grad_x[left_pind_6])/2;
            j_grad_y[pind_0] = (j_grad_y[pind_0] + j_grad_y[left_pind_6])/2;
            u1s_grad_x[pind_0] = (u1s_grad_x[pind_0] + u1s_grad_x[left_pind_6])/2;
            u1s_grad_y[pind_0] = (u1s_grad_y[pind_0] + u1s_grad_y[left_pind_6])/2;
            u2s_grad_x[pind_0] = (u2s_grad_x[pind_0] + u2s_grad_x[left_pind_6])/2;
            u2s_grad_y[pind_0] = (u2s_grad_y[pind_0] + u2s_grad_y[left_pind_6])/2;
            b1s_grad_x[pind_0] = (b1s_grad_x[pind_0] + b1s_grad_x[left_pind_6])/2;
            b1s_grad_y[pind_0] = (b1s_grad_y[pind_0] + b1s_grad_y[left_pind_6])/2;
            b2s_grad_x[pind_0] = (b2s_grad_x[pind_0] + b2s_grad_x[left_pind_6])/2;
            b2s_grad_y[pind_0] = (b2s_grad_y[pind_0] + b2s_grad_y[left_pind_6])/2;

            vorticity_grad_x[pind_1] = (vorticity_grad_x[pind_1] + vorticity_grad_x[left_pind_7])/2;
            vorticity_grad_y[pind_1] = (vorticity_grad_y[pind_1] + vorticity_grad_y[left_pind_7])/2;
            j_grad_x[pind_1] = (j_grad_x[pind_1] + j_grad_x[left_pind_7])/2;
            j_grad_y[pind_1] = (j_grad_y[pind_1] + j_grad_y[left_pind_7])/2;
            u1s_grad_x[pind_1] = (u1s_grad_x[pind_1] + u1s_grad_x[left_pind_7])/2;
            u1s_grad_y[pind_1] = (u1s_grad_y[pind_1] + u1s_grad_y[left_pind_7])/2;
            u2s_grad_x[pind_1] = (u2s_grad_x[pind_1] + u2s_grad_x[left_pind_7])/2;
            u2s_grad_y[pind_1] = (u2s_grad_y[pind_1] + u2s_grad_y[left_pind_7])/2;
            b1s_grad_x[pind_1] = (b1s_grad_x[pind_1] + b1s_grad_x[left_pind_7])/2;
            b1s_grad_y[pind_1] = (b1s_grad_y[pind_1] + b1s_grad_y[left_pind_7])/2;
            b2s_grad_x[pind_1] = (b2s_grad_x[pind_1] + b2s_grad_x[left_pind_7])/2;
            b2s_grad_y[pind_1] = (b2s_grad_y[pind_1] + b2s_grad_y[left_pind_7])/2;

            vorticity_grad_x[pind_2] = (vorticity_grad_x[pind_2] + vorticity_grad_x[left_pind_8])/2;
            vorticity_grad_y[pind_2] = (vorticity_grad_y[pind_2] + vorticity_grad_y[left_pind_8])/2;
            j_grad_x[pind_2] = (j_grad_x[pind_2] + j_grad_x[left_pind_8])/2;
            j_grad_y[pind_2] = (j_grad_y[pind_2] + j_grad_y[left_pind_8])/2;
            u1s_grad_x[pind_2] = (u1s_grad_x[pind_2] + u1s_grad_x[left_pind_8])/2;
            u1s_grad_y[pind_2] = (u1s_grad_y[pind_2] + u1s_grad_y[left_pind_8])/2;
            u2s_grad_x[pind_2] = (u2s_grad_x[pind_2] + u2s_grad_x[left_pind_8])/2;
            u2s_grad_y[pind_2] = (u2s_grad_y[pind_2] + u2s_grad_y[left_pind_8])/2;
            b1s_grad_x[pind_2] = (b1s_grad_x[pind_2] + b1s_grad_x[left_pind_8])/2;
            b1s_grad_y[pind_2] = (b1s_grad_y[pind_2] + b1s_grad_y[left_pind_8])/2;
            b2s_grad_x[pind_2] = (b2s_grad_x[pind_2] + b2s_grad_x[left_pind_8])/2;
            b2s_grad_y[pind_2] = (b2s_grad_y[pind_2] + b2s_grad_y[left_pind_8])/2;


            // right boundary have the same value
            vorticity_grad_x[left_pind_6] = vorticity_grad_x[pind_0];
            vorticity_grad_y[left_pind_6] = vorticity_grad_y[pind_0];
            j_grad_x[left_pind_6] = j_grad_x[pind_0];
            j_grad_y[left_pind_6] = j_grad_y[pind_0];
            u1s_grad_x[left_pind_6] = u1s_grad_x[pind_0];
            u1s_grad_y[left_pind_6] = u1s_grad_y[pind_0];
            u2s_grad_x[left_pind_6] = u2s_grad_x[pind_0];
            u2s_grad_y[left_pind_6] = u2s_grad_y[pind_0];
            b1s_grad_x[left_pind_6] = b1s_grad_x[pind_0];
            b1s_grad_y[left_pind_6] = b1s_grad_y[pind_0];
            b2s_grad_x[left_pind_6] = b2s_grad_x[pind_0];
            b2s_grad_y[left_pind_6] = b2s_grad_y[pind_0];

            vorticity_grad_x[left_pind_7] = vorticity_grad_x[pind_1];
            vorticity_grad_y[left_pind_7] = vorticity_grad_y[pind_1];
            j_grad_x[left_pind_7] = j_grad_x[pind_1];
            j_grad_y[left_pind_7] = j_grad_y[pind_1];
            u1s_grad_x[left_pind_7] = u1s_grad_x[pind_1];
            u1s_grad_y[left_pind_7] = u1s_grad_y[pind_1];
            u2s_grad_x[left_pind_7] = u2s_grad_x[pind_1];
            u2s_grad_y[left_pind_7] = u2s_grad_y[pind_1];
            b1s_grad_x[left_pind_7] = b1s_grad_x[pind_1];
            b1s_grad_y[left_pind_7] = b1s_grad_y[pind_1];
            b2s_grad_x[left_pind_7] = b2s_grad_x[pind_1];
            b2s_grad_y[left_pind_7] = b2s_grad_y[pind_1];

            vorticity_grad_x[left_pind_8] = vorticity_grad_x[pind_2];
            vorticity_grad_y[left_pind_8] = vorticity_grad_y[pind_2];
            j_grad_x[left_pind_8] = j_grad_x[pind_2];
            j_grad_y[left_pind_8] = j_grad_y[pind_2];
            u1s_grad_x[left_pind_8] = u1s_grad_x[pind_2];
            u1s_grad_y[left_pind_8] = u1s_grad_y[pind_2];
            u2s_grad_x[left_pind_8] = u2s_grad_x[pind_2];
            u2s_grad_y[left_pind_8] = u2s_grad_y[pind_2];
            b1s_grad_x[left_pind_8] = b1s_grad_x[pind_2];
            b1s_grad_y[left_pind_8] = b1s_grad_y[pind_2];
            b2s_grad_x[left_pind_8] = b2s_grad_x[pind_2];
            b2s_grad_y[left_pind_8] = b2s_grad_y[pind_2];
        }
    }


    // printing:
    // for (int i = 0; i < static_cast<int>(xs.size()); ++i) {
    //     std::cout
    //         << "i=" << i
    //         << " x=" << xs[i]
    //         << " y=" << ys[i]
    //         << " u1=" << u1s[i]
    //         << " u2=" << u2s[i]
    //         << " b1=" << b1s[i]
    //         << " b2=" << b2s[i]
    //         << " j=" << j0s[i]
    //         << " w=" << w0s[i]
    //         << " w_dx=" << vorticity_grad_x[i]
    //         << " w_dy=" << vorticity_grad_y[i]
    //         << " j_dx=" << j_grad_x[i]
    //         << " j_dy=" << j_grad_y[i]
    //         << "\n";
    // }

    // calculate source terms
    B_dot_grad_j.assign(xs.size(), 0.0);
    for (int i = 0; i < xs.size(); ++i) {
        B_dot_grad_j[i] = b1s[i] * j_grad_x[i] + b2s[i] * j_grad_y[i];
    }

    B_dot_grad_vorticity.assign(xs.size(), 0.0);
    for (int i = 0; i < xs.size(); ++i) {
        B_dot_grad_vorticity[i] = b1s[i] * vorticity_grad_x[i] + b2s[i] * vorticity_grad_y[i];
    }

    B_grad_x_dot_u2_grad.assign(xs.size(), 0.0);
    for (int i = 0; i < xs.size(); ++i) {
        B_grad_x_dot_u2_grad[i] = b1s_grad_x[i] * u2s_grad_x[i] + b2s_grad_x[i] * u2s_grad_y[i];
    }

    B_grad_y_dot_u1_grad.assign(xs.size(), 0.0);
    for (int i = 0; i < xs.size(); ++i) {
        B_grad_y_dot_u1_grad[i] = b1s_grad_y[i] * u1s_grad_x[i] + b2s_grad_y[i] * u1s_grad_y[i];
    }
    // start pushing, particle position, vorticity and current density 
    for (int i = 0; i < xs.size(); i++) {
        xs[i] += dt * u1s[i];
        ys[i] += dt * u2s[i];
        w0s[i] += dt * (nu * vorticity_laplacian[i] + B_dot_grad_j[i]);
        j0s[i] += dt * (mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i]);
    }

    return 0;
}




int AMRStructure::rk4() {
    std::cout << "enter rk4" << std::endl;

    const int xs_size = (int)xs.size();

    // Stage storage: k1, k2, k3, k4 for x, y, w and j
    std::vector<double> k1_dws(xs_size, 0.0), k1_djs(xs_size, 0.0); 
    std::vector<double> k2_u1s(xs_size, 0.0), k2_u2s(xs_size, 0.0), k2_dws(xs_size, 0.0), k2_djs(xs_size, 0.0);
    std::vector<double> k3_u1s(xs_size, 0.0), k3_u2s(xs_size, 0.0), k3_dws(xs_size, 0.0), k3_djs(xs_size, 0.0);
    std::vector<double> k4_u1s(xs_size, 0.0), k4_u2s(xs_size, 0.0), k4_dws(xs_size, 0.0), k4_djs(xs_size, 0.0);

    // Temporary positions
    std::vector<double> xs2(xs_size), ys2(xs_size), ws2(xs_size), js2(xs_size);
    std::vector<double> xs3(xs_size), ys3(xs_size), ws3(xs_size), js3(xs_size);
    std::vector<double> xs4(xs_size), ys4(xs_size), ws4(xs_size), js4(xs_size);

    // -------------------------
    // Stage 1 : k1 = (u1s, u2s, k1_dws, k1_djs)
    // -------------------------
    // U evaluation
    u1s.assign(u1s.size(), 0.0);
    u2s.assign(u2s.size(), 0.0);
    evaluate_u_field(u1s, u2s, xs, ys, u_weights, t);

    // B evaluation
    b1s.assign(b1s.size(), 0.0);
    b2s.assign(b2s.size(), 0.0);
    evaluate_b_field(b1s, b2s, xs, ys, b_weights, t);

    // // u1s_grad, u2s_grad, b1s_grad, b2s_grad evaluation
    // u1s_grad_x.assign(xs.size(), 0.0);
    // u1s_grad_y.assign(xs.size(), 0.0);
    // evaluate_u1s_grad(u1s_grad_x, u1s_grad_y, xs, ys, u_weights, t);


    // u2s_grad_x.assign(xs.size(), 0.0);
    // u2s_grad_y.assign(xs.size(), 0.0);

    // evaluate_u2s_grad(u2s_grad_x, u2s_grad_y, xs, ys, u_weights, t);

    // b1s_grad_x.assign(xs.size(), 0.0);
    // b1s_grad_y.assign(xs.size(), 0.0);

    // evaluate_b1s_grad(b1s_grad_x, b1s_grad_y, xs, ys, b_weights, t);

    // b2s_grad_x.assign(xs.size(), 0.0);
    // b2s_grad_y.assign(xs.size(), 0.0);

    // evaluate_b2s_grad(b2s_grad_x, b2s_grad_y, xs, ys, b_weights, t);

    // // vortex_grad, j_grad evaluation
    // vorticity_grad_x.assign(xs.size(), 0.0);
    // vorticity_grad_y.assign(xs.size(), 0.0);

    // evaluate_vorticity_grad(vorticity_grad_x, vorticity_grad_y, xs, ys, u_weights, t);

    // j_grad_x.assign(xs.size(), 0.0);
    // j_grad_y.assign(xs.size(), 0.0);

    // evaluate_j_grad(j_grad_x, j_grad_y, xs, ys, b_weights, t);

    // // vortex_laplacian, j_laplacian evaluation

    // vorticity_laplacian.assign(xs.size(), 0.0);
    // j_laplacian.assign(xs.size(), 0.0);

    // std::vector<double> vorticity_none_local(xs.size(), 0.0);
    // std::vector<double> j_none_local(xs.size(), 0.0);

    // evaluate_vorticity_laplacian(vorticity_laplacian, vorticity_none_local, xs, ys, u_weights, t);
    // evaluate_j_laplacian(j_laplacian, j_none_local, xs, ys, b_weights, t);

    // // calculate source terms
    // B_dot_grad_j.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_dot_grad_j[i] = b1s[i] * j_grad_x[i] + b2s[i] * j_grad_y[i];
    // }

    // B_dot_grad_vorticity.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_dot_grad_vorticity[i] = b1s[i] * vorticity_grad_x[i] + b2s[i] * vorticity_grad_y[i];
    // }

    // B_grad_x_dot_u2_grad.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_grad_x_dot_u2_grad[i] = b1s_grad_x[i] * u2s_grad_x[i] + b2s_grad_x[i] * u2s_grad_y[i];
    // }

    // B_grad_y_dot_u1_grad.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_grad_y_dot_u1_grad[i] = b1s_grad_y[i] * u1s_grad_x[i] + b2s_grad_y[i] * u1s_grad_y[i];
    // }

    // // rebuild u_weights amd b_weights : first divide by w/j before pushing
    // // then multiply new w/j after pushing 
    // // or just use weights[ii] * ws2, js2
    
    // // // divide w/j before pushing 
    // // for (int i = 0; i < u_weights.size(); i++) {
    // //     u_weights[i] = u_weights[i] / w0s[i];
    // //     b_weights[i] = b_weights[i] / j0s[i];
    // // }

    // // pushing 
    // for (int i = 0; i < xs.size(); i++) {
    //     xs2[i] = xs[i] + 0.5 * dt * u1s[i];
    // }
    // for (int i = 0; i < ys.size(); i++) {
    //     ys2[i] = ys[i] + 0.5 * dt * u2s[i];
    // }
    // for (int i = 0; i < xs.size(); i++) {
    //     k1_dws[i] = nu * vorticity_laplacian[i] + B_dot_grad_j[i]; 
    //     ws2[i] = w0s[i] + 0.5 * dt * k1_dws[i];
    // }
    // for (int i = 0; i < xs.size(); i++) {
    //     k1_djs[i] = mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i];
    //     js2[i] = j0s[i] + 0.5 * dt * k1_djs[i];
    // }

    // // multiply new w/j after pushing 
    // for (int i = 0; i < u_weights.size(); i++) {
    //     u_weights[i] = weights[i] * ws2[i];
    //     b_weights[i] = weights[i] * js2[i];
    // }


    // // -------------------------
    // // Stage 2 : k2_u1s, k2_u2s, k2_dws, k2_djs
    // // -------------------------

    // evaluate_u_field(k2_u1s, k2_u2s, xs2, ys2, u_weights, t + 0.5 * dt);

    // // B evaluation
    // b1s.assign(b1s.size(), 0.0);
    // b2s.assign(b2s.size(), 0.0);
    // evaluate_b_field(b1s, b2s, xs2, ys2, b_weights, t+ 0.5 * dt);

    // // u1s_grad, u2s_grad, b1s_grad, b2s_grad evaluation
    // u1s_grad_x.assign(xs.size(), 0.0);
    // u1s_grad_y.assign(xs.size(), 0.0);
    // evaluate_u1s_grad(u1s_grad_x, u1s_grad_y, xs2, ys2, u_weights, t+ 0.5 * dt);


    // u2s_grad_x.assign(xs.size(), 0.0);
    // u2s_grad_y.assign(xs.size(), 0.0);

    // evaluate_u2s_grad(u2s_grad_x, u2s_grad_y, xs2, ys2, u_weights, t+ 0.5 * dt);

    // b1s_grad_x.assign(xs.size(), 0.0);
    // b1s_grad_y.assign(xs.size(), 0.0);

    // evaluate_b1s_grad(b1s_grad_x, b1s_grad_y, xs2, ys2, b_weights, t+ 0.5 * dt);

    // b2s_grad_x.assign(xs.size(), 0.0);
    // b2s_grad_y.assign(xs.size(), 0.0);

    // evaluate_b2s_grad(b2s_grad_x, b2s_grad_y, xs2, ys2, b_weights, t+ 0.5 * dt);

    // // vortex_grad, j_grad evaluation
    // vorticity_grad_x.assign(xs.size(), 0.0);
    // vorticity_grad_y.assign(xs.size(), 0.0);

    // evaluate_vorticity_grad(vorticity_grad_x, vorticity_grad_y, xs2, ys2, u_weights, t+ 0.5 * dt);

    // j_grad_x.assign(xs.size(), 0.0);
    // j_grad_y.assign(xs.size(), 0.0);

    // evaluate_j_grad(j_grad_x, j_grad_y, xs2, ys2, b_weights, t+ 0.5 * dt);

    // // vortex_laplacian, j_laplacian evaluation

    // vorticity_laplacian.assign(xs.size(), 0.0);
    // j_laplacian.assign(xs.size(), 0.0);

    // vorticity_none_local.assign(xs.size(), 0.0);
    // j_none_local.assign(xs.size(), 0.0);

    // evaluate_vorticity_laplacian(vorticity_laplacian, vorticity_none_local, xs2, ys2, u_weights, t+ 0.5 * dt);
    // evaluate_j_laplacian(j_laplacian, j_none_local, xs2, ys2, b_weights, t+ 0.5 * dt);

    // // calculate source terms
    // B_dot_grad_j.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_dot_grad_j[i] = b1s[i] * j_grad_x[i] + b2s[i] * j_grad_y[i];
    // }

    // B_dot_grad_vorticity.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_dot_grad_vorticity[i] = b1s[i] * vorticity_grad_x[i] + b2s[i] * vorticity_grad_y[i];
    // }

    // B_grad_x_dot_u2_grad.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_grad_x_dot_u2_grad[i] = b1s_grad_x[i] * u2s_grad_x[i] + b2s_grad_x[i] * u2s_grad_y[i];
    // }

    // B_grad_y_dot_u1_grad.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_grad_y_dot_u1_grad[i] = b1s_grad_y[i] * u1s_grad_x[i] + b2s_grad_y[i] * u1s_grad_y[i];
    // }


    // // pushing 
    // for (int i = 0; i < xs.size(); i++) {
    //     xs3[i] = xs[i] + 0.5 * dt * k2_u1s[i];
    // }
    // for (int i = 0; i < ys.size(); i++) {
    //     ys3[i] = ys[i] + 0.5 * dt * k2_u2s[i];
    // }
    // for (int i = 0; i < xs.size(); i++) {
    //     k2_dws[i] = nu * vorticity_laplacian[i] + B_dot_grad_j[i]; 
    //     ws3[i] = w0s[i] + 0.5 * dt * k2_dws[i];
    // }
    // for (int i = 0; i < xs.size(); i++) {
    //     k2_djs[i] = mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i];
    //     js3[i] = j0s[i] + 0.5 * dt * k2_djs[i];
    // }

    // // multiply new w/j after pushing 
    // for (int i = 0; i < u_weights.size(); i++) {
    //     u_weights[i] = weights[i] * ws3[i];
    //     b_weights[i] = weights[i] * js3[i];
    // }




    // // -------------------------
    // // Stage 3: k3_u1s, k3_u2s, k3_dws, k3_djs
    // // -------------------------

    // evaluate_u_field(k3_u1s, k3_u2s, xs3, ys3, u_weights, t + 0.5 * dt);

    // // B evaluation
    // b1s.assign(b1s.size(), 0.0);
    // b2s.assign(b2s.size(), 0.0);
    // evaluate_b_field(b1s, b2s, xs3, ys3, b_weights, t+ 0.5 * dt);

    // // u1s_grad, u2s_grad, b1s_grad, b2s_grad evaluation
    // u1s_grad_x.assign(xs.size(), 0.0);
    // u1s_grad_y.assign(xs.size(), 0.0);
    // evaluate_u1s_grad(u1s_grad_x, u1s_grad_y, xs3, ys3, u_weights, t+ 0.5 * dt);


    // u2s_grad_x.assign(xs.size(), 0.0);
    // u2s_grad_y.assign(xs.size(), 0.0);

    // evaluate_u2s_grad(u2s_grad_x, u2s_grad_y, xs3, ys3, u_weights, t+ 0.5 * dt);

    // b1s_grad_x.assign(xs.size(), 0.0);
    // b1s_grad_y.assign(xs.size(), 0.0);

    // evaluate_b1s_grad(b1s_grad_x, b1s_grad_y, xs3, ys3, b_weights, t+ 0.5 * dt);

    // b2s_grad_x.assign(xs.size(), 0.0);
    // b2s_grad_y.assign(xs.size(), 0.0);

    // evaluate_b2s_grad(b2s_grad_x, b2s_grad_y, xs3, ys3, b_weights, t+ 0.5 * dt);

    // // vortex_grad, j_grad evaluation
    // vorticity_grad_x.assign(xs.size(), 0.0);
    // vorticity_grad_y.assign(xs.size(), 0.0);

    // evaluate_vorticity_grad(vorticity_grad_x, vorticity_grad_y, xs3, ys3, u_weights, t+ 0.5 * dt);

    // j_grad_x.assign(xs.size(), 0.0);
    // j_grad_y.assign(xs.size(), 0.0);

    // evaluate_j_grad(j_grad_x, j_grad_y, xs3, ys3, b_weights, t+ 0.5 * dt);

    // // vortex_laplacian, j_laplacian evaluation

    // vorticity_laplacian.assign(xs.size(), 0.0);
    // j_laplacian.assign(xs.size(), 0.0);

    // vorticity_none_local.assign(xs.size(), 0.0);
    // j_none_local.assign(xs.size(), 0.0);

    // evaluate_vorticity_laplacian(vorticity_laplacian, vorticity_none_local, xs3, ys3, u_weights, t+ 0.5 * dt);
    // evaluate_j_laplacian(j_laplacian, j_none_local, xs3, ys3, b_weights, t+ 0.5 * dt);

    // // calculate source terms
    // B_dot_grad_j.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_dot_grad_j[i] = b1s[i] * j_grad_x[i] + b2s[i] * j_grad_y[i];
    // }

    // B_dot_grad_vorticity.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_dot_grad_vorticity[i] = b1s[i] * vorticity_grad_x[i] + b2s[i] * vorticity_grad_y[i];
    // }

    // B_grad_x_dot_u2_grad.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_grad_x_dot_u2_grad[i] = b1s_grad_x[i] * u2s_grad_x[i] + b2s_grad_x[i] * u2s_grad_y[i];
    // }

    // B_grad_y_dot_u1_grad.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_grad_y_dot_u1_grad[i] = b1s_grad_y[i] * u1s_grad_x[i] + b2s_grad_y[i] * u1s_grad_y[i];
    // }

    // // pushing 
    // for (int i = 0; i < xs.size(); i++) {
    //     xs4[i] = xs[i] + dt * k3_u1s[i];
    // }
    // for (int i = 0; i < ys.size(); i++) {
    //     ys4[i] = ys[i] + dt * k3_u2s[i];
    // }
    // for (int i = 0; i < xs.size(); i++) {
    //     k3_dws[i] = nu * vorticity_laplacian[i] + B_dot_grad_j[i]; 
    //     ws4[i] = w0s[i] + dt * k3_dws[i];
    // }
    // for (int i = 0; i < xs.size(); i++) {
    //     k3_djs[i] = mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i];
    //     js4[i] = j0s[i] + dt * k3_djs[i];
    // }

    // // multiply new w/j after pushing 
    // for (int i = 0; i < u_weights.size(); i++) {
    //     u_weights[i] = weights[i] * ws4[i];
    //     b_weights[i] = weights[i] * js4[i];
    // }




    // // -------------------------
    // // Stage 4: k4_u1s, k4_u2s, k4_dws, k4_djs
    // // -------------------------
    // // evaluate_u_field(k4x, k4y, x4, y4, u_weights, t + dt);

    // evaluate_u_field(k4_u1s, k4_u2s, xs4, ys4, u_weights, t + dt);

    // // B evaluation
    // b1s.assign(b1s.size(), 0.0);
    // b2s.assign(b2s.size(), 0.0);
    // evaluate_b_field(b1s, b2s, xs4, ys4, b_weights, t + dt);

    // // u1s_grad, u2s_grad, b1s_grad, b2s_grad evaluation
    // u1s_grad_x.assign(xs.size(), 0.0);
    // u1s_grad_y.assign(xs.size(), 0.0);
    // evaluate_u1s_grad(u1s_grad_x, u1s_grad_y, xs4, ys4, u_weights, t+ dt);


    // u2s_grad_x.assign(xs.size(), 0.0);
    // u2s_grad_y.assign(xs.size(), 0.0);

    // evaluate_u2s_grad(u2s_grad_x, u2s_grad_y, xs4, ys4, u_weights, t+ dt);

    // b1s_grad_x.assign(xs.size(), 0.0);
    // b1s_grad_y.assign(xs.size(), 0.0);

    // evaluate_b1s_grad(b1s_grad_x, b1s_grad_y, xs4, ys4, b_weights, t+ dt);

    // b2s_grad_x.assign(xs.size(), 0.0);
    // b2s_grad_y.assign(xs.size(), 0.0);

    // evaluate_b2s_grad(b2s_grad_x, b2s_grad_y, xs4, ys4, b_weights, t+ dt);

    // // vortex_grad, j_grad evaluation
    // vorticity_grad_x.assign(xs.size(), 0.0);
    // vorticity_grad_y.assign(xs.size(), 0.0);

    // evaluate_vorticity_grad(vorticity_grad_x, vorticity_grad_y, xs4, ys4, u_weights, t+ dt);

    // j_grad_x.assign(xs.size(), 0.0);
    // j_grad_y.assign(xs.size(), 0.0);

    // evaluate_j_grad(j_grad_x, j_grad_y, xs4, ys4, b_weights, t+ dt);

    // // vortex_laplacian, j_laplacian evaluation

    // vorticity_laplacian.assign(xs.size(), 0.0);
    // j_laplacian.assign(xs.size(), 0.0);

    // vorticity_none_local.assign(xs.size(), 0.0);
    // j_none_local.assign(xs.size(), 0.0);

    // evaluate_vorticity_laplacian(vorticity_laplacian, vorticity_none_local, xs4, ys4, u_weights, t+ dt);
    // evaluate_j_laplacian(j_laplacian, j_none_local, xs4, ys4, b_weights, t+ dt);

    // // calculate source terms
    // B_dot_grad_j.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_dot_grad_j[i] = b1s[i] * j_grad_x[i] + b2s[i] * j_grad_y[i];
    // }

    // B_dot_grad_vorticity.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_dot_grad_vorticity[i] = b1s[i] * vorticity_grad_x[i] + b2s[i] * vorticity_grad_y[i];
    // }

    // B_grad_x_dot_u2_grad.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_grad_x_dot_u2_grad[i] = b1s_grad_x[i] * u2s_grad_x[i] + b2s_grad_x[i] * u2s_grad_y[i];
    // }

    // B_grad_y_dot_u1_grad.assign(xs.size(), 0.0);
    // for (int i = 0; i < xs.size(); ++i) {
    //     B_grad_y_dot_u1_grad[i] = b1s_grad_y[i] * u1s_grad_x[i] + b2s_grad_y[i] * u1s_grad_y[i];
    // }

    // // source terms: 
    // for (int i = 0; i < xs.size(); i++) {
    //     k4_dws[i] = nu * vorticity_laplacian[i] + B_dot_grad_j[i]; 
    // }
    // for (int i = 0; i < xs.size(); i++) {
    //     k4_djs[i] = mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i];
    // }


    // // -------------------------
    // // Combine to update positions
    // // x_{n+1} = x_n + dt/6*(k1 + 2k2 + 2k3 + k4)
    // // (u1s, u2s, k1_dws, k1_djs), (k2_u1s, k2_u2s, k2_dws, k2_djs), (k3_u1s, k3_u2s, k3_dws, k3_djs), (k4_u1s, k4_u2s, k4_dws, k4_djs)
    // // -------------------------
    // for (int i = 0; i < xs_size; ++i) {
    //     xs[i] += (dt / 6.0) * (u1s[i] + 2.0 * k2_u1s[i] + 2.0 * k3_u1s[i] + k4_u1s[i]);
    //     ys[i] += (dt / 6.0) * (u2s[i] + 2.0 * k2_u2s[i] + 2.0 * k3_u2s[i] + k4_u2s[i]);
    //     w0s[i] += (dt / 6.0) * (k1_dws[i] + 2.0 * k2_dws[i] + 2.0 * k3_dws[i] + k4_dws[i]);
    //     j0s[i] += (dt / 6.0) * (k1_djs[i] + 2.0 * k2_djs[i] + 2.0 * k3_djs[i] + k4_djs[i]);
    // }


    return 0;
}
