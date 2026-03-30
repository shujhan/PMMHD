#include "AMRStructure.hpp"
#include <iomanip>
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




// int AMRStructure::euler() {
//     cout << "enter euler" << endl;
//     u1s.assign(u1s.size(), 0.0);
//     u2s.assign(u2s.size(), 0.0);
//     evaluate_u_field(u1s, u2s, xs, ys, u_weights, t);

//     #ifdef DEBUG
//         cout << "u1s/u2s first 5:" << endl;
//         for (int i = 0; i < std::min<int>(5, (int)u1s.size()); ++i) {
//             cout << i << " u1=" << u1s[i] << " u2=" << u2s[i] << endl;
//         }
//     #endif

//     // B evaluation
//     b1s.assign(b1s.size(), 0.0);
//     b2s.assign(b2s.size(), 0.0);
//     evaluate_b_field(b1s, b2s, xs, ys, b_weights, t);

//     // u1s_grad, u2s_grad, b1s_grad, b2s_grad evaluation

//     u1s_grad_x.assign(xs.size(), 0.0);
//     u1s_grad_y.assign(xs.size(), 0.0);

//     evaluate_u1s_grad(u1s_grad_x, u1s_grad_y, xs, ys, u_weights, t);

//     u2s_grad_x.assign(xs.size(), 0.0);
//     u2s_grad_y.assign(xs.size(), 0.0);

//     evaluate_u2s_grad(u2s_grad_x, u2s_grad_y, xs, ys, u_weights, t);

//     b1s_grad_x.assign(xs.size(), 0.0);
//     b1s_grad_y.assign(xs.size(), 0.0);

//     evaluate_b1s_grad(b1s_grad_x, b1s_grad_y, xs, ys, b_weights, t);

//     b2s_grad_x.assign(xs.size(), 0.0);
//     b2s_grad_y.assign(xs.size(), 0.0);

//     evaluate_b2s_grad(b2s_grad_x, b2s_grad_y, xs, ys, b_weights, t);

//     // vortex_grad, j_grad evaluation
//     vorticity_grad_x.assign(xs.size(), 0.0);
//     vorticity_grad_y.assign(xs.size(), 0.0);

//     evaluate_vorticity_grad(vorticity_grad_x, vorticity_grad_y, xs, ys, u_weights, t);

//     j_grad_x.assign(xs.size(), 0.0);
//     j_grad_y.assign(xs.size(), 0.0);

//     evaluate_j_grad(j_grad_x, j_grad_y, xs, ys, b_weights, t);
//     for (int i = 0; i < xs.size(); ++i) {
//         cout << i << " j_grad_x=" << j_grad_x[i] << " j_grad_y=" << j_grad_y[i] << endl;
//     }

//     // vortex_laplacian, j_laplacian evaluation

//     vorticity_laplacian.assign(xs.size(), 0.0);
//     j_laplacian.assign(xs.size(), 0.0);

//     std::vector<double> vorticity_none_local(xs.size(), 0.0);
//     std::vector<double> j_none_local(xs.size(), 0.0);

//     evaluate_vorticity_laplacian(vorticity_laplacian, vorticity_none_local, xs, ys, u_weights, t);
//     evaluate_j_laplacian(j_laplacian, j_none_local, xs, ys, b_weights, t);


//     // calculate source terms
//     B_dot_grad_j.assign(xs.size(), 0.0);
//     for (int i = 0; i < xs.size(); ++i) {
//         B_dot_grad_j[i] = b1s[i] * j_grad_x[i] + b2s[i] * j_grad_y[i];
//     }

//     B_dot_grad_vorticity.assign(xs.size(), 0.0);
//     for (int i = 0; i < xs.size(); ++i) {
//         B_dot_grad_vorticity[i] = b1s[i] * vorticity_grad_x[i] + b2s[i] * vorticity_grad_y[i];
//     }

//     B_grad_x_dot_u2_grad.assign(xs.size(), 0.0);
//     for (int i = 0; i < xs.size(); ++i) {
//         B_grad_x_dot_u2_grad[i] = b1s_grad_x[i] * u2s_grad_x[i] + b2s_grad_x[i] * u2s_grad_y[i];
//     }

//     B_grad_y_dot_u1_grad.assign(xs.size(), 0.0);
//     for (int i = 0; i < xs.size(); ++i) {
//         B_grad_y_dot_u1_grad[i] = b1s_grad_y[i] * u1s_grad_x[i] + b2s_grad_y[i] * u1s_grad_y[i];
//     }
//     // start pushing, particle position, vorticity and current density 
//     for (int i = 0; i < xs.size(); i++) {
//         xs[i] += dt * u1s[i];
//     }

//     for (int i = 0; i < ys.size(); i++) {
//         ys[i] += dt * u2s[i];
//     }

//     for (int i = 0; i < xs.size(); i++) {
//         w0s[i] += dt * (nu * vorticity_laplacian[i] + B_dot_grad_j[i]);
//     }

//     for (int i = 0; i < xs.size(); i++) {
//         j0s[i] += dt * (mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i]);
//     }

//     return 0;
// }







int AMRStructure::euler() {
    cout << "enter euler" << endl;
    u1s.assign(xs.size(), 0.0);
    u2s.assign(xs.size(), 0.0);
    evaluate_u_field(u1s, u2s, xs, ys, u_weights, t);

    // B evaluation
    b1s.assign(xs.size(), 0.0);
    b2s.assign(xs.size(), 0.0);
    evaluate_b_field(b1s, b2s, xs, ys, b_weights, t);

    //external field for alfven wave
    for (size_t i = 0; i < b1s.size(); ++i) {
        b1s[i] +=  1.0;
    }

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


    // for each leaf panel, do centered finite difference
    for (int panel_ind = 0; panel_ind < panels.size(); panel_ind++) {
        Panel* panel = &(panels[panel_ind]);
        // only use leaf panels 
        if (panel->child_inds_start > -1) {
            continue;
        }
        // cout << "panel " << panel_ind << ": " << endl;
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

        Panel* left_panel = &(panels[panel->left_nbr_ind]);
        const int* left_panel_point_inds = left_panel->point_inds;
        double left_panel_w0s[9], left_panel_j0s[9];
        double left_panel_u1s[9], left_panel_u2s[9];
        double left_panel_b1s[9], left_panel_b2s[9];

        for (int ii = 0; ii < 9; ++ii) {
            int pind = left_panel_point_inds[ii];
            left_panel_w0s[ii] = w0s[pind];
            left_panel_j0s[ii] = j0s[pind];
            left_panel_u1s[ii] = u1s[pind];
            left_panel_u2s[ii] = u2s[pind];
            left_panel_b1s[ii] = b1s[pind];
            left_panel_b2s[ii] = b2s[pind];
        }

        Panel* right_panel = &(panels[panel->right_nbr_ind]);
        const int* right_panel_point_inds = right_panel->point_inds;
        double right_panel_w0s[9], right_panel_j0s[9];
        double right_panel_u1s[9], right_panel_u2s[9];
        double right_panel_b1s[9], right_panel_b2s[9];

        for (int ii = 0; ii < 9; ++ii) {
            int pind = right_panel_point_inds[ii];
            right_panel_w0s[ii] = w0s[pind];
            right_panel_j0s[ii] = j0s[pind];
            right_panel_u1s[ii] = u1s[pind];
            right_panel_u2s[ii] = u2s[pind];
            right_panel_b1s[ii] = b1s[pind];
            right_panel_b2s[ii] = b2s[pind];
        }

        double top_panel_w0s[9], top_panel_j0s[9];
        double top_panel_u1s[9], top_panel_u2s[9];
        double top_panel_b1s[9], top_panel_b2s[9];
        if (panel->top_nbr_ind > -2) {
            Panel* top_panel = &(panels[panel->top_nbr_ind]);
            const int* top_panel_point_inds = top_panel->point_inds;
            for (int ii = 0; ii < 9; ++ii) {
                int pind = top_panel_point_inds[ii];
                top_panel_w0s[ii] = w0s[pind];
                top_panel_j0s[ii] = j0s[pind];
                top_panel_u1s[ii] = u1s[pind];
                top_panel_u2s[ii] = u2s[pind];
                top_panel_b1s[ii] = b1s[pind];
                top_panel_b2s[ii] = b2s[pind];
            }
        }

        double bottom_panel_w0s[9], bottom_panel_j0s[9];
        double bottom_panel_u1s[9], bottom_panel_u2s[9];
        double bottom_panel_b1s[9], bottom_panel_b2s[9];
        if (panel->bottom_nbr_ind > -2) {
            Panel* bottom_panel = &(panels[panel->bottom_nbr_ind]);
            const int* bottom_panel_point_inds = bottom_panel->point_inds;
            for (int ii = 0; ii < 9; ++ii) {
                int pind = bottom_panel_point_inds[ii];
                bottom_panel_w0s[ii] = w0s[pind];
                bottom_panel_j0s[ii] = j0s[pind];
                bottom_panel_u1s[ii] = u1s[pind];
                bottom_panel_u2s[ii] = u2s[pind];
                bottom_panel_b1s[ii] = b1s[pind];
                bottom_panel_b2s[ii] = b2s[pind];
            }
        }


        vector<double> dx_j0(9, 0.0);
        vector<double> dx_w0(9, 0.0);
        vector<double> dx_u1s(9, 0.0);
        vector<double> dx_u2s(9, 0.0);
        vector<double> dx_b1s(9, 0.0);
        vector<double> dx_b2s(9, 0.0);

        vector<double> dy_j0(9, 0.0);
        vector<double> dy_w0(9, 0.0);
        vector<double> dy_u1s(9, 0.0);
        vector<double> dy_u2s(9, 0.0);
        vector<double> dy_b1s(9, 0.0);
        vector<double> dy_b2s(9, 0.0);

        /////////////// j0 /////////////
        double hx = panel_xs[3]- panel_xs[0];
        dx_j0[0] = (panel_j0s[3] - left_panel_j0s[3])/(2*hx);
        dx_j0[3] = (panel_j0s[6] - panel_j0s[0])/(2*hx);
        dx_j0[6] = (right_panel_j0s[3] - panel_j0s[3])/(2*hx);

        dx_j0[1] = (panel_j0s[4] - left_panel_j0s[4])/(2*hx);
        dx_j0[4] = (panel_j0s[7] - panel_j0s[1])/(2*hx);
        dx_j0[7] = (right_panel_j0s[4] - panel_j0s[4])/(2*hx);

        dx_j0[2] = (panel_j0s[5] - left_panel_j0s[5])/(2*hx);
        dx_j0[5] = (panel_j0s[8] - panel_j0s[2])/(2*hx);
        dx_j0[8] = (right_panel_j0s[5] - panel_j0s[5])/(2*hx);

        // std::cout << std::fixed << std::setprecision(17);
        // cout << "ys: " << endl;
        // cout << panel_ys[0] << ", " << panel_ys[3] << ", " << panel_ys[6] << endl;
        // cout << panel_ys[1] << ", " << panel_ys[4] << ", " << panel_ys[7] << endl;
        // cout << panel_ys[2] << ", " << panel_ys[5] << ", " << panel_ys[8] << endl;
        // cout << "j0s: " << endl;
        // cout << panel_j0s[0] << ", " << panel_j0s[3] << ", " << panel_j0s[6] << endl;
        // cout << panel_j0s[1] << ", " << panel_j0s[4] << ", " << panel_j0s[7] << endl;
        // cout << panel_j0s[2] << ", " << panel_j0s[5] << ", " << panel_j0s[8] << endl;
        // cout << (panel_j0s[3] - left_panel_j0s[3])<< endl;
        // cout << (panel_j0s[6] - panel_j0s[0]) << endl;
        // cout << (right_panel_j0s[3] - panel_j0s[3]) << endl;

        // cout << (panel_j0s[4] - left_panel_j0s[4]) << endl;
        // cout << (panel_j0s[7] - panel_j0s[1]) << endl;
        // cout << (right_panel_j0s[4] - panel_j0s[4]) << endl;

        // cout << (panel_j0s[5] - left_panel_j0s[5]) << endl;
        // cout << (panel_j0s[8] - panel_j0s[2]) << endl;
        // cout << (right_panel_j0s[5] - panel_j0s[5]) << endl;


        double hy = panel_ys[1]- panel_ys[0];
        if (panel->bottom_nbr_ind> -2) {
            dy_j0[0] = (panel_j0s[1] - bottom_panel_j0s[1])/(2*hy);
            dy_j0[3] = (panel_j0s[4] - bottom_panel_j0s[4])/(2*hy);
            dy_j0[6] = (panel_j0s[7] - bottom_panel_j0s[7])/(2*hy);
        }
        else {
            dy_j0[0] = (-3*(panel_j0s[0] - panel_j0s[1]) + (panel_j0s[1] - panel_j0s[2]))/(2*hy);
            dy_j0[3] = (-3*(panel_j0s[3] - panel_j0s[4]) + (panel_j0s[4] - panel_j0s[5]))/(2*hy);
            dy_j0[6] = (-3*(panel_j0s[6] - panel_j0s[7]) + (panel_j0s[7] - panel_j0s[8]))/(2*hy);
        }

        dy_j0[1] = (panel_j0s[2] - panel_j0s[0])/(2*hy);
        dy_j0[4] = (panel_j0s[5] - panel_j0s[3])/(2*hy);
        dy_j0[7] = (panel_j0s[8] - panel_j0s[6])/(2*hy);

        if (panel->top_nbr_ind> -2) {
            dy_j0[2] = (top_panel_j0s[1] - panel_j0s[1])/(2*hy);
            dy_j0[5] = (top_panel_j0s[4] - panel_j0s[4])/(2*hy);
            dy_j0[8] = (top_panel_j0s[7] - panel_j0s[7])/(2*hy);
        }
        else {
            dy_j0[2] = (-3*(panel_j0s[0] - panel_j0s[1]) + (panel_j0s[1] - panel_j0s[2]))/(2*hy);
            dy_j0[5] = (-3*(panel_j0s[3] - panel_j0s[4]) + (panel_j0s[4] - panel_j0s[5]))/(2*hy);
            dy_j0[8] = (-3*(panel_j0s[6] - panel_j0s[7]) + (panel_j0s[7] - panel_j0s[8]))/(2*hy);
        }


        //////////// w0 ///////////////
        dx_w0[0] = (panel_w0s[3] - left_panel_w0s[3])/(2*hx);
        dx_w0[3] = (panel_w0s[6] - panel_w0s[0])/(2*hx);
        dx_w0[6] = (right_panel_w0s[3] - panel_w0s[3])/(2*hx);

        dx_w0[1] = (panel_w0s[4] - left_panel_w0s[4])/(2*hx);
        dx_w0[4] = (panel_w0s[7] - panel_w0s[1])/(2*hx);
        dx_w0[7] = (right_panel_w0s[4] - panel_w0s[4])/(2*hx);

        dx_w0[2] = (panel_w0s[5] - left_panel_w0s[5])/(2*hx);
        dx_w0[5] = (panel_w0s[8] - panel_w0s[2])/(2*hx);
        dx_w0[8] = (right_panel_w0s[5] - panel_w0s[5])/(2*hx);
    
        if (panel->bottom_nbr_ind> -2) {
            dy_w0[0] = (panel_w0s[1] - bottom_panel_w0s[1])/(2*hy);
            dy_w0[3] = (panel_w0s[4] - bottom_panel_w0s[4])/(2*hy);
            dy_w0[6] = (panel_w0s[7] - bottom_panel_w0s[7])/(2*hy);
        }
        else {
            dy_w0[0] = (-3*(panel_w0s[0] - panel_w0s[1]) + (panel_w0s[1] - panel_w0s[2]))/(2*hy);
            dy_w0[3] = (-3*(panel_w0s[3] - panel_w0s[4]) + (panel_w0s[4] - panel_w0s[5]))/(2*hy);
            dy_w0[6] = (-3*(panel_w0s[6] - panel_w0s[7]) + (panel_w0s[7] - panel_w0s[8]))/(2*hy);
        }

        dy_w0[1] = (panel_w0s[2] - panel_w0s[0])/(2*hy);
        dy_w0[4] = (panel_w0s[5] - panel_w0s[3])/(2*hy);
        dy_w0[7] = (panel_w0s[8] - panel_w0s[6])/(2*hy);

        if (panel->top_nbr_ind> -2) {
            dy_w0[2] = (top_panel_w0s[1] - panel_w0s[1])/(2*hy);
            dy_w0[5] = (top_panel_w0s[4] - panel_w0s[4])/(2*hy);
            dy_w0[8] = (top_panel_w0s[7] - panel_w0s[7])/(2*hy);
        }
        else {
            dy_w0[2] = (-3*(panel_w0s[0] - panel_w0s[1]) + (panel_w0s[1] - panel_w0s[2]))/(2*hy);
            dy_w0[5] = (-3*(panel_w0s[3] - panel_w0s[4]) + (panel_w0s[4] - panel_w0s[5]))/(2*hy);
            dy_w0[8] = (-3*(panel_w0s[6] - panel_w0s[7]) + (panel_w0s[7] - panel_w0s[8]))/(2*hy);
        }

    
        
        //////////// u1s ///////////////
        dx_u1s[0] = (panel_u1s[3] - left_panel_u1s[3])/(2*hx);
        dx_u1s[3] = (panel_u1s[6] - panel_u1s[0])/(2*hx);
        dx_u1s[6] = (right_panel_u1s[3] - panel_u1s[3])/(2*hx);

        dx_u1s[1] = (panel_u1s[4] - left_panel_u1s[4])/(2*hx);
        dx_u1s[4] = (panel_u1s[7] - panel_u1s[1])/(2*hx);
        dx_u1s[7] = (right_panel_u1s[4] - panel_u1s[4])/(2*hx);

        dx_u1s[2] = (panel_u1s[5] - left_panel_u1s[5])/(2*hx);
        dx_u1s[5] = (panel_u1s[8] - panel_u1s[2])/(2*hx);
        dx_u1s[8] = (right_panel_u1s[5] - panel_u1s[5])/(2*hx);
    
        if (panel->bottom_nbr_ind> -2) {
            dy_u1s[0] = (panel_u1s[1] - bottom_panel_u1s[1])/(2*hy);
            dy_u1s[3] = (panel_u1s[4] - bottom_panel_u1s[4])/(2*hy);
            dy_u1s[6] = (panel_u1s[7] - bottom_panel_u1s[7])/(2*hy);
        }
        else {
            dy_u1s[0] = (-3*(panel_u1s[0] - panel_u1s[1]) + (panel_u1s[1] - panel_u1s[2]))/(2*hy);
            dy_u1s[3] = (-3*(panel_u1s[3] - panel_u1s[4]) + (panel_u1s[4] - panel_u1s[5]))/(2*hy);
            dy_u1s[6] = (-3*(panel_u1s[6] - panel_u1s[7]) + (panel_u1s[7] - panel_u1s[8]))/(2*hy);
        }

        dy_u1s[1] = (panel_u1s[2] - panel_u1s[0])/(2*hy);
        dy_u1s[4] = (panel_u1s[5] - panel_u1s[3])/(2*hy);
        dy_u1s[7] = (panel_u1s[8] - panel_u1s[6])/(2*hy);

        if (panel->top_nbr_ind> -2) {
            dy_u1s[2] = (top_panel_u1s[1] - panel_u1s[1])/(2*hy);
            dy_u1s[5] = (top_panel_u1s[4] - panel_u1s[4])/(2*hy);
            dy_u1s[8] = (top_panel_u1s[7] - panel_u1s[7])/(2*hy);
        }
        else {
            dy_u1s[2] = (-3*(panel_u1s[0] - panel_u1s[1]) + (panel_u1s[1] - panel_u1s[2]))/(2*hy);
            dy_u1s[5] = (-3*(panel_u1s[3] - panel_u1s[4]) + (panel_u1s[4] - panel_u1s[5]))/(2*hy);
            dy_u1s[8] = (-3*(panel_u1s[6] - panel_u1s[7]) + (panel_u1s[7] - panel_u1s[8]))/(2*hy);
        }
        
        
        //////////// u2s ///////////////
        dx_u2s[0] = (panel_u2s[3] - left_panel_u2s[3])/(2*hx);
        dx_u2s[3] = (panel_u2s[6] - panel_u2s[0])/(2*hx);
        dx_u2s[6] = (right_panel_u2s[3] - panel_u2s[3])/(2*hx);

        dx_u2s[1] = (panel_u2s[4] - left_panel_u2s[4])/(2*hx);
        dx_u2s[4] = (panel_u2s[7] - panel_u2s[1])/(2*hx);
        dx_u2s[7] = (right_panel_u2s[4] - panel_u2s[4])/(2*hx);

        dx_u2s[2] = (panel_u2s[5] - left_panel_u2s[5])/(2*hx);
        dx_u2s[5] = (panel_u2s[8] - panel_u2s[2])/(2*hx);
        dx_u2s[8] = (right_panel_u2s[5] - panel_u2s[5])/(2*hx);
    
        if (panel->bottom_nbr_ind> -2) {
            dy_u2s[0] = (panel_u2s[1] - bottom_panel_u2s[1])/(2*hy);
            dy_u2s[3] = (panel_u2s[4] - bottom_panel_u2s[4])/(2*hy);
            dy_u2s[6] = (panel_u2s[7] - bottom_panel_u2s[7])/(2*hy);
        }
        else {
            dy_u2s[0] = (-3*(panel_u2s[0] - panel_u2s[1]) + (panel_u2s[1] - panel_u2s[2]))/(2*hy);
            dy_u2s[3] = (-3*(panel_u2s[3] - panel_u2s[4]) + (panel_u2s[4] - panel_u2s[5]))/(2*hy);
            dy_u2s[6] = (-3*(panel_u2s[6] - panel_u2s[7]) + (panel_u2s[7] - panel_u2s[8]))/(2*hy);
        }

        dy_u2s[1] = (panel_u2s[2] - panel_u2s[0])/(2*hy);
        dy_u2s[4] = (panel_u2s[5] - panel_u2s[3])/(2*hy);
        dy_u2s[7] = (panel_u2s[8] - panel_u2s[6])/(2*hy);

        if (panel->top_nbr_ind> -2) {
            dy_u2s[2] = (top_panel_u2s[1] - panel_u2s[1])/(2*hy);
            dy_u2s[5] = (top_panel_u2s[4] - panel_u2s[4])/(2*hy);
            dy_u2s[8] = (top_panel_u2s[7] - panel_u2s[7])/(2*hy);
        }
        else {
            dy_u2s[2] = (-3*(panel_u2s[0] - panel_u2s[1]) + (panel_u2s[1] - panel_u2s[2]))/(2*hy);
            dy_u2s[5] = (-3*(panel_u2s[3] - panel_u2s[4]) + (panel_u2s[4] - panel_u2s[5]))/(2*hy);
            dy_u2s[8] = (-3*(panel_u2s[6] - panel_u2s[7]) + (panel_u2s[7] - panel_u2s[8]))/(2*hy);
        }



        //////////// b1s ///////////////
        dx_b1s[0] = (panel_b1s[3] - left_panel_b1s[3])/(2*hx);
        dx_b1s[3] = (panel_b1s[6] - panel_b1s[0])/(2*hx);
        dx_b1s[6] = (right_panel_b1s[3] - panel_b1s[3])/(2*hx);

        dx_b1s[1] = (panel_b1s[4] - left_panel_b1s[4])/(2*hx);
        dx_b1s[4] = (panel_b1s[7] - panel_b1s[1])/(2*hx);
        dx_b1s[7] = (right_panel_b1s[4] - panel_b1s[4])/(2*hx);

        dx_b1s[2] = (panel_b1s[5] - left_panel_b1s[5])/(2*hx);
        dx_b1s[5] = (panel_b1s[8] - panel_b1s[2])/(2*hx);
        dx_b1s[8] = (right_panel_b1s[5] - panel_b1s[5])/(2*hx);
    
        if (panel->bottom_nbr_ind> -2) {
            dy_b1s[0] = (panel_b1s[1] - bottom_panel_b1s[1])/(2*hy);
            dy_b1s[3] = (panel_b1s[4] - bottom_panel_b1s[4])/(2*hy);
            dy_b1s[6] = (panel_b1s[7] - bottom_panel_b1s[7])/(2*hy);
        }
        else {
            dy_b1s[0] = (-3*(panel_b1s[0] - panel_b1s[1]) + (panel_b1s[1] - panel_b1s[2]))/(2*hy);
            dy_b1s[3] = (-3*(panel_b1s[3] - panel_b1s[4]) + (panel_b1s[4] - panel_b1s[5]))/(2*hy);
            dy_b1s[6] = (-3*(panel_b1s[6] - panel_b1s[7]) + (panel_b1s[7] - panel_b1s[8]))/(2*hy);
        }

        dy_b1s[1] = (panel_b1s[2] - panel_b1s[0])/(2*hy);
        dy_b1s[4] = (panel_b1s[5] - panel_b1s[3])/(2*hy);
        dy_b1s[7] = (panel_b1s[8] - panel_b1s[6])/(2*hy);

        if (panel->top_nbr_ind> -2) {
            dy_b1s[2] = (top_panel_b1s[1] - panel_b1s[1])/(2*hy);
            dy_b1s[5] = (top_panel_b1s[4] - panel_b1s[4])/(2*hy);
            dy_b1s[8] = (top_panel_b1s[7] - panel_b1s[7])/(2*hy);
        }
        else {
            dy_b1s[2] = (-3*(panel_b1s[0] - panel_b1s[1]) + (panel_b1s[1] - panel_b1s[2]))/(2*hy);
            dy_b1s[5] = (-3*(panel_b1s[3] - panel_b1s[4]) + (panel_b1s[4] - panel_b1s[5]))/(2*hy);
            dy_b1s[8] = (-3*(panel_b1s[6] - panel_b1s[7]) + (panel_b1s[7] - panel_b1s[8]))/(2*hy);
        }


        //////////// b2s ///////////////
        dx_b2s[0] = (panel_b2s[3] - left_panel_b2s[3])/(2*hx);
        dx_b2s[3] = (panel_b2s[6] - panel_b2s[0])/(2*hx);
        dx_b2s[6] = (right_panel_b2s[3] - panel_b2s[3])/(2*hx);

        dx_b2s[1] = (panel_b2s[4] - left_panel_b2s[4])/(2*hx);
        dx_b2s[4] = (panel_b2s[7] - panel_b2s[1])/(2*hx);
        dx_b2s[7] = (right_panel_b2s[4] - panel_b2s[4])/(2*hx);

        dx_b2s[2] = (panel_b2s[5] - left_panel_b2s[5])/(2*hx);
        dx_b2s[5] = (panel_b2s[8] - panel_b2s[2])/(2*hx);
        dx_b2s[8] = (right_panel_b2s[5] - panel_b2s[5])/(2*hx);
    
        if (panel->bottom_nbr_ind> -2) {
            dy_b2s[0] = (panel_b2s[1] - bottom_panel_b2s[1])/(2*hy);
            dy_b2s[3] = (panel_b2s[4] - bottom_panel_b2s[4])/(2*hy);
            dy_b2s[6] = (panel_b2s[7] - bottom_panel_b2s[7])/(2*hy);
        }
        else {
            dy_b2s[0] = (-3*(panel_b2s[0] - panel_b2s[1]) + (panel_b2s[1] - panel_b2s[2]))/(2*hy);
            dy_b2s[3] = (-3*(panel_b2s[3] - panel_b2s[4]) + (panel_b2s[4] - panel_b2s[5]))/(2*hy);
            dy_b2s[6] = (-3*(panel_b2s[6] - panel_b2s[7]) + (panel_b2s[7] - panel_b2s[8]))/(2*hy);
        }

        dy_b2s[1] = (panel_b2s[2] - panel_b2s[0])/(2*hy);
        dy_b2s[4] = (panel_b2s[5] - panel_b2s[3])/(2*hy);
        dy_b2s[7] = (panel_b2s[8] - panel_b2s[6])/(2*hy);

        if (panel->top_nbr_ind> -2) {
            dy_b2s[2] = (top_panel_b2s[1] - panel_b2s[1])/(2*hy);
            dy_b2s[5] = (top_panel_b2s[4] - panel_b2s[4])/(2*hy);
            dy_b2s[8] = (top_panel_b2s[7] - panel_b2s[7])/(2*hy);
        }
        else {
            dy_b2s[2] = (-3*(panel_b2s[0] - panel_b2s[1]) + (panel_b2s[1] - panel_b2s[2]))/(2*hy);
            dy_b2s[5] = (-3*(panel_b2s[3] - panel_b2s[4]) + (panel_b2s[4] - panel_b2s[5]))/(2*hy);
            dy_b2s[8] = (-3*(panel_b2s[6] - panel_b2s[7]) + (panel_b2s[7] - panel_b2s[8]))/(2*hy);
        }

        for (int ii = 0; ii < 9; ++ii) {
            int pind = panel_point_inds[ii];
            vorticity_grad_x[pind] = dx_w0[ii];
            j_grad_x[pind]         = dx_j0[ii];
            u1s_grad_x[pind]       = dx_u1s[ii];
            u2s_grad_x[pind]       = dx_u2s[ii];
            b1s_grad_x[pind]       = dx_b1s[ii];
            b2s_grad_x[pind]       = dx_b2s[ii];

            vorticity_grad_y[pind] = dy_w0[ii];
            j_grad_y[pind]         = dy_j0[ii];
            u1s_grad_y[pind]       = dy_u1s[ii];
            u2s_grad_y[pind]       = dy_u2s[ii]; 
            b1s_grad_y[pind]       = dy_b1s[ii];
            b2s_grad_y[pind]       = dy_b2s[ii];
        }

    }

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
        // cout << "B_dot_grad_j: " << B_dot_grad_j[i] << endl;
        // cout << " B_dot_grad_vorticity: " <<  B_dot_grad_vorticity[i] << endl;
    
    }

    return 0;
}





int AMRStructure::rk4() {
    cout << "enter rk4" << endl;
    const int xs_size = (int)xs.size();

    // Stage storage: k1, k2, k3, k4 for x, y, w and j
    std::vector<double> k1_u1s(xs_size, 0.0), k1_u2s(xs_size, 0.0), k1_dws(xs_size, 0.0), k1_djs(xs_size, 0.0);
    std::vector<double> k2_u1s(xs_size, 0.0), k2_u2s(xs_size, 0.0), k2_dws(xs_size, 0.0), k2_djs(xs_size, 0.0);
    std::vector<double> k3_u1s(xs_size, 0.0), k3_u2s(xs_size, 0.0), k3_dws(xs_size, 0.0), k3_djs(xs_size, 0.0);
    std::vector<double> k4_u1s(xs_size, 0.0), k4_u2s(xs_size, 0.0), k4_dws(xs_size, 0.0), k4_djs(xs_size, 0.0);

    // Temporary positions
    std::vector<double> xs2(xs_size), ys2(xs_size), ws2(xs_size), js2(xs_size);
    std::vector<double> xs3(xs_size), ys3(xs_size), ws3(xs_size), js3(xs_size);
    std::vector<double> xs4(xs_size), ys4(xs_size), ws4(xs_size), js4(xs_size);

    // ------------------------------------------------------------
    // Stage 1 : k1 = (u1s, u2s, k1_dws, k1_djs)
    // ------------------------------------------------------------
    
    compute_rhs_state(xs, ys, w0s, j0s, t, k1_u1s, k1_u2s, k1_dws, k1_djs);

    
    // ------------------------------------------------------------
    // Stage 2
    // xs2 = xs + dt/2 * k1_u1s
    // ys2 = ys + dt/2 * k1_u2s
    // ws2 = ws + dt/2 * k1_w
    // js2 = js + dt/2 * k1_j
    // ------------------------------------------------------------
    for (int i = 0; i < xs_size; ++i) {
        xs2[i] = xs[i] + 0.5 * dt * k1_u1s[i];
        ys2[i] = ys[i] + 0.5 * dt * k1_u2s[i];
        ws2[i] = w0s[i] + 0.5 * dt * k1_dws[i];
        js2[i] = j0s[i] + 0.5 * dt * k1_djs[i];
    }

    // rebuild u_weight and b_weight terms 
    // multiply new w/j after pushing 
    for (int i = 0; i < u_weights.size(); i++) {
        u_weights[i] = weights[i] * ws2[i];
        b_weights[i] = weights[i] * js2[i];
    }

    compute_rhs_state(xs2, ys2, ws2, js2, t + 0.5 * dt,
                            k2_u1s, k2_u2s, k2_dws, k2_djs);

    // ------------------------------------------------------------
    // Stage 3
    // ------------------------------------------------------------
    for (int i = 0; i < xs_size; ++i) {
        xs3[i] = xs[i] + 0.5 * dt * k2_u1s[i];
        ys3[i] = ys[i] + 0.5 * dt * k2_u2s[i];
        ws3[i] = w0s[i] + 0.5 * dt * k2_dws[i];
        js3[i] = j0s[i] + 0.5 * dt * k2_djs[i];
    }

    for (int i = 0; i < u_weights.size(); i++) {
        u_weights[i] = weights[i] * ws3[i];
        b_weights[i] = weights[i] * js3[i];
    }

    compute_rhs_state(xs3, ys3, ws3, js3, t + 0.5 * dt,
                            k3_u1s, k3_u2s, k3_dws, k3_djs);

    // ------------------------------------------------------------
    // Stage 4
    // ------------------------------------------------------------
    for (int i = 0; i < xs_size; ++i) {
        xs4[i] = xs[i] + dt * k3_u1s[i];
        ys4[i] = ys[i] + dt * k3_u2s[i];
        ws4[i] = w0s[i] + dt * k3_dws[i];
        js4[i] = j0s[i] + dt * k3_djs[i];
    }

    for (int i = 0; i < u_weights.size(); i++) {
        u_weights[i] = weights[i] * ws4[i];
        b_weights[i] = weights[i] * js4[i];
    }

    compute_rhs_state(xs4, ys4, ws4, js4, t + dt,
                            k4_u1s, k4_u2s, k4_dws, k4_djs);

    // ------------------------------------------------------------
    // Final RK4 update
    // x uses (k1_u1s, k2_u1s, k3_u1s, k4_u1s)
    // y uses (k1_u2s, k2_u2s, k3_u2s, k4_u2s)
    // w uses (k1_dws, k2_dws, k3_dws, k4_dws)
    // j uses (k1_djs, k2_djs, k3_djs, k4_djs)
    // ------------------------------------------------------------
    for (int i = 0; i < xs_size; ++i) {
        xs[i] += (dt / 6.0) * (k1_u1s[i] + 2.0 * k2_u1s[i] + 2.0 * k3_u1s[i] + k4_u1s[i]);
        ys[i] += (dt / 6.0) * (k1_u2s[i] + 2.0 * k2_u2s[i] + 2.0 * k3_u2s[i] + k4_u2s[i]);

        w0s[i] += (dt / 6.0) * (k1_dws[i] + 2.0 * k2_dws[i] + 2.0 * k3_dws[i] + k4_dws[i]);
        j0s[i] += (dt / 6.0) * (k1_djs[i] + 2.0 * k2_djs[i] + 2.0 * k3_djs[i] + k4_djs[i]);
    }

    return 0;
}




int AMRStructure::compute_rhs_state(
    std::vector<double>& xs_in,
    std::vector<double>& ys_in,
    std::vector<double>& w0s_in,
    std::vector<double>& j0s_in,
    double t_in,
    std::vector<double>& dxs_dt,
    std::vector<double>& dys_dt,
    std::vector<double>& dw0s_dt,
    std::vector<double>& dj0s_dt
) {
    cout << "enter compute_rhs_state" << endl;

    // ------------------------------------------------------------
    // backup current member state
    // ------------------------------------------------------------
    // std::vector<double> xs_backup  = xs;
    // std::vector<double> ys_backup  = ys;
    // std::vector<double> w0s_backup = w0s;
    // std::vector<double> j0s_backup = j0s;
    // double t_backup = t;

    // ------------------------------------------------------------
    // load stage state into member arrays
    // ------------------------------------------------------------
    // xs  = xs_in;
    // ys  = ys_in;
    // w0s = w0s_in;
    // j0s = j0s_in;
    // t   = t_in;

    // ------------------------------------------------------------
    // IMPORTANT:
    // If your vorticity_laplacian and j_laplacian depend on the
    // current stage state, recompute them HERE before using them.
    //
    // Example (replace with your actual routine if needed):
    // compute_laplacian_terms();
    // ------------------------------------------------------------

    // velocity evaluation
    u1s.assign(xs.size(), 0.0);
    u2s.assign(xs.size(), 0.0);
    evaluate_u_field(u1s, u2s, xs_in, ys_in, u_weights, t_in);

    // B evaluation
    b1s.assign(xs.size(), 0.0);
    b2s.assign(xs.size(), 0.0);
    evaluate_b_field(b1s, b2s, xs_in, ys_in, b_weights, t_in);
    for (size_t i = 0; i < b1s.size(); ++i) {
        b1s[i] += 1.0;
    }

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

    // for each leaf panel, do centered finite difference
    for (int panel_ind = 0; panel_ind < panels.size(); panel_ind++) {
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
            panel_xs[ii]  = xs[pind];
            panel_ys[ii]  = ys[pind];
            panel_w0s[ii] = w0s[pind];
            panel_j0s[ii] = j0s[pind];
            panel_u1s[ii] = u1s[pind];
            panel_u2s[ii] = u2s[pind];
            panel_b1s[ii] = b1s[pind];
            panel_b2s[ii] = b2s[pind];
        }

        Panel* left_panel = &(panels[panel->left_nbr_ind]);
        const int* left_panel_point_inds = left_panel->point_inds;
        double left_panel_w0s[9], left_panel_j0s[9];
        double left_panel_u1s[9], left_panel_u2s[9];
        double left_panel_b1s[9], left_panel_b2s[9];

        for (int ii = 0; ii < 9; ++ii) {
            int pind = left_panel_point_inds[ii];
            left_panel_w0s[ii] = w0s[pind];
            left_panel_j0s[ii] = j0s[pind];
            left_panel_u1s[ii] = u1s[pind];
            left_panel_u2s[ii] = u2s[pind];
            left_panel_b1s[ii] = b1s[pind];
            left_panel_b2s[ii] = b2s[pind];
        }

        Panel* right_panel = &(panels[panel->right_nbr_ind]);
        const int* right_panel_point_inds = right_panel->point_inds;
        double right_panel_w0s[9], right_panel_j0s[9];
        double right_panel_u1s[9], right_panel_u2s[9];
        double right_panel_b1s[9], right_panel_b2s[9];

        for (int ii = 0; ii < 9; ++ii) {
            int pind = right_panel_point_inds[ii];
            right_panel_w0s[ii] = w0s[pind];
            right_panel_j0s[ii] = j0s[pind];
            right_panel_u1s[ii] = u1s[pind];
            right_panel_u2s[ii] = u2s[pind];
            right_panel_b1s[ii] = b1s[pind];
            right_panel_b2s[ii] = b2s[pind];
        }

        double top_panel_w0s[9], top_panel_j0s[9];
        double top_panel_u1s[9], top_panel_u2s[9];
        double top_panel_b1s[9], top_panel_b2s[9];
        if (panel->top_nbr_ind > -2) {
            Panel* top_panel = &(panels[panel->top_nbr_ind]);
            const int* top_panel_point_inds = top_panel->point_inds;
            for (int ii = 0; ii < 9; ++ii) {
                int pind = top_panel_point_inds[ii];
                top_panel_w0s[ii] = w0s[pind];
                top_panel_j0s[ii] = j0s[pind];
                top_panel_u1s[ii] = u1s[pind];
                top_panel_u2s[ii] = u2s[pind];
                top_panel_b1s[ii] = b1s[pind];
                top_panel_b2s[ii] = b2s[pind];
            }
        }

        double bottom_panel_w0s[9], bottom_panel_j0s[9];
        double bottom_panel_u1s[9], bottom_panel_u2s[9];
        double bottom_panel_b1s[9], bottom_panel_b2s[9];
        if (panel->bottom_nbr_ind > -2) {
            Panel* bottom_panel = &(panels[panel->bottom_nbr_ind]);
            const int* bottom_panel_point_inds = bottom_panel->point_inds;
            for (int ii = 0; ii < 9; ++ii) {
                int pind = bottom_panel_point_inds[ii];
                bottom_panel_w0s[ii] = w0s[pind];
                bottom_panel_j0s[ii] = j0s[pind];
                bottom_panel_u1s[ii] = u1s[pind];
                bottom_panel_u2s[ii] = u2s[pind];
                bottom_panel_b1s[ii] = b1s[pind];
                bottom_panel_b2s[ii] = b2s[pind];
            }
        }

        std::vector<double> dx_j0(9, 0.0);
        std::vector<double> dx_w0(9, 0.0);
        std::vector<double> dx_u1s(9, 0.0);
        std::vector<double> dx_u2s(9, 0.0);
        std::vector<double> dx_b1s(9, 0.0);
        std::vector<double> dx_b2s(9, 0.0);

        std::vector<double> dy_j0(9, 0.0);
        std::vector<double> dy_w0(9, 0.0);
        std::vector<double> dy_u1s(9, 0.0);
        std::vector<double> dy_u2s(9, 0.0);
        std::vector<double> dy_b1s(9, 0.0);
        std::vector<double> dy_b2s(9, 0.0);

        /////////////// j0 /////////////
        double hx = panel_xs[3] - panel_xs[0];
        dx_j0[0] = (panel_j0s[3] - left_panel_j0s[3]) / (2 * hx);
        dx_j0[3] = (panel_j0s[6] - panel_j0s[0]) / (2 * hx);
        dx_j0[6] = (right_panel_j0s[3] - panel_j0s[3]) / (2 * hx);

        dx_j0[1] = (panel_j0s[4] - left_panel_j0s[4]) / (2 * hx);
        dx_j0[4] = (panel_j0s[7] - panel_j0s[1]) / (2 * hx);
        dx_j0[7] = (right_panel_j0s[4] - panel_j0s[4]) / (2 * hx);

        dx_j0[2] = (panel_j0s[5] - left_panel_j0s[5]) / (2 * hx);
        dx_j0[5] = (panel_j0s[8] - panel_j0s[2]) / (2 * hx);
        dx_j0[8] = (right_panel_j0s[5] - panel_j0s[5]) / (2 * hx);

        double hy = panel_ys[1] - panel_ys[0];
        if (panel->bottom_nbr_ind > -2) {
            dy_j0[0] = (panel_j0s[1] - bottom_panel_j0s[1]) / (2 * hy);
            dy_j0[3] = (panel_j0s[4] - bottom_panel_j0s[4]) / (2 * hy);
            dy_j0[6] = (panel_j0s[7] - bottom_panel_j0s[7]) / (2 * hy);
        }
        else {
            dy_j0[0] = (-3 * (panel_j0s[0] - panel_j0s[1]) + (panel_j0s[1] - panel_j0s[2])) / (2 * hy);
            dy_j0[3] = (-3 * (panel_j0s[3] - panel_j0s[4]) + (panel_j0s[4] - panel_j0s[5])) / (2 * hy);
            dy_j0[6] = (-3 * (panel_j0s[6] - panel_j0s[7]) + (panel_j0s[7] - panel_j0s[8])) / (2 * hy);
        }

        dy_j0[1] = (panel_j0s[2] - panel_j0s[0]) / (2 * hy);
        dy_j0[4] = (panel_j0s[5] - panel_j0s[3]) / (2 * hy);
        dy_j0[7] = (panel_j0s[8] - panel_j0s[6]) / (2 * hy);

        if (panel->top_nbr_ind > -2) {
            dy_j0[2] = (top_panel_j0s[1] - panel_j0s[1]) / (2 * hy);
            dy_j0[5] = (top_panel_j0s[4] - panel_j0s[4]) / (2 * hy);
            dy_j0[8] = (top_panel_j0s[7] - panel_j0s[7]) / (2 * hy);
        }
        else {
            dy_j0[2] = (-3 * (panel_j0s[0] - panel_j0s[1]) + (panel_j0s[1] - panel_j0s[2])) / (2 * hy);
            dy_j0[5] = (-3 * (panel_j0s[3] - panel_j0s[4]) + (panel_j0s[4] - panel_j0s[5])) / (2 * hy);
            dy_j0[8] = (-3 * (panel_j0s[6] - panel_j0s[7]) + (panel_j0s[7] - panel_j0s[8])) / (2 * hy);
        }

        //////////// w0 ///////////////
        dx_w0[0] = (panel_w0s[3] - left_panel_w0s[3]) / (2 * hx);
        dx_w0[3] = (panel_w0s[6] - panel_w0s[0]) / (2 * hx);
        dx_w0[6] = (right_panel_w0s[3] - panel_w0s[3]) / (2 * hx);

        dx_w0[1] = (panel_w0s[4] - left_panel_w0s[4]) / (2 * hx);
        dx_w0[4] = (panel_w0s[7] - panel_w0s[1]) / (2 * hx);
        dx_w0[7] = (right_panel_w0s[4] - panel_w0s[4]) / (2 * hx);

        dx_w0[2] = (panel_w0s[5] - left_panel_w0s[5]) / (2 * hx);
        dx_w0[5] = (panel_w0s[8] - panel_w0s[2]) / (2 * hx);
        dx_w0[8] = (right_panel_w0s[5] - panel_w0s[5]) / (2 * hx);

        if (panel->bottom_nbr_ind > -2) {
            dy_w0[0] = (panel_w0s[1] - bottom_panel_w0s[1]) / (2 * hy);
            dy_w0[3] = (panel_w0s[4] - bottom_panel_w0s[4]) / (2 * hy);
            dy_w0[6] = (panel_w0s[7] - bottom_panel_w0s[7]) / (2 * hy);
        }
        else {
            dy_w0[0] = (-3 * (panel_w0s[0] - panel_w0s[1]) + (panel_w0s[1] - panel_w0s[2])) / (2 * hy);
            dy_w0[3] = (-3 * (panel_w0s[3] - panel_w0s[4]) + (panel_w0s[4] - panel_w0s[5])) / (2 * hy);
            dy_w0[6] = (-3 * (panel_w0s[6] - panel_w0s[7]) + (panel_w0s[7] - panel_w0s[8])) / (2 * hy);
        }

        dy_w0[1] = (panel_w0s[2] - panel_w0s[0]) / (2 * hy);
        dy_w0[4] = (panel_w0s[5] - panel_w0s[3]) / (2 * hy);
        dy_w0[7] = (panel_w0s[8] - panel_w0s[6]) / (2 * hy);

        if (panel->top_nbr_ind > -2) {
            dy_w0[2] = (top_panel_w0s[1] - panel_w0s[1]) / (2 * hy);
            dy_w0[5] = (top_panel_w0s[4] - panel_w0s[4]) / (2 * hy);
            dy_w0[8] = (top_panel_w0s[7] - panel_w0s[7]) / (2 * hy);
        }
        else {
            dy_w0[2] = (-3 * (panel_w0s[0] - panel_w0s[1]) + (panel_w0s[1] - panel_w0s[2])) / (2 * hy);
            dy_w0[5] = (-3 * (panel_w0s[3] - panel_w0s[4]) + (panel_w0s[4] - panel_w0s[5])) / (2 * hy);
            dy_w0[8] = (-3 * (panel_w0s[6] - panel_w0s[7]) + (panel_w0s[7] - panel_w0s[8])) / (2 * hy);
        }

        //////////// u1s ///////////////
        dx_u1s[0] = (panel_u1s[3] - left_panel_u1s[3]) / (2 * hx);
        dx_u1s[3] = (panel_u1s[6] - panel_u1s[0]) / (2 * hx);
        dx_u1s[6] = (right_panel_u1s[3] - panel_u1s[3]) / (2 * hx);

        dx_u1s[1] = (panel_u1s[4] - left_panel_u1s[4]) / (2 * hx);
        dx_u1s[4] = (panel_u1s[7] - panel_u1s[1]) / (2 * hx);
        dx_u1s[7] = (right_panel_u1s[4] - panel_u1s[4]) / (2 * hx);

        dx_u1s[2] = (panel_u1s[5] - left_panel_u1s[5]) / (2 * hx);
        dx_u1s[5] = (panel_u1s[8] - panel_u1s[2]) / (2 * hx);
        dx_u1s[8] = (right_panel_u1s[5] - panel_u1s[5]) / (2 * hx);

        if (panel->bottom_nbr_ind > -2) {
            dy_u1s[0] = (panel_u1s[1] - bottom_panel_u1s[1]) / (2 * hy);
            dy_u1s[3] = (panel_u1s[4] - bottom_panel_u1s[4]) / (2 * hy);
            dy_u1s[6] = (panel_u1s[7] - bottom_panel_u1s[7]) / (2 * hy);
        }
        else {
            dy_u1s[0] = (-3 * (panel_u1s[0] - panel_u1s[1]) + (panel_u1s[1] - panel_u1s[2])) / (2 * hy);
            dy_u1s[3] = (-3 * (panel_u1s[3] - panel_u1s[4]) + (panel_u1s[4] - panel_u1s[5])) / (2 * hy);
            dy_u1s[6] = (-3 * (panel_u1s[6] - panel_u1s[7]) + (panel_u1s[7] - panel_u1s[8])) / (2 * hy);
        }

        dy_u1s[1] = (panel_u1s[2] - panel_u1s[0]) / (2 * hy);
        dy_u1s[4] = (panel_u1s[5] - panel_u1s[3]) / (2 * hy);
        dy_u1s[7] = (panel_u1s[8] - panel_u1s[6]) / (2 * hy);

        if (panel->top_nbr_ind > -2) {
            dy_u1s[2] = (top_panel_u1s[1] - panel_u1s[1]) / (2 * hy);
            dy_u1s[5] = (top_panel_u1s[4] - panel_u1s[4]) / (2 * hy);
            dy_u1s[8] = (top_panel_u1s[7] - panel_u1s[7]) / (2 * hy);
        }
        else {
            dy_u1s[2] = (-3 * (panel_u1s[0] - panel_u1s[1]) + (panel_u1s[1] - panel_u1s[2])) / (2 * hy);
            dy_u1s[5] = (-3 * (panel_u1s[3] - panel_u1s[4]) + (panel_u1s[4] - panel_u1s[5])) / (2 * hy);
            dy_u1s[8] = (-3 * (panel_u1s[6] - panel_u1s[7]) + (panel_u1s[7] - panel_u1s[8])) / (2 * hy);
        }

        //////////// u2s ///////////////
        dx_u2s[0] = (panel_u2s[3] - left_panel_u2s[3]) / (2 * hx);
        dx_u2s[3] = (panel_u2s[6] - panel_u2s[0]) / (2 * hx);
        dx_u2s[6] = (right_panel_u2s[3] - panel_u2s[3]) / (2 * hx);

        dx_u2s[1] = (panel_u2s[4] - left_panel_u2s[4]) / (2 * hx);
        dx_u2s[4] = (panel_u2s[7] - panel_u2s[1]) / (2 * hx);
        dx_u2s[7] = (right_panel_u2s[4] - panel_u2s[4]) / (2 * hx);

        dx_u2s[2] = (panel_u2s[5] - left_panel_u2s[5]) / (2 * hx);
        dx_u2s[5] = (panel_u2s[8] - panel_u2s[2]) / (2 * hx);
        dx_u2s[8] = (right_panel_u2s[5] - panel_u2s[5]) / (2 * hx);

        if (panel->bottom_nbr_ind > -2) {
            dy_u2s[0] = (panel_u2s[1] - bottom_panel_u2s[1]) / (2 * hy);
            dy_u2s[3] = (panel_u2s[4] - bottom_panel_u2s[4]) / (2 * hy);
            dy_u2s[6] = (panel_u2s[7] - bottom_panel_u2s[7]) / (2 * hy);
        }
        else {
            dy_u2s[0] = (-3 * (panel_u2s[0] - panel_u2s[1]) + (panel_u2s[1] - panel_u2s[2])) / (2 * hy);
            dy_u2s[3] = (-3 * (panel_u2s[3] - panel_u2s[4]) + (panel_u2s[4] - panel_u2s[5])) / (2 * hy);
            dy_u2s[6] = (-3 * (panel_u2s[6] - panel_u2s[7]) + (panel_u2s[7] - panel_u2s[8])) / (2 * hy);
        }

        dy_u2s[1] = (panel_u2s[2] - panel_u2s[0]) / (2 * hy);
        dy_u2s[4] = (panel_u2s[5] - panel_u2s[3]) / (2 * hy);
        dy_u2s[7] = (panel_u2s[8] - panel_u2s[6]) / (2 * hy);

        if (panel->top_nbr_ind > -2) {
            dy_u2s[2] = (top_panel_u2s[1] - panel_u2s[1]) / (2 * hy);
            dy_u2s[5] = (top_panel_u2s[4] - panel_u2s[4]) / (2 * hy);
            dy_u2s[8] = (top_panel_u2s[7] - panel_u2s[7]) / (2 * hy);
        }
        else {
            dy_u2s[2] = (-3 * (panel_u2s[0] - panel_u2s[1]) + (panel_u2s[1] - panel_u2s[2])) / (2 * hy);
            dy_u2s[5] = (-3 * (panel_u2s[3] - panel_u2s[4]) + (panel_u2s[4] - panel_u2s[5])) / (2 * hy);
            dy_u2s[8] = (-3 * (panel_u2s[6] - panel_u2s[7]) + (panel_u2s[7] - panel_u2s[8])) / (2 * hy);
        }

        //////////// b1s ///////////////
        dx_b1s[0] = (panel_b1s[3] - left_panel_b1s[3]) / (2 * hx);
        dx_b1s[3] = (panel_b1s[6] - panel_b1s[0]) / (2 * hx);
        dx_b1s[6] = (right_panel_b1s[3] - panel_b1s[3]) / (2 * hx);

        dx_b1s[1] = (panel_b1s[4] - left_panel_b1s[4]) / (2 * hx);
        dx_b1s[4] = (panel_b1s[7] - panel_b1s[1]) / (2 * hx);
        dx_b1s[7] = (right_panel_b1s[4] - panel_b1s[4]) / (2 * hx);

        dx_b1s[2] = (panel_b1s[5] - left_panel_b1s[5]) / (2 * hx);
        dx_b1s[5] = (panel_b1s[8] - panel_b1s[2]) / (2 * hx);
        dx_b1s[8] = (right_panel_b1s[5] - panel_b1s[5]) / (2 * hx);

        if (panel->bottom_nbr_ind > -2) {
            dy_b1s[0] = (panel_b1s[1] - bottom_panel_b1s[1]) / (2 * hy);
            dy_b1s[3] = (panel_b1s[4] - bottom_panel_b1s[4]) / (2 * hy);
            dy_b1s[6] = (panel_b1s[7] - bottom_panel_b1s[7]) / (2 * hy);
        }
        else {
            dy_b1s[0] = (-3 * (panel_b1s[0] - panel_b1s[1]) + (panel_b1s[1] - panel_b1s[2])) / (2 * hy);
            dy_b1s[3] = (-3 * (panel_b1s[3] - panel_b1s[4]) + (panel_b1s[4] - panel_b1s[5])) / (2 * hy);
            dy_b1s[6] = (-3 * (panel_b1s[6] - panel_b1s[7]) + (panel_b1s[7] - panel_b1s[8])) / (2 * hy);
        }

        dy_b1s[1] = (panel_b1s[2] - panel_b1s[0]) / (2 * hy);
        dy_b1s[4] = (panel_b1s[5] - panel_b1s[3]) / (2 * hy);
        dy_b1s[7] = (panel_b1s[8] - panel_b1s[6]) / (2 * hy);

        if (panel->top_nbr_ind > -2) {
            dy_b1s[2] = (top_panel_b1s[1] - panel_b1s[1]) / (2 * hy);
            dy_b1s[5] = (top_panel_b1s[4] - panel_b1s[4]) / (2 * hy);
            dy_b1s[8] = (top_panel_b1s[7] - panel_b1s[7]) / (2 * hy);
        }
        else {
            dy_b1s[2] = (-3 * (panel_b1s[0] - panel_b1s[1]) + (panel_b1s[1] - panel_b1s[2])) / (2 * hy);
            dy_b1s[5] = (-3 * (panel_b1s[3] - panel_b1s[4]) + (panel_b1s[4] - panel_b1s[5])) / (2 * hy);
            dy_b1s[8] = (-3 * (panel_b1s[6] - panel_b1s[7]) + (panel_b1s[7] - panel_b1s[8])) / (2 * hy);
        }

        //////////// b2s ///////////////
        dx_b2s[0] = (panel_b2s[3] - left_panel_b2s[3]) / (2 * hx);
        dx_b2s[3] = (panel_b2s[6] - panel_b2s[0]) / (2 * hx);
        dx_b2s[6] = (right_panel_b2s[3] - panel_b2s[3]) / (2 * hx);

        dx_b2s[1] = (panel_b2s[4] - left_panel_b2s[4]) / (2 * hx);
        dx_b2s[4] = (panel_b2s[7] - panel_b2s[1]) / (2 * hx);
        dx_b2s[7] = (right_panel_b2s[4] - panel_b2s[4]) / (2 * hx);

        dx_b2s[2] = (panel_b2s[5] - left_panel_b2s[5]) / (2 * hx);
        dx_b2s[5] = (panel_b2s[8] - panel_b2s[2]) / (2 * hx);
        dx_b2s[8] = (right_panel_b2s[5] - panel_b2s[5]) / (2 * hx);

        if (panel->bottom_nbr_ind > -2) {
            dy_b2s[0] = (panel_b2s[1] - bottom_panel_b2s[1]) / (2 * hy);
            dy_b2s[3] = (panel_b2s[4] - bottom_panel_b2s[4]) / (2 * hy);
            dy_b2s[6] = (panel_b2s[7] - bottom_panel_b2s[7]) / (2 * hy);
        }
        else {
            dy_b2s[0] = (-3 * (panel_b2s[0] - panel_b2s[1]) + (panel_b2s[1] - panel_b2s[2])) / (2 * hy);
            dy_b2s[3] = (-3 * (panel_b2s[3] - panel_b2s[4]) + (panel_b2s[4] - panel_b2s[5])) / (2 * hy);
            dy_b2s[6] = (-3 * (panel_b2s[6] - panel_b2s[7]) + (panel_b2s[7] - panel_b2s[8])) / (2 * hy);
        }

        dy_b2s[1] = (panel_b2s[2] - panel_b2s[0]) / (2 * hy);
        dy_b2s[4] = (panel_b2s[5] - panel_b2s[3]) / (2 * hy);
        dy_b2s[7] = (panel_b2s[8] - panel_b2s[6]) / (2 * hy);

        if (panel->top_nbr_ind > -2) {
            dy_b2s[2] = (top_panel_b2s[1] - panel_b2s[1]) / (2 * hy);
            dy_b2s[5] = (top_panel_b2s[4] - panel_b2s[4]) / (2 * hy);
            dy_b2s[8] = (top_panel_b2s[7] - panel_b2s[7]) / (2 * hy);
        }
        else {
            dy_b2s[2] = (-3 * (panel_b2s[0] - panel_b2s[1]) + (panel_b2s[1] - panel_b2s[2])) / (2 * hy);
            dy_b2s[5] = (-3 * (panel_b2s[3] - panel_b2s[4]) + (panel_b2s[4] - panel_b2s[5])) / (2 * hy);
            dy_b2s[8] = (-3 * (panel_b2s[6] - panel_b2s[7]) + (panel_b2s[7] - panel_b2s[8])) / (2 * hy);
        }

        for (int ii = 0; ii < 9; ++ii) {
            int pind = panel_point_inds[ii];
            vorticity_grad_x[pind] = dx_w0[ii];
            j_grad_x[pind]         = dx_j0[ii];
            u1s_grad_x[pind]       = dx_u1s[ii];
            u2s_grad_x[pind]       = dx_u2s[ii];
            b1s_grad_x[pind]       = dx_b1s[ii];
            b2s_grad_x[pind]       = dx_b2s[ii];

            vorticity_grad_y[pind] = dy_w0[ii];
            j_grad_y[pind]         = dy_j0[ii];
            u1s_grad_y[pind]       = dy_u1s[ii];
            u2s_grad_y[pind]       = dy_u2s[ii];
            b1s_grad_y[pind]       = dy_b1s[ii];
            b2s_grad_y[pind]       = dy_b2s[ii];
        }
    }

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

    // RHS
    dxs_dt.assign(xs.size(), 0.0);
    dys_dt.assign(xs.size(), 0.0);
    dw0s_dt.assign(xs.size(), 0.0);
    dj0s_dt.assign(xs.size(), 0.0);

    for (int i = 0; i < xs.size(); ++i) {
        dxs_dt[i]  = u1s[i];
        dys_dt[i]  = u2s[i];
        dw0s_dt[i] = nu * vorticity_laplacian[i] + B_dot_grad_j[i];
        dj0s_dt[i] = mu * j_laplacian[i]
                   + B_dot_grad_vorticity[i]
                   + 2.0 * B_grad_x_dot_u2_grad[i]
                   - 2.0 * B_grad_y_dot_u1_grad[i];
    }

    // ------------------------------------------------------------
    // restore original member state
    // ------------------------------------------------------------
    // xs  = xs_backup;
    // ys  = ys_backup;
    // w0s = w0s_backup;
    // j0s = j0s_backup;
    // t   = t_backup;

    return 0;
}