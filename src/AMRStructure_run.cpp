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


// int AMRStructure::euler() {
//     cout << "enter euler" << endl;
//     u1s.assign(xs.size(), 0.0);
//     u2s.assign(xs.size(), 0.0);
//     evaluate_u_field(u1s, u2s, xs, ys, u_weights, t);

//     #ifdef DEBUG
//         cout << "u1s/u2s first 5:" << endl;
//         for (int i = 0; i < std::min<int>(5, (int)u1s.size()); ++i) {
//             cout << i << " u1=" << u1s[i] << " u2=" << u2s[i] << endl;
//         }
//     #endif

//     // B evaluation
//     b1s.assign(xs.size(), 0.0);
//     b2s.assign(xs.size(), 0.0);
//     evaluate_b_field(b1s, b2s, xs, ys, b_weights, t);

//     // start interpolation
//     u1s_grad_x.assign(xs.size(), 0.0);
//     u1s_grad_y.assign(xs.size(), 0.0);
//     u2s_grad_x.assign(xs.size(), 0.0);
//     u2s_grad_y.assign(xs.size(), 0.0);
//     b1s_grad_x.assign(xs.size(), 0.0);
//     b1s_grad_y.assign(xs.size(), 0.0);
//     b2s_grad_x.assign(xs.size(), 0.0);
//     b2s_grad_y.assign(xs.size(), 0.0);

//     vorticity_grad_x.assign(xs.size(), 0.0);
//     vorticity_grad_y.assign(xs.size(), 0.0);
//     j_grad_x.assign(xs.size(), 0.0);
//     j_grad_y.assign(xs.size(), 0.0);


//     // to average: 
//     std::vector<int> grad_count(xs.size(), 0);
//     for (int panel_ind = 0; panel_ind < panels.size(); panel_ind++) {
//         // interpolate_from_panel_to_points
//         Panel* panel = &(panels[panel_ind]);
//         // only use leaf panels 
//         if (panel->child_inds_start > -1) {
//             continue;
//         }
//         const int* panel_point_inds = panel->point_inds;
//         double panel_xs[9], panel_ys[9];
//         double panel_w0s[9], panel_j0s[9];
//         double panel_u1s[9], panel_u2s[9];
//         double panel_b1s[9], panel_b2s[9];

//         for (int ii = 0; ii < 9; ++ii) {
//             int pind = panel_point_inds[ii];
//             panel_xs[ii] = xs[pind];
//             panel_ys[ii] = ys[pind];
//             panel_w0s[ii] = w0s[pind];
//             panel_j0s[ii] = j0s[pind];
//             panel_u1s[ii] = u1s[pind];
//             panel_u2s[ii] = u2s[pind];
//             panel_b1s[ii] = b1s[pind];
//             panel_b2s[ii] = b2s[pind];
//         }

//         double panel_dx[9], panel_dy[9];
//         for (int ii = 0; ii < 9; ++ii) {
//             panel_dx[ii] = panel_xs[ii] - panel_xs[4];
//             panel_dy[ii] = panel_ys[ii] - panel_ys[4];
//         }

//         vector<double> dx_j0(9, 0.0);
//         vector<double> dx_w0(9, 0.0);
//         vector<double> dx_u1s(9, 0.0);
//         vector<double> dx_u2s(9, 0.0);
//         vector<double> dx_b1s(9, 0.0);
//         vector<double> dx_b2s(9, 0.0);

//         vector<double> dy_j0(9, 0.0);
//         vector<double> dy_w0(9, 0.0);
//         vector<double> dy_u1s(9, 0.0);
//         vector<double> dy_u2s(9, 0.0);
//         vector<double> dy_b1s(9, 0.0);
//         vector<double> dy_b2s(9, 0.0);

//         /////////////// j0 /////////////
//         double hx = panel_xs[3]- panel_xs[0];
//         // dx_j0[0] = (-3*(panel_j0s[0] - panel_j0s[3]) + (panel_j0s[3] - panel_j0s[6]))/(2*hx);
//         dx_j0[0] = (panel_j0s[3] - panel_j0s[0])/(hx);
//         dx_j0[3] = (panel_j0s[6] - panel_j0s[0])/(2*hx);
//         // dx_j0[6] = (3*(panel_j0s[6] - panel_j0s[3]) + (panel_j0s[0] - panel_j0s[3]))/(2*hx);
//         dx_j0[6] = (panel_j0s[6] - panel_j0s[3])/(hx);

//         // dx_j0[1] = (-3*(panel_j0s[1] - panel_j0s[4]) + (panel_j0s[4] - panel_j0s[7]))/(2*hx);
//         dx_j0[1] = (panel_j0s[4] - panel_j0s[1])/(hx);
//         dx_j0[4] = (panel_j0s[7] - panel_j0s[1])/(2*hx);
//         // dx_j0[7] = (3*(panel_j0s[7] - panel_j0s[4]) + (panel_j0s[1] - panel_j0s[4]))/(2*hx);
//         dx_j0[7] = (panel_j0s[7] - panel_j0s[4])/(hx);

//         // dx_j0[2] = (-3*(panel_j0s[2] - panel_j0s[5]) + (panel_j0s[5] - panel_j0s[8]))/(2*hx);
//         dx_j0[2] = (panel_j0s[5] - panel_j0s[2])/(hx);
//         dx_j0[5] = (panel_j0s[8] - panel_j0s[2])/(2*hx);
//         // dx_j0[8] = (3*(panel_j0s[8] - panel_j0s[5]) + (panel_j0s[2] - panel_j0s[5]))/(2*hx);
//         dx_j0[8] = (panel_j0s[8] - panel_j0s[5])/(hx);

//         cout << panel_j0s[0] << ", " << panel_j0s[3] << ", " << panel_j0s[6] << endl;
//         cout << panel_j0s[1] << ", " << panel_j0s[4] << ", " << panel_j0s[7] << endl;
//         cout << panel_j0s[2] << ", " << panel_j0s[5] << ", " << panel_j0s[8] << endl;
//         cout << (panel_j0s[3] - panel_j0s[0])<< endl;
//         cout << (panel_j0s[6] - panel_j0s[0]) << endl;
//         cout << (panel_j0s[6] - panel_j0s[3]) << endl;

//         cout << (panel_j0s[4] - panel_j0s[1]) << endl;
//         cout << (panel_j0s[7] - panel_j0s[1]) << endl;
//         cout << (panel_j0s[7] - panel_j0s[4]) << endl;

//         cout << (panel_j0s[5] - panel_j0s[2]) << endl;
//         cout << (panel_j0s[8] - panel_j0s[2]) << endl;
//         cout << (panel_j0s[8] - panel_j0s[5]) << endl;


//         double hy = panel_ys[1]- panel_ys[0];
//         // dy_j0[0] = (-3*(panel_j0s[0] - panel_j0s[1]) + (panel_j0s[1] - panel_j0s[2]))/(2*hy);
//         // dy_j0[1] = (panel_j0s[2] - panel_j0s[0])/(2*hy);
//         // dy_j0[2] = (3*(panel_j0s[2] - panel_j0s[1]) + (panel_j0s[0] - panel_j0s[1]))/(2*hy);

//         // dy_j0[3] = (-3*(panel_j0s[3] - panel_j0s[4]) + (panel_j0s[4] - panel_j0s[5]))/(2*hy);
//         // dy_j0[4] = (panel_j0s[5] - panel_j0s[3])/(2*hy);
//         // dy_j0[5] = (3*(panel_j0s[5] - panel_j0s[4]) + (panel_j0s[3] - panel_j0s[4]))/(2*hy);

//         // dy_j0[6] = (-3*(panel_j0s[6] - panel_j0s[7]) + (panel_j0s[7] - panel_j0s[8]))/(2*hy);
//         // dy_j0[7] = (panel_j0s[8] - panel_j0s[6])/(2*hy);
//         // dy_j0[8] = (3*(panel_j0s[8] - panel_j0s[7]) + (panel_j0s[6] - panel_j0s[7]))/(2*hy);

//         dy_j0[0] = (panel_j0s[1] - panel_j0s[0])/(hy);
//         dy_j0[1] = (panel_j0s[2] - panel_j0s[0])/(2*hy);
//         dy_j0[2] = (panel_j0s[2] - panel_j0s[1])/(hy);

//         dy_j0[3] = (panel_j0s[4] - panel_j0s[3])/(hy);
//         dy_j0[4] = (panel_j0s[5] - panel_j0s[3])/(2*hy);
//         dy_j0[5] = (panel_j0s[5] - panel_j0s[4])/(hy);

//         dy_j0[6] = (panel_j0s[7] - panel_j0s[6])/(hy);
//         dy_j0[7] = (panel_j0s[8] - panel_j0s[6])/(2*hy);
//         dy_j0[8] = (panel_j0s[7] - panel_j0s[6])/(hy);

//         //////////// w0 ///////////////
//         // dx_w0[0] = (-3*(panel_w0s[0] - panel_w0s[3]) + (panel_w0s[3] - panel_w0s[6]))/(2*hx);
//         // dx_w0[3] = (panel_w0s[6] - panel_w0s[0])/(2*hx);
//         // dx_w0[6] = (3*(panel_w0s[6] - panel_w0s[3]) + (panel_w0s[0] - panel_w0s[3]))/(2*hx);

//         // dx_w0[1] = (-3*(panel_w0s[1] - panel_w0s[4]) + (panel_w0s[4] - panel_w0s[7]))/(2*hx);
//         // dx_w0[4] = (panel_w0s[7] - panel_w0s[1])/(2*hx);
//         // dx_w0[7] = (3*(panel_w0s[7] - panel_w0s[4]) + (panel_w0s[1] - panel_w0s[4]))/(2*hx);

//         // dx_w0[2] = (-3*(panel_w0s[2] - panel_w0s[5]) + (panel_w0s[5] - panel_w0s[8]))/(2*hx);
//         // dx_w0[5] = (panel_w0s[8] - panel_w0s[2])/(2*hx);
//         // dx_w0[8] = (3*(panel_w0s[8] - panel_w0s[5]) + (panel_w0s[2] - panel_w0s[5]))/(2*hx);

        
//         // dy_w0[0] = (-3*(panel_w0s[0] - panel_w0s[1]) + (panel_w0s[1] - panel_w0s[2]))/(2*hy);
//         // dy_w0[1] = (panel_w0s[2] - panel_w0s[0])/(2*hy);
//         // dy_w0[2] = (3*(panel_w0s[2] - panel_w0s[1]) + (panel_w0s[0] - panel_w0s[1]))/(2*hy);

//         // dy_w0[3] = (-3*(panel_w0s[3] - panel_w0s[4]) + (panel_w0s[4] - panel_w0s[5]))/(2*hy);
//         // dy_w0[4] = (panel_w0s[5] - panel_w0s[3])/(2*hy);
//         // dy_w0[5] = (3*(panel_w0s[5] - panel_w0s[4]) + (panel_w0s[3] - panel_w0s[4]))/(2*hy);

//         // dy_w0[6] = (-3*(panel_w0s[6] - panel_w0s[7]) + (panel_w0s[7] - panel_w0s[8]))/(2*hy);
//         // dy_w0[7] = (panel_w0s[8] - panel_w0s[6])/(2*hy);
//         // dy_w0[8] = (3*(panel_w0s[8] - panel_w0s[7]) + (panel_w0s[6] - panel_w0s[7]))/(2*hy);
    
        
//         dx_w0[0] = (panel_w0s[3] - panel_w0s[0])/(hx);
//         dx_w0[3] = (panel_w0s[6] - panel_w0s[0])/(2*hx);
//         dx_w0[6] = (panel_w0s[6] - panel_w0s[3])/(hx);

//         dx_w0[1] = (panel_w0s[4] - panel_w0s[1])/(hx);
//         dx_w0[4] = (panel_w0s[7] - panel_w0s[1])/(2*hx);
//         dx_w0[7] = (panel_w0s[7] - panel_w0s[4])/(hx);

//         dx_w0[2] = (panel_w0s[5] - panel_w0s[2])/(hx);
//         dx_w0[5] = (panel_w0s[8] - panel_w0s[2])/(2*hx);
//         dx_w0[8] = (panel_w0s[8] - panel_w0s[5])/(hx);

//         dy_w0[0] = (panel_w0s[1] - panel_w0s[0])/(hy);
//         dy_w0[1] = (panel_w0s[2] - panel_w0s[0])/(2*hy);
//         dy_w0[2] = (panel_w0s[2] - panel_w0s[1])/(hy);

//         dy_w0[3] = (panel_w0s[4] - panel_w0s[3])/(hy);
//         dy_w0[4] = (panel_w0s[5] - panel_w0s[3])/(2*hy);
//         dy_w0[5] = (panel_w0s[5] - panel_w0s[4])/(hy);

//         dy_w0[6] = (panel_w0s[7] - panel_w0s[6])/(hy);
//         dy_w0[7] = (panel_w0s[8] - panel_w0s[6])/(2*hy);
//         dy_w0[8] = (panel_w0s[7] - panel_w0s[6])/(hy);


    
        
//         //////////// u1s ///////////////
//         // dx_u1s[0] = (-3*(panel_u1s[0] - panel_u1s[3]) + (panel_u1s[3] - panel_u1s[6]))/(2*hx);
//         // dx_u1s[3] = (panel_u1s[6] - panel_u1s[0])/(2*hx);
//         // dx_u1s[6] = (3*(panel_u1s[6] - panel_u1s[3]) + (panel_u1s[0] - panel_u1s[3]))/(2*hx);

//         // dx_u1s[1] = (-3*(panel_u1s[1] - panel_u1s[4]) + (panel_u1s[4] - panel_u1s[7]))/(2*hx);
//         // dx_u1s[4] = (panel_u1s[7] - panel_u1s[1])/(2*hx);
//         // dx_u1s[7] = (3*(panel_u1s[7] - panel_u1s[4]) + (panel_u1s[1] - panel_u1s[4]))/(2*hx);

//         // dx_u1s[2] = (-3*(panel_u1s[2] - panel_u1s[5]) + (panel_u1s[5] - panel_u1s[8]))/(2*hx);
//         // dx_u1s[5] = (panel_u1s[8] - panel_u1s[2])/(2*hx);
//         // dx_u1s[8] = (3*(panel_u1s[8] - panel_u1s[5]) + (panel_u1s[2] - panel_u1s[5]))/(2*hx);

        
//         // dy_u1s[0] = (-3*(panel_u1s[0] - panel_u1s[1]) + (panel_u1s[1] - panel_u1s[2]))/(2*hy);
//         // dy_u1s[1] = (panel_u1s[2] - panel_u1s[0])/(2*hy);
//         // dy_u1s[2] = (3*(panel_u1s[2] - panel_u1s[1]) + (panel_u1s[0] - panel_u1s[1]))/(2*hy);

//         // dy_u1s[3] = (-3*(panel_u1s[3] - panel_u1s[4]) + (panel_u1s[4] - panel_u1s[5]))/(2*hy);
//         // dy_u1s[4] = (panel_u1s[5] - panel_u1s[3])/(2*hy);
//         // dy_u1s[5] = (3*(panel_u1s[5] - panel_u1s[4]) + (panel_u1s[3] - panel_u1s[4]))/(2*hy);

//         // dy_u1s[6] = (-3*(panel_u1s[6] - panel_u1s[7]) + (panel_u1s[7] - panel_u1s[8]))/(2*hy);
//         // dy_u1s[7] = (panel_u1s[8] - panel_u1s[6])/(2*hy);
//         // dy_u1s[8] = (3*(panel_u1s[8] - panel_u1s[7]) + (panel_u1s[6] - panel_u1s[7]))/(2*hy);
    

//         dx_u1s[0] = (panel_u1s[3] - panel_u1s[0])/(hx);
//         dx_u1s[3] = (panel_u1s[6] - panel_u1s[0])/(2*hx);
//         dx_u1s[6] = (panel_u1s[6] - panel_u1s[3])/(hx);

//         dx_u1s[1] = (panel_u1s[4] - panel_u1s[1])/(hx);
//         dx_u1s[4] = (panel_u1s[7] - panel_u1s[1])/(2*hx);
//         dx_u1s[7] = (panel_u1s[7] - panel_u1s[4])/(hx);

//         dx_u1s[2] = (panel_u1s[5] - panel_u1s[2])/(hx);
//         dx_u1s[5] = (panel_u1s[8] - panel_u1s[2])/(2*hx);
//         dx_u1s[8] = (panel_u1s[8] - panel_u1s[5])/(hx);

//         dy_u1s[0] = (panel_u1s[1] - panel_u1s[0])/(hy);
//         dy_u1s[1] = (panel_u1s[2] - panel_u1s[0])/(2*hy);
//         dy_u1s[2] = (panel_u1s[2] - panel_u1s[1])/(hy);

//         dy_u1s[3] = (panel_u1s[4] - panel_u1s[3])/(hy);
//         dy_u1s[4] = (panel_u1s[5] - panel_u1s[3])/(2*hy);
//         dy_u1s[5] = (panel_u1s[5] - panel_u1s[4])/(hy);

//         dy_u1s[6] = (panel_u1s[7] - panel_u1s[6])/(hy);
//         dy_u1s[7] = (panel_u1s[8] - panel_u1s[6])/(2*hy);
//         dy_u1s[8] = (panel_u1s[7] - panel_u1s[6])/(hy);
        
        
//         //////////// u2s ///////////////
//         // dx_u2s[0] = (-3*(panel_u2s[0] - panel_u2s[3]) + (panel_u2s[3] - panel_u2s[6]))/(2*hx);
//         // dx_u2s[3] = (panel_u2s[6] - panel_u2s[0])/(2*hx);
//         // dx_u2s[6] = (3*(panel_u2s[6] - panel_u2s[3]) + (panel_u2s[0] - panel_u2s[3]))/(2*hx);

//         // dx_u2s[1] = (-3*(panel_u2s[1] - panel_u2s[4]) + (panel_u2s[4] - panel_u2s[7]))/(2*hx);
//         // dx_u2s[4] = (panel_u2s[7] - panel_u2s[1])/(2*hx);
//         // dx_u2s[7] = (3*(panel_u2s[7] - panel_u2s[4]) + (panel_u2s[1] - panel_u2s[4]))/(2*hx);

//         // dx_u2s[2] = (-3*(panel_u2s[2] - panel_u2s[5]) + (panel_u2s[5] - panel_u2s[8]))/(2*hx);
//         // dx_u2s[5] = (panel_u2s[8] - panel_u2s[2])/(2*hx);
//         // dx_u2s[8] = (3*(panel_u2s[8] - panel_u2s[5]) + (panel_u2s[2] - panel_u2s[5]))/(2*hx);

        
//         // dy_u2s[0] = (-3*(panel_u2s[0] - panel_u2s[1]) + (panel_u2s[1] - panel_u2s[2]))/(2*hy);
//         // dy_u2s[1] = (panel_u2s[2] - panel_u2s[0])/(2*hy);
//         // dy_u2s[2] = (3*(panel_u2s[2] - panel_u2s[1]) + (panel_u2s[0] - panel_u2s[1]))/(2*hy);

//         // dy_u2s[3] = (-3*(panel_u2s[3] - panel_u2s[4]) + (panel_u2s[4] - panel_u2s[5]))/(2*hy);
//         // dy_u2s[4] = (panel_u2s[5] - panel_u2s[3])/(2*hy);
//         // dy_u2s[5] = (3*(panel_u2s[5] - panel_u2s[4]) + (panel_u2s[3] - panel_u2s[4]))/(2*hy);

//         // dy_u2s[6] = (-3*(panel_u2s[6] - panel_u2s[7]) + (panel_u2s[7] - panel_u2s[8]))/(2*hy);
//         // dy_u2s[7] = (panel_u2s[8] - panel_u2s[6])/(2*hy);
//         // dy_u2s[8] = (3*(panel_u2s[8] - panel_u2s[7]) + (panel_u2s[6] - panel_u2s[7]))/(2*hy);
    

//         dx_u2s[0] = (panel_u2s[3] - panel_u2s[0])/(hx);
//         dx_u2s[3] = (panel_u2s[6] - panel_u2s[0])/(2*hx);
//         dx_u2s[6] = (panel_u2s[6] - panel_u2s[3])/(hx);

//         dx_u2s[1] = (panel_u2s[4] - panel_u2s[1])/(hx);
//         dx_u2s[4] = (panel_u2s[7] - panel_u2s[1])/(2*hx);
//         dx_u2s[7] = (panel_u2s[7] - panel_u2s[4])/(hx);

//         dx_u2s[2] = (panel_u2s[5] - panel_u2s[2])/(hx);
//         dx_u2s[5] = (panel_u2s[8] - panel_u2s[2])/(2*hx);
//         dx_u2s[8] = (panel_u2s[8] - panel_u2s[5])/(hx);

//         dy_u2s[0] = (panel_u2s[1] - panel_u2s[0])/(hy);
//         dy_u2s[1] = (panel_u2s[2] - panel_u2s[0])/(2*hy);
//         dy_u2s[2] = (panel_u2s[2] - panel_u2s[1])/(hy);

//         dy_u2s[3] = (panel_u2s[4] - panel_u2s[3])/(hy);
//         dy_u2s[4] = (panel_u2s[5] - panel_u2s[3])/(2*hy);
//         dy_u2s[5] = (panel_u2s[5] - panel_u2s[4])/(hy);

//         dy_u2s[6] = (panel_u2s[7] - panel_u2s[6])/(hy);
//         dy_u2s[7] = (panel_u2s[8] - panel_u2s[6])/(2*hy);
//         dy_u2s[8] = (panel_u2s[7] - panel_u2s[6])/(hy);




//         //////////// b1s ///////////////
//         // dx_b1s[0] = (-3*(panel_b1s[0] - panel_b1s[3]) + (panel_b1s[3] - panel_b1s[6]))/(2*hx);
//         // dx_b1s[3] = (panel_b1s[6] - panel_b1s[0])/(2*hx);
//         // dx_b1s[6] = (3*(panel_b1s[6] - panel_b1s[3]) + (panel_b1s[0] - panel_b1s[3]))/(2*hx);

//         // dx_b1s[1] = (-3*(panel_b1s[1] - panel_b1s[4]) + (panel_b1s[4] - panel_b1s[7]))/(2*hx);
//         // dx_b1s[4] = (panel_b1s[7] - panel_b1s[1])/(2*hx);
//         // dx_b1s[7] = (3*(panel_b1s[7] - panel_b1s[4]) + (panel_b1s[1] - panel_b1s[4]))/(2*hx);

//         // dx_b1s[2] = (-3*(panel_b1s[2] - panel_b1s[5]) + (panel_b1s[5] - panel_b1s[8]))/(2*hx);
//         // dx_b1s[5] = (panel_b1s[8] - panel_b1s[2])/(2*hx);
//         // dx_b1s[8] = (3*(panel_b1s[8] - panel_b1s[5]) + (panel_b1s[2] - panel_b1s[5]))/(2*hx);

        
//         // dy_b1s[0] = (-3*(panel_b1s[0] - panel_b1s[1]) + (panel_b1s[1] - panel_b1s[2]))/(2*hy);
//         // dy_b1s[1] = (panel_b1s[2] - panel_b1s[0])/(2*hy);
//         // dy_b1s[2] = (3*(panel_b1s[2] - panel_b1s[1]) + (panel_b1s[0] - panel_b1s[1]))/(2*hy);

//         // dy_b1s[3] = (-3*(panel_b1s[3] - panel_b1s[4]) + (panel_b1s[4] - panel_b1s[5]))/(2*hy);
//         // dy_b1s[4] = (panel_b1s[5] - panel_b1s[3])/(2*hy);
//         // dy_b1s[5] = (3*(panel_b1s[5] - panel_b1s[4]) + (panel_b1s[3] - panel_b1s[4]))/(2*hy);

//         // dy_b1s[6] = (-3*(panel_b1s[6] - panel_b1s[7]) + (panel_b1s[7] - panel_b1s[8]))/(2*hy);
//         // dy_b1s[7] = (panel_b1s[8] - panel_b1s[6])/(2*hy);
//         // dy_b1s[8] = (3*(panel_b1s[8] - panel_b1s[7]) + (panel_b1s[6] - panel_b1s[7]))/(2*hy);
    

//         dx_b1s[0] = (panel_b1s[3] - panel_b1s[0])/(hx);
//         dx_b1s[3] = (panel_b1s[6] - panel_b1s[0])/(2*hx);
//         dx_b1s[6] = (panel_b1s[6] - panel_b1s[3])/(hx);

//         dx_b1s[1] = (panel_b1s[4] - panel_b1s[1])/(hx);
//         dx_b1s[4] = (panel_b1s[7] - panel_b1s[1])/(2*hx);
//         dx_b1s[7] = (panel_b1s[7] - panel_b1s[4])/(hx);

//         dx_b1s[2] = (panel_b1s[5] - panel_b1s[2])/(hx);
//         dx_b1s[5] = (panel_b1s[8] - panel_b1s[2])/(2*hx);
//         dx_b1s[8] = (panel_b1s[8] - panel_b1s[5])/(hx);

//         dy_b1s[0] = (panel_b1s[1] - panel_b1s[0])/(hy);
//         dy_b1s[1] = (panel_b1s[2] - panel_b1s[0])/(2*hy);
//         dy_b1s[2] = (panel_b1s[2] - panel_b1s[1])/(hy);

//         dy_b1s[3] = (panel_b1s[4] - panel_b1s[3])/(hy);
//         dy_b1s[4] = (panel_b1s[5] - panel_b1s[3])/(2*hy);
//         dy_b1s[5] = (panel_b1s[5] - panel_b1s[4])/(hy);

//         dy_b1s[6] = (panel_b1s[7] - panel_b1s[6])/(hy);
//         dy_b1s[7] = (panel_b1s[8] - panel_b1s[6])/(2*hy);
//         dy_b1s[8] = (panel_b1s[7] - panel_b1s[6])/(hy);


//         //////////// b2s ///////////////
//         // dx_b2s[0] = (-3*(panel_b2s[0] - panel_b2s[3]) + (panel_b2s[3] - panel_b2s[6]))/(2*hx);
//         // dx_b2s[3] = (panel_b2s[6] - panel_b2s[0])/(2*hx);
//         // dx_b2s[6] = (3*(panel_b2s[6] - panel_b2s[3]) + (panel_b2s[0] - panel_b2s[3]))/(2*hx);

//         // dx_b2s[1] = (-3*(panel_b2s[1] - panel_b2s[4]) + (panel_b2s[4] - panel_b2s[7]))/(2*hx);
//         // dx_b2s[4] = (panel_b2s[7] - panel_b2s[1])/(2*hx);
//         // dx_b2s[7] = (3*(panel_b2s[7] - panel_b2s[4]) + (panel_b2s[1] - panel_b2s[4]))/(2*hx);

//         // dx_b2s[2] = (-3*(panel_b2s[2] - panel_b2s[5]) + (panel_b2s[5] - panel_b2s[8]))/(2*hx);
//         // dx_b2s[5] = (panel_b2s[8] - panel_b2s[2])/(2*hx);
//         // dx_b2s[8] = (3*(panel_b2s[8] - panel_b2s[5]) + (panel_b2s[2] - panel_b2s[5]))/(2*hx);

        
//         // dy_b2s[0] = (-3*(panel_b2s[0] - panel_b2s[1]) + (panel_b2s[1] - panel_b2s[2]))/(2*hy);
//         // dy_b2s[1] = (panel_b2s[2] - panel_b2s[0])/(2*hy);
//         // dy_b2s[2] = (3*(panel_b2s[2] - panel_b2s[1]) + (panel_b2s[0] - panel_b2s[1]))/(2*hy);

//         // dy_b2s[3] = (-3*(panel_b2s[3] - panel_b2s[4]) + (panel_b2s[4] - panel_b2s[5]))/(2*hy);
//         // dy_b2s[4] = (panel_b2s[5] - panel_b2s[3])/(2*hy);
//         // dy_b2s[5] = (3*(panel_b2s[5] - panel_b2s[4]) + (panel_b2s[3] - panel_b2s[4]))/(2*hy);

//         // dy_b2s[6] = (-3*(panel_b2s[6] - panel_b2s[7]) + (panel_b2s[7] - panel_b2s[8]))/(2*hy);
//         // dy_b2s[7] = (panel_b2s[8] - panel_b2s[6])/(2*hy);
//         // dy_b2s[8] = (3*(panel_b2s[8] - panel_b2s[7]) + (panel_b2s[6] - panel_b2s[7]))/(2*hy);
    

//         dx_b2s[0] = (panel_b2s[3] - panel_b2s[0])/(hx);
//         dx_b2s[3] = (panel_b2s[6] - panel_b2s[0])/(2*hx);
//         dx_b2s[6] = (panel_b2s[6] - panel_b2s[3])/(hx);

//         dx_b2s[1] = (panel_b2s[4] - panel_b2s[1])/(hx);
//         dx_b2s[4] = (panel_b2s[7] - panel_b2s[1])/(2*hx);
//         dx_b2s[7] = (panel_b2s[7] - panel_b2s[4])/(hx);

//         dx_b2s[2] = (panel_b2s[5] - panel_b2s[2])/(hx);
//         dx_b2s[5] = (panel_b2s[8] - panel_b2s[2])/(2*hx);
//         dx_b2s[8] = (panel_b2s[8] - panel_b2s[5])/(hx);

//         dy_b2s[0] = (panel_b2s[1] - panel_b2s[0])/(hy);
//         dy_b2s[1] = (panel_b2s[2] - panel_b2s[0])/(2*hy);
//         dy_b2s[2] = (panel_b2s[2] - panel_b2s[1])/(hy);

//         dy_b2s[3] = (panel_b2s[4] - panel_b2s[3])/(hy);
//         dy_b2s[4] = (panel_b2s[5] - panel_b2s[3])/(2*hy);
//         dy_b2s[5] = (panel_b2s[5] - panel_b2s[4])/(hy);

//         dy_b2s[6] = (panel_b2s[7] - panel_b2s[6])/(hy);
//         dy_b2s[7] = (panel_b2s[8] - panel_b2s[6])/(2*hy);
//         dy_b2s[8] = (panel_b2s[7] - panel_b2s[6])/(hy);



//         for (int ii = 0; ii < 9; ++ii) {
//             int pind = panel_point_inds[ii];
//             vorticity_grad_x[pind] += dx_w0[ii];
//             j_grad_x[pind]         += dx_j0[ii];
//             u1s_grad_x[pind]       += dx_u1s[ii];
//             u2s_grad_x[pind]       += dx_u2s[ii];
//             b1s_grad_x[pind]       += dx_b1s[ii];
//             b2s_grad_x[pind]       += dx_b2s[ii];

//             vorticity_grad_y[pind] += dy_w0[ii];
//             j_grad_y[pind]         += dy_j0[ii];
//             u1s_grad_y[pind]       += dy_u1s[ii];
//             u2s_grad_y[pind]       += dy_u2s[ii]; 
//             b1s_grad_y[pind]       += dy_b1s[ii];
//             b2s_grad_y[pind]       += dy_b2s[ii];

//             grad_count[pind] += 1;
//         }


//         // print out the gradient for each panel 
//         // cout << "leaf panel: " << panel->panel_ind <<endl;
//         // for (int ii = 0; ii < 9; ++ii) {
//         //     int pind = panel_point_inds[ii];
//         //     cout << "i=" << ii
//         //     << " x=" << xs[pind]
//         //     << " y=" << ys[pind]
//         //     << " b1=" << b1s[pind]
//         //     << " b2=" << b2s[pind]
//         //     << " j=" << j0s[pind]
//         //     // << " w=" << w0s[pind]
//         //     // << " w_dx=" << dx_w0[ii]
//         //     // << " w_dy=" << dy_w0[ii]

//         //     // << " u1_dx=" << dx_u1s[ii]
//         //     // << " u1_dy=" << dy_u1s[ii]
//         //     // << " u2_dx=" << dx_u2s[ii]
//         //     // << " u2_dy=" << dy_u2s[ii]

//         //     // << " b1_dx=" << dx_u1s[ii]
//         //     // << " b1_dy=" << dy_u1s[ii]
//         //     // << " b2_dx=" << dx_u2s[ii]
//         //     // << " b2_dy=" << dy_u2s[ii]
//         //     << " j_dx=" << dx_j0[ii]
//         //     << " j_dy=" << dy_j0[ii]
//         //     << "\n";
//         // }

//     }

//     // average gradients
//     for (int i = 0; i < xs.size(); ++i) {
//         if (grad_count[i] > 0) {
//             double inv = 1.0 / grad_count[i];
//             vorticity_grad_x[i] *= inv; vorticity_grad_y[i] *= inv;
//             j_grad_x[i] *= inv;         j_grad_y[i] *= inv;
//             u1s_grad_x[i] *= inv;       u1s_grad_y[i] *= inv;
//             u2s_grad_x[i] *= inv;       u2s_grad_y[i] *= inv;
//             b1s_grad_x[i] *= inv;       b1s_grad_y[i] *= inv;
//             b2s_grad_x[i] *= inv;       b2s_grad_y[i] *= inv;
//         }
//     }

//     // treat boundary
//     for (int panel_ind = 0; panel_ind < panels.size(); panel_ind++) {
//         Panel* panel = &(panels[panel_ind]);
//         // only use leaf panels 
//         if (panel->child_inds_start > -1) {
//             continue;
//         }
//         if (panel->is_left_bdry) {
//             // cout << " left bdry panel: " << panel->panel_ind << endl;
//             const int* panel_point_inds = panel->point_inds;
//             Panel* left_panel = &(panels[panel->left_nbr_ind]);
//             const int* left_panel_point_inds = left_panel->point_inds;
//             // three left boundary points of the panel, point index in the panel is 0,1,2
//             // three right boundary points of the left panel, point index in the panel is 6,7,8
//             int pind_0 = panel_point_inds[0];
//             int pind_1 = panel_point_inds[1];
//             int pind_2 = panel_point_inds[2];

//             int left_pind_6 = left_panel_point_inds[6];
//             int left_pind_7 = left_panel_point_inds[7];
//             int left_pind_8 = left_panel_point_inds[8];

//             vorticity_grad_x[pind_0] = (vorticity_grad_x[pind_0] + vorticity_grad_x[left_pind_6])/2;
//             vorticity_grad_y[pind_0] = (vorticity_grad_y[pind_0] + vorticity_grad_y[left_pind_6])/2;
//             j_grad_x[pind_0] = (j_grad_x[pind_0] + j_grad_x[left_pind_6])/2;
//             j_grad_y[pind_0] = (j_grad_y[pind_0] + j_grad_y[left_pind_6])/2;
//             u1s_grad_x[pind_0] = (u1s_grad_x[pind_0] + u1s_grad_x[left_pind_6])/2;
//             u1s_grad_y[pind_0] = (u1s_grad_y[pind_0] + u1s_grad_y[left_pind_6])/2;
//             u2s_grad_x[pind_0] = (u2s_grad_x[pind_0] + u2s_grad_x[left_pind_6])/2;
//             u2s_grad_y[pind_0] = (u2s_grad_y[pind_0] + u2s_grad_y[left_pind_6])/2;
//             b1s_grad_x[pind_0] = (b1s_grad_x[pind_0] + b1s_grad_x[left_pind_6])/2;
//             b1s_grad_y[pind_0] = (b1s_grad_y[pind_0] + b1s_grad_y[left_pind_6])/2;
//             b2s_grad_x[pind_0] = (b2s_grad_x[pind_0] + b2s_grad_x[left_pind_6])/2;
//             b2s_grad_y[pind_0] = (b2s_grad_y[pind_0] + b2s_grad_y[left_pind_6])/2;

//             vorticity_grad_x[pind_1] = (vorticity_grad_x[pind_1] + vorticity_grad_x[left_pind_7])/2;
//             vorticity_grad_y[pind_1] = (vorticity_grad_y[pind_1] + vorticity_grad_y[left_pind_7])/2;
//             j_grad_x[pind_1] = (j_grad_x[pind_1] + j_grad_x[left_pind_7])/2;
//             j_grad_y[pind_1] = (j_grad_y[pind_1] + j_grad_y[left_pind_7])/2;
//             u1s_grad_x[pind_1] = (u1s_grad_x[pind_1] + u1s_grad_x[left_pind_7])/2;
//             u1s_grad_y[pind_1] = (u1s_grad_y[pind_1] + u1s_grad_y[left_pind_7])/2;
//             u2s_grad_x[pind_1] = (u2s_grad_x[pind_1] + u2s_grad_x[left_pind_7])/2;
//             u2s_grad_y[pind_1] = (u2s_grad_y[pind_1] + u2s_grad_y[left_pind_7])/2;
//             b1s_grad_x[pind_1] = (b1s_grad_x[pind_1] + b1s_grad_x[left_pind_7])/2;
//             b1s_grad_y[pind_1] = (b1s_grad_y[pind_1] + b1s_grad_y[left_pind_7])/2;
//             b2s_grad_x[pind_1] = (b2s_grad_x[pind_1] + b2s_grad_x[left_pind_7])/2;
//             b2s_grad_y[pind_1] = (b2s_grad_y[pind_1] + b2s_grad_y[left_pind_7])/2;

//             vorticity_grad_x[pind_2] = (vorticity_grad_x[pind_2] + vorticity_grad_x[left_pind_8])/2;
//             vorticity_grad_y[pind_2] = (vorticity_grad_y[pind_2] + vorticity_grad_y[left_pind_8])/2;
//             j_grad_x[pind_2] = (j_grad_x[pind_2] + j_grad_x[left_pind_8])/2;
//             j_grad_y[pind_2] = (j_grad_y[pind_2] + j_grad_y[left_pind_8])/2;
//             u1s_grad_x[pind_2] = (u1s_grad_x[pind_2] + u1s_grad_x[left_pind_8])/2;
//             u1s_grad_y[pind_2] = (u1s_grad_y[pind_2] + u1s_grad_y[left_pind_8])/2;
//             u2s_grad_x[pind_2] = (u2s_grad_x[pind_2] + u2s_grad_x[left_pind_8])/2;
//             u2s_grad_y[pind_2] = (u2s_grad_y[pind_2] + u2s_grad_y[left_pind_8])/2;
//             b1s_grad_x[pind_2] = (b1s_grad_x[pind_2] + b1s_grad_x[left_pind_8])/2;
//             b1s_grad_y[pind_2] = (b1s_grad_y[pind_2] + b1s_grad_y[left_pind_8])/2;
//             b2s_grad_x[pind_2] = (b2s_grad_x[pind_2] + b2s_grad_x[left_pind_8])/2;
//             b2s_grad_y[pind_2] = (b2s_grad_y[pind_2] + b2s_grad_y[left_pind_8])/2;


//             // right boundary have the same value
//             vorticity_grad_x[left_pind_6] = vorticity_grad_x[pind_0];
//             vorticity_grad_y[left_pind_6] = vorticity_grad_y[pind_0];
//             j_grad_x[left_pind_6] = j_grad_x[pind_0];
//             j_grad_y[left_pind_6] = j_grad_y[pind_0];
//             u1s_grad_x[left_pind_6] = u1s_grad_x[pind_0];
//             u1s_grad_y[left_pind_6] = u1s_grad_y[pind_0];
//             u2s_grad_x[left_pind_6] = u2s_grad_x[pind_0];
//             u2s_grad_y[left_pind_6] = u2s_grad_y[pind_0];
//             b1s_grad_x[left_pind_6] = b1s_grad_x[pind_0];
//             b1s_grad_y[left_pind_6] = b1s_grad_y[pind_0];
//             b2s_grad_x[left_pind_6] = b2s_grad_x[pind_0];
//             b2s_grad_y[left_pind_6] = b2s_grad_y[pind_0];

//             vorticity_grad_x[left_pind_7] = vorticity_grad_x[pind_1];
//             vorticity_grad_y[left_pind_7] = vorticity_grad_y[pind_1];
//             j_grad_x[left_pind_7] = j_grad_x[pind_1];
//             j_grad_y[left_pind_7] = j_grad_y[pind_1];
//             u1s_grad_x[left_pind_7] = u1s_grad_x[pind_1];
//             u1s_grad_y[left_pind_7] = u1s_grad_y[pind_1];
//             u2s_grad_x[left_pind_7] = u2s_grad_x[pind_1];
//             u2s_grad_y[left_pind_7] = u2s_grad_y[pind_1];
//             b1s_grad_x[left_pind_7] = b1s_grad_x[pind_1];
//             b1s_grad_y[left_pind_7] = b1s_grad_y[pind_1];
//             b2s_grad_x[left_pind_7] = b2s_grad_x[pind_1];
//             b2s_grad_y[left_pind_7] = b2s_grad_y[pind_1];

//             vorticity_grad_x[left_pind_8] = vorticity_grad_x[pind_2];
//             vorticity_grad_y[left_pind_8] = vorticity_grad_y[pind_2];
//             j_grad_x[left_pind_8] = j_grad_x[pind_2];
//             j_grad_y[left_pind_8] = j_grad_y[pind_2];
//             u1s_grad_x[left_pind_8] = u1s_grad_x[pind_2];
//             u1s_grad_y[left_pind_8] = u1s_grad_y[pind_2];
//             u2s_grad_x[left_pind_8] = u2s_grad_x[pind_2];
//             u2s_grad_y[left_pind_8] = u2s_grad_y[pind_2];
//             b1s_grad_x[left_pind_8] = b1s_grad_x[pind_2];
//             b1s_grad_y[left_pind_8] = b1s_grad_y[pind_2];
//             b2s_grad_x[left_pind_8] = b2s_grad_x[pind_2];
//             b2s_grad_y[left_pind_8] = b2s_grad_y[pind_2];
//         }
//     }


//     // printing:
//     // for (int i = 0; i < static_cast<int>(xs.size()); ++i) {
//     //     std::cout
//     //         << "i=" << i
//     //         << " x=" << xs[i]
//     //         << " y=" << ys[i]
//     //         << " u1=" << u1s[i]
//     //         << " u2=" << u2s[i]
//     //         << " b1=" << b1s[i]
//     //         << " b2=" << b2s[i]
//     //         << " j=" << j0s[i]
//     //         << " w=" << w0s[i]
//     //         << " w_dx=" << vorticity_grad_x[i]
//     //         << " w_dy=" << vorticity_grad_y[i]
//     //         << " j_dx=" << j_grad_x[i]
//     //         << " j_dy=" << j_grad_y[i]
//     //         << "\n";
//     // }

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
//         ys[i] += dt * u2s[i];
//         w0s[i] += dt * (nu * vorticity_laplacian[i] + B_dot_grad_j[i]);
//         j0s[i] += dt * (mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i]);
//     }

//     // // Lax–Friedrichs, loop through panels 
//     // std::vector<int> average_count(xs.size(), 0);

//     // std::vector<double> w0s_average_state(xs.size(), 0);
//     // std::vector<double> j0s_average_state(xs.size(), 0);

//     // for (int panel_ind = 0; panel_ind < panels.size(); panel_ind++) {
//     //     Panel* panel = &(panels[panel_ind]);
//     //     if (panel->child_inds_start > -1) {
//     //         continue;
//     //     }
//     //     const int* panel_point_inds = panel->point_inds;
//     //     double panel_w0s[9], panel_j0s[9];

//     //     for (int ii = 0; ii < 9; ++ii) {
//     //         int pind = panel_point_inds[ii];
//     //         panel_w0s[ii] = w0s[pind];
//     //         panel_j0s[ii] = j0s[pind];
//     //     }

//     //     // go through each index 0 to 8
//     //     w0s_average_state[panel_point_inds[0]] += panel_w0s[1] + panel_w0s[3];
//     //     j0s_average_state[panel_point_inds[0]] += panel_j0s[1] + panel_j0s[3];
//     //     average_count[panel_point_inds[0]] += 2;
//     //     w0s_average_state[panel_point_inds[1]] += panel_w0s[0] + panel_w0s[2] + panel_w0s[4];
//     //     j0s_average_state[panel_point_inds[1]] += panel_j0s[0] + panel_j0s[2] + panel_j0s[4];
//     //     average_count[panel_point_inds[1]] += 3;
//     //     w0s_average_state[panel_point_inds[2]] += panel_w0s[1] + panel_w0s[5];
//     //     j0s_average_state[panel_point_inds[2]] += panel_j0s[1] + panel_j0s[5];
//     //     average_count[panel_point_inds[2]] += 2;
//     //     w0s_average_state[panel_point_inds[3]] += panel_w0s[0] + panel_w0s[4] + panel_w0s[6];
//     //     j0s_average_state[panel_point_inds[3]] += panel_j0s[0] + panel_j0s[4] + panel_j0s[6];
//     //     average_count[panel_point_inds[3]] += 3;
//     //     w0s_average_state[panel_point_inds[4]] += panel_w0s[1] + panel_w0s[3] + panel_w0s[5] + panel_w0s[7];
//     //     j0s_average_state[panel_point_inds[4]] += panel_j0s[1] + panel_j0s[3] + panel_j0s[5] + panel_j0s[7];
//     //     average_count[panel_point_inds[4]] += 4;
//     //     w0s_average_state[panel_point_inds[5]] += panel_w0s[2] + panel_w0s[4] + panel_w0s[8];
//     //     j0s_average_state[panel_point_inds[5]] += panel_j0s[2] + panel_j0s[4] + panel_j0s[8];
//     //     average_count[panel_point_inds[5]] += 3;
//     //     w0s_average_state[panel_point_inds[6]] += panel_w0s[3] + panel_w0s[7];
//     //     j0s_average_state[panel_point_inds[6]] += panel_j0s[3] + panel_j0s[7];
//     //     average_count[panel_point_inds[6]] += 2;
//     //     w0s_average_state[panel_point_inds[7]] += panel_w0s[4] + panel_w0s[6] + panel_w0s[8];
//     //     j0s_average_state[panel_point_inds[7]] += panel_j0s[4] + panel_j0s[6] + panel_j0s[8];
//     //     average_count[panel_point_inds[7]] += 3;
//     //     w0s_average_state[panel_point_inds[8]] += panel_w0s[5] + panel_w0s[7];
//     //     j0s_average_state[panel_point_inds[8]] += panel_j0s[5] + panel_j0s[7];
//     //     average_count[panel_point_inds[8]] += 2;


//     //     if (panel->is_left_bdry) {
//     //         Panel* left_panel = &(panels[panel->left_nbr_ind]);
//     //         const int* left_panel_point_inds = left_panel->point_inds;
//     //         // three left boundary points of the panel, point index in the panel is 0,1,2
//     //         // three right boundary points of the left panel, point index in the panel is 6,7,8
//     //         int pind_0 = panel_point_inds[0];
//     //         int pind_1 = panel_point_inds[1];
//     //         int pind_2 = panel_point_inds[2];

//     //         // cout << "debug here for Lax Friedrichs 5 " << endl;

//     //         int left_pind_3 = left_panel_point_inds[3];
//     //         int left_pind_4 = left_panel_point_inds[4];
//     //         int left_pind_5 = left_panel_point_inds[5];
//     //         w0s_average_state[pind_0] += w0s[left_pind_3];
//     //         w0s_average_state[pind_1] += w0s[left_pind_4];
//     //         w0s_average_state[pind_2] += w0s[left_pind_5];

//     //         j0s_average_state[pind_0] += j0s[left_pind_3];
//     //         j0s_average_state[pind_1] += j0s[left_pind_4];
//     //         j0s_average_state[pind_2] += j0s[left_pind_5];

//     //         average_count[pind_0] += 1;
//     //         average_count[pind_1] += 1;
//     //         average_count[pind_2] += 1;
//     //     }

//     //     if (panel->is_right_bdry) {
//     //         // cout << "debug here for Lax Friedrichs 6 " << endl;
//     //         Panel* right_panel = &(panels[panel->right_nbr_ind]);
//     //         const int* right_panel_point_inds = right_panel->point_inds;
//     //         int pind_6 = panel_point_inds[6];
//     //         int pind_7 = panel_point_inds[7];
//     //         int pind_8 = panel_point_inds[8];

//     //         // cout << "debug here for Lax Friedrichs 7 " << endl;

//     //         int right_pind_3 = right_panel_point_inds[3];
//     //         int right_pind_4 = right_panel_point_inds[4];
//     //         int right_pind_5 = right_panel_point_inds[5];
//     //         w0s_average_state[pind_6] += w0s[right_pind_3];
//     //         w0s_average_state[pind_7] += w0s[right_pind_4];
//     //         w0s_average_state[pind_8] += w0s[right_pind_5];

//     //         j0s_average_state[pind_6] += j0s[right_pind_3];
//     //         j0s_average_state[pind_7] += j0s[right_pind_4];
//     //         j0s_average_state[pind_8] += j0s[right_pind_5];

//     //         average_count[pind_6] += 1;
//     //         average_count[pind_7] += 1;
//     //         average_count[pind_8] += 1;
//     //     }
//     // }

//     // for (int i = 0; i < xs.size(); ++i) {
//     //     if (average_count[i] > 0) {
//     //         // cout << "average_count: " << average_count[i] << endl;
//     //         // cout << "w0s_average_state: " << w0s_average_state[i] << endl;
//     //         double inv = 1.0 / average_count[i];
//     //         w0s_average_state[i] *= inv;
//     //         j0s_average_state[i] *= inv;
//     //     }
//     // }

//     // // start pushing, particle position, vorticity and current density 
//     // for (int i = 0; i < xs.size(); i++) {
//     //     xs[i] += dt * u1s[i];
//     //     ys[i] += dt * u2s[i];
//     //     w0s[i] = w0s_average_state[i] + dt * (nu * vorticity_laplacian[i] + B_dot_grad_j[i]);
//     //     j0s[i] = j0s_average_state[i] + dt * (mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i]);
//     // }

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

    // // Temporary positions
    // std::vector<double> xs2(xs_size), ys2(xs_size), ws2(xs_size), js2(xs_size);
    // std::vector<double> xs3(xs_size), ys3(xs_size), ws3(xs_size), js3(xs_size);
    // std::vector<double> xs4(xs_size), ys4(xs_size), ws4(xs_size), js4(xs_size);

    // // -------------------------
    // // Stage 1 : k1 = (u1s, u2s, k1_dws, k1_djs)
    // // -------------------------
    // // U evaluation
    // u1s.assign(u1s.size(), 0.0);
    // u2s.assign(u2s.size(), 0.0);
    // evaluate_u_field(u1s, u2s, xs, ys, u_weights, t);

    // // B evaluation
    // b1s.assign(b1s.size(), 0.0);
    // b2s.assign(b2s.size(), 0.0);
    // evaluate_b_field(b1s, b2s, xs, ys, b_weights, t);

    // // solve derivative terms 
    // u1s_grad_x.assign(xs.size(), 0.0);
    // u1s_grad_y.assign(xs.size(), 0.0);
    // u2s_grad_x.assign(xs.size(), 0.0);
    // u2s_grad_y.assign(xs.size(), 0.0);
    // b1s_grad_x.assign(xs.size(), 0.0);
    // b1s_grad_y.assign(xs.size(), 0.0);
    // b2s_grad_x.assign(xs.size(), 0.0);
    // b2s_grad_y.assign(xs.size(), 0.0);
    // vorticity_grad_x.assign(xs.size(), 0.0);
    // vorticity_grad_y.assign(xs.size(), 0.0);
    // j_grad_x.assign(xs.size(), 0.0);
    // j_grad_y.assign(xs.size(), 0.0);

    // std::vector<int> grad_count(xs.size(), 0);
    // for (int panel_ind = 0; panel_ind < panels.size(); panel_ind++) {
    //     // interpolate_from_panel_to_points
    //     Panel* panel = &(panels[panel_ind]);
    //     // only use leaf panels 
    //     if (panel->child_inds_start > -1) {
    //         continue;
    //     }
    //     const int* panel_point_inds = panel->point_inds;
    //     double panel_xs[9], panel_ys[9];
    //     double panel_w0s[9], panel_j0s[9];
    //     double panel_u1s[9], panel_u2s[9];
    //     double panel_b1s[9], panel_b2s[9];

    //     for (int ii = 0; ii < 9; ++ii) {
    //         int pind = panel_point_inds[ii];
    //         panel_xs[ii] = xs[pind];
    //         panel_ys[ii] = ys[pind];
    //         panel_w0s[ii] = w0s[pind];
    //         panel_j0s[ii] = j0s[pind];
    //         panel_u1s[ii] = u1s[pind];
    //         panel_u2s[ii] = u2s[pind];
    //         panel_b1s[ii] = b1s[pind];
    //         panel_b2s[ii] = b2s[pind];
    //     }

    //     double panel_dx[9], panel_dy[9];
    //     for (int ii = 0; ii < 9; ++ii) {
    //         panel_dx[ii] = panel_xs[ii] - panel_xs[4];
    //         panel_dy[ii] = panel_ys[ii] - panel_ys[4];
    //     }

    //     vector<double> dx_j0(9, 0.0);
    //     vector<double> dx_w0(9, 0.0);
    //     vector<double> dx_u1s(9, 0.0);
    //     vector<double> dx_u2s(9, 0.0);
    //     vector<double> dx_b1s(9, 0.0);
    //     vector<double> dx_b2s(9, 0.0);

    //     vector<double> dy_j0(9, 0.0);
    //     vector<double> dy_w0(9, 0.0);
    //     vector<double> dy_u1s(9, 0.0);
    //     vector<double> dy_u2s(9, 0.0);
    //     vector<double> dy_b1s(9, 0.0);
    //     vector<double> dy_b2s(9, 0.0);

    //     /////////////// j0 /////////////
    //     double hx = panel_xs[3]- panel_xs[0];
    //     dx_j0[0] = (-3*(panel_j0s[0] - panel_j0s[3]) + (panel_j0s[3] - panel_j0s[6]))/(2*hx);
    //     dx_j0[3] = (panel_j0s[6] - panel_j0s[0])/(2*hx);
    //     dx_j0[6] = (3*(panel_j0s[6] - panel_j0s[3]) + (panel_j0s[0] - panel_j0s[3]))/(2*hx);

    //     dx_j0[1] = (-3*(panel_j0s[1] - panel_j0s[4]) + (panel_j0s[4] - panel_j0s[7]))/(2*hx);
    //     dx_j0[4] = (panel_j0s[7] - panel_j0s[1])/(2*hx);
    //     dx_j0[7] = (3*(panel_j0s[7] - panel_j0s[4]) + (panel_j0s[1] - panel_j0s[4]))/(2*hx);

    //     dx_j0[2] = (-3*(panel_j0s[2] - panel_j0s[5]) + (panel_j0s[5] - panel_j0s[8]))/(2*hx);
    //     dx_j0[5] = (panel_j0s[8] - panel_j0s[2])/(2*hx);
    //     dx_j0[8] = (3*(panel_j0s[8] - panel_j0s[5]) + (panel_j0s[2] - panel_j0s[5]))/(2*hx);

    //     double hy = panel_ys[1]- panel_ys[0];
    //     dy_j0[0] = (-3*(panel_j0s[0] - panel_j0s[1]) + (panel_j0s[1] - panel_j0s[2]))/(2*hy);
    //     dy_j0[1] = (panel_j0s[2] - panel_j0s[0])/(2*hy);
    //     dy_j0[2] = (3*(panel_j0s[2] - panel_j0s[1]) + (panel_j0s[0] - panel_j0s[1]))/(2*hy);

    //     dy_j0[3] = (-3*(panel_j0s[3] - panel_j0s[4]) + (panel_j0s[4] - panel_j0s[5]))/(2*hy);
    //     dy_j0[4] = (panel_j0s[5] - panel_j0s[3])/(2*hy);
    //     dy_j0[5] = (3*(panel_j0s[5] - panel_j0s[4]) + (panel_j0s[3] - panel_j0s[4]))/(2*hy);

    //     dy_j0[6] = (-3*(panel_j0s[6] - panel_j0s[7]) + (panel_j0s[7] - panel_j0s[8]))/(2*hy);
    //     dy_j0[7] = (panel_j0s[8] - panel_j0s[6])/(2*hy);
    //     dy_j0[8] = (3*(panel_j0s[8] - panel_j0s[7]) + (panel_j0s[6] - panel_j0s[7]))/(2*hy);



    //     //////////// w0 ///////////////
    //     dx_w0[0] = (-3*(panel_w0s[0] - panel_w0s[3]) + (panel_w0s[3] - panel_w0s[6]))/(2*hx);
    //     dx_w0[3] = (panel_w0s[6] - panel_w0s[0])/(2*hx);
    //     dx_w0[6] = (3*(panel_w0s[6] - panel_w0s[3]) + (panel_w0s[0] - panel_w0s[3]))/(2*hx);

    //     dx_w0[1] = (-3*(panel_w0s[1] - panel_w0s[4]) + (panel_w0s[4] - panel_w0s[7]))/(2*hx);
    //     dx_w0[4] = (panel_w0s[7] - panel_w0s[1])/(2*hx);
    //     dx_w0[7] = (3*(panel_w0s[7] - panel_w0s[4]) + (panel_w0s[1] - panel_w0s[4]))/(2*hx);

    //     dx_w0[2] = (-3*(panel_w0s[2] - panel_w0s[5]) + (panel_w0s[5] - panel_w0s[8]))/(2*hx);
    //     dx_w0[5] = (panel_w0s[8] - panel_w0s[2])/(2*hx);
    //     dx_w0[8] = (3*(panel_w0s[8] - panel_w0s[5]) + (panel_w0s[2] - panel_w0s[5]))/(2*hx);

        
    //     dy_w0[0] = (-3*(panel_w0s[0] - panel_w0s[1]) + (panel_w0s[1] - panel_w0s[2]))/(2*hy);
    //     dy_w0[1] = (panel_w0s[2] - panel_w0s[0])/(2*hy);
    //     dy_w0[2] = (3*(panel_w0s[2] - panel_w0s[1]) + (panel_w0s[0] - panel_w0s[1]))/(2*hy);

    //     dy_w0[3] = (-3*(panel_w0s[3] - panel_w0s[4]) + (panel_w0s[4] - panel_w0s[5]))/(2*hy);
    //     dy_w0[4] = (panel_w0s[5] - panel_w0s[3])/(2*hy);
    //     dy_w0[5] = (3*(panel_w0s[5] - panel_w0s[4]) + (panel_w0s[3] - panel_w0s[4]))/(2*hy);

    //     dy_w0[6] = (-3*(panel_w0s[6] - panel_w0s[7]) + (panel_w0s[7] - panel_w0s[8]))/(2*hy);
    //     dy_w0[7] = (panel_w0s[8] - panel_w0s[6])/(2*hy);
    //     dy_w0[8] = (3*(panel_w0s[8] - panel_w0s[7]) + (panel_w0s[6] - panel_w0s[7]))/(2*hy);
    
        
        
    
        
    //     //////////// u1s ///////////////
    //     dx_u1s[0] = (-3*(panel_u1s[0] - panel_u1s[3]) + (panel_u1s[3] - panel_u1s[6]))/(2*hx);
    //     dx_u1s[3] = (panel_u1s[6] - panel_u1s[0])/(2*hx);
    //     dx_u1s[6] = (3*(panel_u1s[6] - panel_u1s[3]) + (panel_u1s[0] - panel_u1s[3]))/(2*hx);

    //     dx_u1s[1] = (-3*(panel_u1s[1] - panel_u1s[4]) + (panel_u1s[4] - panel_u1s[7]))/(2*hx);
    //     dx_u1s[4] = (panel_u1s[7] - panel_u1s[1])/(2*hx);
    //     dx_u1s[7] = (3*(panel_u1s[7] - panel_u1s[4]) + (panel_u1s[1] - panel_u1s[4]))/(2*hx);

    //     dx_u1s[2] = (-3*(panel_u1s[2] - panel_u1s[5]) + (panel_u1s[5] - panel_u1s[8]))/(2*hx);
    //     dx_u1s[5] = (panel_u1s[8] - panel_u1s[2])/(2*hx);
    //     dx_u1s[8] = (3*(panel_u1s[8] - panel_u1s[5]) + (panel_u1s[2] - panel_u1s[5]))/(2*hx);

        
    //     dy_u1s[0] = (-3*(panel_u1s[0] - panel_u1s[1]) + (panel_u1s[1] - panel_u1s[2]))/(2*hy);
    //     dy_u1s[1] = (panel_u1s[2] - panel_u1s[0])/(2*hy);
    //     dy_u1s[2] = (3*(panel_u1s[2] - panel_u1s[1]) + (panel_u1s[0] - panel_u1s[1]))/(2*hy);

    //     dy_u1s[3] = (-3*(panel_u1s[3] - panel_u1s[4]) + (panel_u1s[4] - panel_u1s[5]))/(2*hy);
    //     dy_u1s[4] = (panel_u1s[5] - panel_u1s[3])/(2*hy);
    //     dy_u1s[5] = (3*(panel_u1s[5] - panel_u1s[4]) + (panel_u1s[3] - panel_u1s[4]))/(2*hy);

    //     dy_u1s[6] = (-3*(panel_u1s[6] - panel_u1s[7]) + (panel_u1s[7] - panel_u1s[8]))/(2*hy);
    //     dy_u1s[7] = (panel_u1s[8] - panel_u1s[6])/(2*hy);
    //     dy_u1s[8] = (3*(panel_u1s[8] - panel_u1s[7]) + (panel_u1s[6] - panel_u1s[7]))/(2*hy);
    
        
        
    //     //////////// u2s ///////////////
    //     dx_u2s[0] = (-3*(panel_u2s[0] - panel_u2s[3]) + (panel_u2s[3] - panel_u2s[6]))/(2*hx);
    //     dx_u2s[3] = (panel_u2s[6] - panel_u2s[0])/(2*hx);
    //     dx_u2s[6] = (3*(panel_u2s[6] - panel_u2s[3]) + (panel_u2s[0] - panel_u2s[3]))/(2*hx);

    //     dx_u2s[1] = (-3*(panel_u2s[1] - panel_u2s[4]) + (panel_u2s[4] - panel_u2s[7]))/(2*hx);
    //     dx_u2s[4] = (panel_u2s[7] - panel_u2s[1])/(2*hx);
    //     dx_u2s[7] = (3*(panel_u2s[7] - panel_u2s[4]) + (panel_u2s[1] - panel_u2s[4]))/(2*hx);

    //     dx_u2s[2] = (-3*(panel_u2s[2] - panel_u2s[5]) + (panel_u2s[5] - panel_u2s[8]))/(2*hx);
    //     dx_u2s[5] = (panel_u2s[8] - panel_u2s[2])/(2*hx);
    //     dx_u2s[8] = (3*(panel_u2s[8] - panel_u2s[5]) + (panel_u2s[2] - panel_u2s[5]))/(2*hx);

        
    //     dy_u2s[0] = (-3*(panel_u2s[0] - panel_u2s[1]) + (panel_u2s[1] - panel_u2s[2]))/(2*hy);
    //     dy_u2s[1] = (panel_u2s[2] - panel_u2s[0])/(2*hy);
    //     dy_u2s[2] = (3*(panel_u2s[2] - panel_u2s[1]) + (panel_u2s[0] - panel_u2s[1]))/(2*hy);

    //     dy_u2s[3] = (-3*(panel_u2s[3] - panel_u2s[4]) + (panel_u2s[4] - panel_u2s[5]))/(2*hy);
    //     dy_u2s[4] = (panel_u2s[5] - panel_u2s[3])/(2*hy);
    //     dy_u2s[5] = (3*(panel_u2s[5] - panel_u2s[4]) + (panel_u2s[3] - panel_u2s[4]))/(2*hy);

    //     dy_u2s[6] = (-3*(panel_u2s[6] - panel_u2s[7]) + (panel_u2s[7] - panel_u2s[8]))/(2*hy);
    //     dy_u2s[7] = (panel_u2s[8] - panel_u2s[6])/(2*hy);
    //     dy_u2s[8] = (3*(panel_u2s[8] - panel_u2s[7]) + (panel_u2s[6] - panel_u2s[7]))/(2*hy);
    
        
    //     //////////// b1s ///////////////
    //     dx_b1s[0] = (-3*(panel_b1s[0] - panel_b1s[3]) + (panel_b1s[3] - panel_b1s[6]))/(2*hx);
    //     dx_b1s[3] = (panel_b1s[6] - panel_b1s[0])/(2*hx);
    //     dx_b1s[6] = (3*(panel_b1s[6] - panel_b1s[3]) + (panel_b1s[0] - panel_b1s[3]))/(2*hx);

    //     dx_b1s[1] = (-3*(panel_b1s[1] - panel_b1s[4]) + (panel_b1s[4] - panel_b1s[7]))/(2*hx);
    //     dx_b1s[4] = (panel_b1s[7] - panel_b1s[1])/(2*hx);
    //     dx_b1s[7] = (3*(panel_b1s[7] - panel_b1s[4]) + (panel_b1s[1] - panel_b1s[4]))/(2*hx);

    //     dx_b1s[2] = (-3*(panel_b1s[2] - panel_b1s[5]) + (panel_b1s[5] - panel_b1s[8]))/(2*hx);
    //     dx_b1s[5] = (panel_b1s[8] - panel_b1s[2])/(2*hx);
    //     dx_b1s[8] = (3*(panel_b1s[8] - panel_b1s[5]) + (panel_b1s[2] - panel_b1s[5]))/(2*hx);

        
    //     dy_b1s[0] = (-3*(panel_b1s[0] - panel_b1s[1]) + (panel_b1s[1] - panel_b1s[2]))/(2*hy);
    //     dy_b1s[1] = (panel_b1s[2] - panel_b1s[0])/(2*hy);
    //     dy_b1s[2] = (3*(panel_b1s[2] - panel_b1s[1]) + (panel_b1s[0] - panel_b1s[1]))/(2*hy);

    //     dy_b1s[3] = (-3*(panel_b1s[3] - panel_b1s[4]) + (panel_b1s[4] - panel_b1s[5]))/(2*hy);
    //     dy_b1s[4] = (panel_b1s[5] - panel_b1s[3])/(2*hy);
    //     dy_b1s[5] = (3*(panel_b1s[5] - panel_b1s[4]) + (panel_b1s[3] - panel_b1s[4]))/(2*hy);

    //     dy_b1s[6] = (-3*(panel_b1s[6] - panel_b1s[7]) + (panel_b1s[7] - panel_b1s[8]))/(2*hy);
    //     dy_b1s[7] = (panel_b1s[8] - panel_b1s[6])/(2*hy);
    //     dy_b1s[8] = (3*(panel_b1s[8] - panel_b1s[7]) + (panel_b1s[6] - panel_b1s[7]))/(2*hy);
    



    //     //////////// b2s ///////////////
    //     dx_b2s[0] = (-3*(panel_b2s[0] - panel_b2s[3]) + (panel_b2s[3] - panel_b2s[6]))/(2*hx);
    //     dx_b2s[3] = (panel_b2s[6] - panel_b2s[0])/(2*hx);
    //     dx_b2s[6] = (3*(panel_b2s[6] - panel_b2s[3]) + (panel_b2s[0] - panel_b2s[3]))/(2*hx);

    //     dx_b2s[1] = (-3*(panel_b2s[1] - panel_b2s[4]) + (panel_b2s[4] - panel_b2s[7]))/(2*hx);
    //     dx_b2s[4] = (panel_b2s[7] - panel_b2s[1])/(2*hx);
    //     dx_b2s[7] = (3*(panel_b2s[7] - panel_b2s[4]) + (panel_b2s[1] - panel_b2s[4]))/(2*hx);

    //     dx_b2s[2] = (-3*(panel_b2s[2] - panel_b2s[5]) + (panel_b2s[5] - panel_b2s[8]))/(2*hx);
    //     dx_b2s[5] = (panel_b2s[8] - panel_b2s[2])/(2*hx);
    //     dx_b2s[8] = (3*(panel_b2s[8] - panel_b2s[5]) + (panel_b2s[2] - panel_b2s[5]))/(2*hx);

        
    //     dy_b2s[0] = (-3*(panel_b2s[0] - panel_b2s[1]) + (panel_b2s[1] - panel_b2s[2]))/(2*hy);
    //     dy_b2s[1] = (panel_b2s[2] - panel_b2s[0])/(2*hy);
    //     dy_b2s[2] = (3*(panel_b2s[2] - panel_b2s[1]) + (panel_b2s[0] - panel_b2s[1]))/(2*hy);

    //     dy_b2s[3] = (-3*(panel_b2s[3] - panel_b2s[4]) + (panel_b2s[4] - panel_b2s[5]))/(2*hy);
    //     dy_b2s[4] = (panel_b2s[5] - panel_b2s[3])/(2*hy);
    //     dy_b2s[5] = (3*(panel_b2s[5] - panel_b2s[4]) + (panel_b2s[3] - panel_b2s[4]))/(2*hy);

    //     dy_b2s[6] = (-3*(panel_b2s[6] - panel_b2s[7]) + (panel_b2s[7] - panel_b2s[8]))/(2*hy);
    //     dy_b2s[7] = (panel_b2s[8] - panel_b2s[6])/(2*hy);
    //     dy_b2s[8] = (3*(panel_b2s[8] - panel_b2s[7]) + (panel_b2s[6] - panel_b2s[7]))/(2*hy);
    


    //     for (int ii = 0; ii < 9; ++ii) {
    //         int pind = panel_point_inds[ii];
    //         vorticity_grad_x[pind] += dx_w0[ii];
    //         j_grad_x[pind]         += dx_j0[ii];
    //         u1s_grad_x[pind]       += dx_u1s[ii];
    //         u2s_grad_x[pind]       += dx_u2s[ii];
    //         b1s_grad_x[pind]       += dx_b1s[ii];
    //         b2s_grad_x[pind]       += dx_b2s[ii];

    //         vorticity_grad_y[pind] += dy_w0[ii];
    //         j_grad_y[pind]         += dy_j0[ii];
    //         u1s_grad_y[pind]       += dy_u1s[ii];
    //         u2s_grad_y[pind]       += dy_u2s[ii]; 
    //         b1s_grad_y[pind]       += dy_b1s[ii];
    //         b2s_grad_y[pind]       += dy_b2s[ii];

    //         grad_count[pind] += 1;
    //     }
    // }

    // // average gradients
    // for (int i = 0; i < xs.size(); ++i) {
    //     if (grad_count[i] > 0) {
    //         double inv = 1.0 / grad_count[i];
    //         vorticity_grad_x[i] *= inv; vorticity_grad_y[i] *= inv;
    //         j_grad_x[i] *= inv;         j_grad_y[i] *= inv;
    //         u1s_grad_x[i] *= inv;       u1s_grad_y[i] *= inv;
    //         u2s_grad_x[i] *= inv;       u2s_grad_y[i] *= inv;
    //         b1s_grad_x[i] *= inv;       b1s_grad_y[i] *= inv;
    //         b2s_grad_x[i] *= inv;       b2s_grad_y[i] *= inv;
    //     }
    // }

    // // treat boundary
    // for (int panel_ind = 0; panel_ind < panels.size(); panel_ind++) {
    //     Panel* panel = &(panels[panel_ind]);
    //     // only use leaf panels 
    //     if (panel->child_inds_start > -1) {
    //         continue;
    //     }
    //     if (panel->is_left_bdry) {
    //         // cout << " left bdry panel: " << panel->panel_ind << endl;
    //         const int* panel_point_inds = panel->point_inds;
    //         Panel* left_panel = &(panels[panel->left_nbr_ind]);
    //         const int* left_panel_point_inds = left_panel->point_inds;
    //         // three left boundary points of the panel, point index in the panel is 0,1,2
    //         // three right boundary points of the left panel, point index in the panel is 6,7,8
    //         int pind_0 = panel_point_inds[0];
    //         int pind_1 = panel_point_inds[1];
    //         int pind_2 = panel_point_inds[2];

    //         int left_pind_6 = left_panel_point_inds[6];
    //         int left_pind_7 = left_panel_point_inds[7];
    //         int left_pind_8 = left_panel_point_inds[8];

    //         vorticity_grad_x[pind_0] = (vorticity_grad_x[pind_0] + vorticity_grad_x[left_pind_6])/2;
    //         vorticity_grad_y[pind_0] = (vorticity_grad_y[pind_0] + vorticity_grad_y[left_pind_6])/2;
    //         j_grad_x[pind_0] = (j_grad_x[pind_0] + j_grad_x[left_pind_6])/2;
    //         j_grad_y[pind_0] = (j_grad_y[pind_0] + j_grad_y[left_pind_6])/2;
    //         u1s_grad_x[pind_0] = (u1s_grad_x[pind_0] + u1s_grad_x[left_pind_6])/2;
    //         u1s_grad_y[pind_0] = (u1s_grad_y[pind_0] + u1s_grad_y[left_pind_6])/2;
    //         u2s_grad_x[pind_0] = (u2s_grad_x[pind_0] + u2s_grad_x[left_pind_6])/2;
    //         u2s_grad_y[pind_0] = (u2s_grad_y[pind_0] + u2s_grad_y[left_pind_6])/2;
    //         b1s_grad_x[pind_0] = (b1s_grad_x[pind_0] + b1s_grad_x[left_pind_6])/2;
    //         b1s_grad_y[pind_0] = (b1s_grad_y[pind_0] + b1s_grad_y[left_pind_6])/2;
    //         b2s_grad_x[pind_0] = (b2s_grad_x[pind_0] + b2s_grad_x[left_pind_6])/2;
    //         b2s_grad_y[pind_0] = (b2s_grad_y[pind_0] + b2s_grad_y[left_pind_6])/2;

    //         vorticity_grad_x[pind_1] = (vorticity_grad_x[pind_1] + vorticity_grad_x[left_pind_7])/2;
    //         vorticity_grad_y[pind_1] = (vorticity_grad_y[pind_1] + vorticity_grad_y[left_pind_7])/2;
    //         j_grad_x[pind_1] = (j_grad_x[pind_1] + j_grad_x[left_pind_7])/2;
    //         j_grad_y[pind_1] = (j_grad_y[pind_1] + j_grad_y[left_pind_7])/2;
    //         u1s_grad_x[pind_1] = (u1s_grad_x[pind_1] + u1s_grad_x[left_pind_7])/2;
    //         u1s_grad_y[pind_1] = (u1s_grad_y[pind_1] + u1s_grad_y[left_pind_7])/2;
    //         u2s_grad_x[pind_1] = (u2s_grad_x[pind_1] + u2s_grad_x[left_pind_7])/2;
    //         u2s_grad_y[pind_1] = (u2s_grad_y[pind_1] + u2s_grad_y[left_pind_7])/2;
    //         b1s_grad_x[pind_1] = (b1s_grad_x[pind_1] + b1s_grad_x[left_pind_7])/2;
    //         b1s_grad_y[pind_1] = (b1s_grad_y[pind_1] + b1s_grad_y[left_pind_7])/2;
    //         b2s_grad_x[pind_1] = (b2s_grad_x[pind_1] + b2s_grad_x[left_pind_7])/2;
    //         b2s_grad_y[pind_1] = (b2s_grad_y[pind_1] + b2s_grad_y[left_pind_7])/2;

    //         vorticity_grad_x[pind_2] = (vorticity_grad_x[pind_2] + vorticity_grad_x[left_pind_8])/2;
    //         vorticity_grad_y[pind_2] = (vorticity_grad_y[pind_2] + vorticity_grad_y[left_pind_8])/2;
    //         j_grad_x[pind_2] = (j_grad_x[pind_2] + j_grad_x[left_pind_8])/2;
    //         j_grad_y[pind_2] = (j_grad_y[pind_2] + j_grad_y[left_pind_8])/2;
    //         u1s_grad_x[pind_2] = (u1s_grad_x[pind_2] + u1s_grad_x[left_pind_8])/2;
    //         u1s_grad_y[pind_2] = (u1s_grad_y[pind_2] + u1s_grad_y[left_pind_8])/2;
    //         u2s_grad_x[pind_2] = (u2s_grad_x[pind_2] + u2s_grad_x[left_pind_8])/2;
    //         u2s_grad_y[pind_2] = (u2s_grad_y[pind_2] + u2s_grad_y[left_pind_8])/2;
    //         b1s_grad_x[pind_2] = (b1s_grad_x[pind_2] + b1s_grad_x[left_pind_8])/2;
    //         b1s_grad_y[pind_2] = (b1s_grad_y[pind_2] + b1s_grad_y[left_pind_8])/2;
    //         b2s_grad_x[pind_2] = (b2s_grad_x[pind_2] + b2s_grad_x[left_pind_8])/2;
    //         b2s_grad_y[pind_2] = (b2s_grad_y[pind_2] + b2s_grad_y[left_pind_8])/2;


    //         // right boundary have the same value
    //         vorticity_grad_x[left_pind_6] = vorticity_grad_x[pind_0];
    //         vorticity_grad_y[left_pind_6] = vorticity_grad_y[pind_0];
    //         j_grad_x[left_pind_6] = j_grad_x[pind_0];
    //         j_grad_y[left_pind_6] = j_grad_y[pind_0];
    //         u1s_grad_x[left_pind_6] = u1s_grad_x[pind_0];
    //         u1s_grad_y[left_pind_6] = u1s_grad_y[pind_0];
    //         u2s_grad_x[left_pind_6] = u2s_grad_x[pind_0];
    //         u2s_grad_y[left_pind_6] = u2s_grad_y[pind_0];
    //         b1s_grad_x[left_pind_6] = b1s_grad_x[pind_0];
    //         b1s_grad_y[left_pind_6] = b1s_grad_y[pind_0];
    //         b2s_grad_x[left_pind_6] = b2s_grad_x[pind_0];
    //         b2s_grad_y[left_pind_6] = b2s_grad_y[pind_0];

    //         vorticity_grad_x[left_pind_7] = vorticity_grad_x[pind_1];
    //         vorticity_grad_y[left_pind_7] = vorticity_grad_y[pind_1];
    //         j_grad_x[left_pind_7] = j_grad_x[pind_1];
    //         j_grad_y[left_pind_7] = j_grad_y[pind_1];
    //         u1s_grad_x[left_pind_7] = u1s_grad_x[pind_1];
    //         u1s_grad_y[left_pind_7] = u1s_grad_y[pind_1];
    //         u2s_grad_x[left_pind_7] = u2s_grad_x[pind_1];
    //         u2s_grad_y[left_pind_7] = u2s_grad_y[pind_1];
    //         b1s_grad_x[left_pind_7] = b1s_grad_x[pind_1];
    //         b1s_grad_y[left_pind_7] = b1s_grad_y[pind_1];
    //         b2s_grad_x[left_pind_7] = b2s_grad_x[pind_1];
    //         b2s_grad_y[left_pind_7] = b2s_grad_y[pind_1];

    //         vorticity_grad_x[left_pind_8] = vorticity_grad_x[pind_2];
    //         vorticity_grad_y[left_pind_8] = vorticity_grad_y[pind_2];
    //         j_grad_x[left_pind_8] = j_grad_x[pind_2];
    //         j_grad_y[left_pind_8] = j_grad_y[pind_2];
    //         u1s_grad_x[left_pind_8] = u1s_grad_x[pind_2];
    //         u1s_grad_y[left_pind_8] = u1s_grad_y[pind_2];
    //         u2s_grad_x[left_pind_8] = u2s_grad_x[pind_2];
    //         u2s_grad_y[left_pind_8] = u2s_grad_y[pind_2];
    //         b1s_grad_x[left_pind_8] = b1s_grad_x[pind_2];
    //         b1s_grad_y[left_pind_8] = b1s_grad_y[pind_2];
    //         b2s_grad_x[left_pind_8] = b2s_grad_x[pind_2];
    //         b2s_grad_y[left_pind_8] = b2s_grad_y[pind_2];
    //     }
    // }

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


    // // -------------------------
    // // Stage 2 : k2_u1s, k2_u2s, k2_dws, k2_djs
    // // -------------------------
    // evaluate_u_field(k2_u1s, k2_u2s, xs2, ys2, u_weights, t + 0.5 * dt);
    // // B evaluation
    // b1s.assign(b1s.size(), 0.0);
    // b2s.assign(b2s.size(), 0.0);
    // evaluate_b_field(b1s, b2s, xs2, ys2, b_weights, t+ 0.5 * dt);





    // // -------------------------
    // // Stage 3: k3_u1s, k3_u2s, k3_dws, k3_djs
    // // -------------------------
    // evaluate_u_field(k3_u1s, k3_u2s, xs3, ys3, u_weights, t + 0.5 * dt);
    // b1s.assign(b1s.size(), 0.0);
    // b2s.assign(b2s.size(), 0.0);
    // evaluate_b_field(b1s, b2s, xs3, ys3, b_weights, t+ 0.5 * dt);

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



    // // -------------------------
    // // Stage 4: k4_u1s, k4_u2s, k4_dws, k4_djs
    // // -------------------------
    // // evaluate_u_field(k4x, k4y, x4, y4, u_weights, t + dt);

    // evaluate_u_field(k4_u1s, k4_u2s, xs4, ys4, u_weights, t + dt);

    // // B evaluation
    // b1s.assign(b1s.size(), 0.0);
    // b2s.assign(b2s.size(), 0.0);
    // evaluate_b_field(b1s, b2s, xs4, ys4, b_weights, t + dt);



    return 0;
}
