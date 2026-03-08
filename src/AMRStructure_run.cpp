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

        // // start barycentric lagrange polynomial
        // double wx_0 = 1/((panel_xs[0]- panel_xs[3]) * (panel_xs[0]- panel_xs[6]));
        // double wx_1 = 1/((panel_xs[3]- panel_xs[0]) * (panel_xs[3]- panel_xs[6]));
        // double wx_2 = 1/((panel_xs[6]- panel_xs[0]) * (panel_xs[6]- panel_xs[3]));

        // double wy_0 = 1/((panel_ys[0]- panel_ys[1]) * (panel_ys[0]- panel_ys[2]));
        // double wy_1 = 1/((panel_ys[1]- panel_ys[0]) * (panel_ys[1]- panel_ys[2]));
        // double wy_2 = 1/((panel_ys[2]- panel_ys[0]) * (panel_ys[2]- panel_ys[1]));


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
        // first row
        // dx_j0[0] = wx_1/wx_0 * (panel_j0s[3] - panel_j0s[0])/(panel_xs[0] - panel_xs[3]) + wx_2/wx_0 * (panel_j0s[6] - panel_j0s[0])/(panel_xs[0] - panel_xs[6]);
        // dx_j0[3] = wx_0/wx_1 * (panel_j0s[0] - panel_j0s[3])/(panel_xs[3] - panel_xs[0]) + wx_2/wx_1 * (panel_j0s[6] - panel_j0s[3])/(panel_xs[3] - panel_xs[6]);
        // dx_j0[6] = wx_0/wx_2 * (panel_j0s[0] - panel_j0s[6])/(panel_xs[6] - panel_xs[0]) + wx_1/wx_2 * (panel_j0s[3] - panel_j0s[6])/(panel_xs[6] - panel_xs[3]);
        
        // // second row
        // dx_j0[1] = wx_1/wx_0 * (panel_j0s[4] - panel_j0s[1])/(panel_xs[1] - panel_xs[4]) + wx_2/wx_0 * (panel_j0s[7] - panel_j0s[1])/(panel_xs[1] - panel_xs[7]);
        // dx_j0[4] = wx_0/wx_1 * (panel_j0s[1] - panel_j0s[4])/(panel_xs[4] - panel_xs[1]) + wx_2/wx_1 * (panel_j0s[7] - panel_j0s[4])/(panel_xs[4] - panel_xs[7]);
        // dx_j0[7] = wx_0/wx_2 * (panel_j0s[1] - panel_j0s[7])/(panel_xs[7] - panel_xs[1]) + wx_1/wx_2 * (panel_j0s[4] - panel_j0s[7])/(panel_xs[7] - panel_xs[4]);
        
        // // third row
        // dx_j0[2] = wx_1/wx_0 * (panel_j0s[5] - panel_j0s[2])/(panel_xs[2] - panel_xs[5]) + wx_2/wx_0 * (panel_j0s[8] - panel_j0s[2])/(panel_xs[2] - panel_xs[8]);
        // dx_j0[5] = wx_0/wx_1 * (panel_j0s[2] - panel_j0s[5])/(panel_xs[5] - panel_xs[2]) + wx_2/wx_1 * (panel_j0s[8] - panel_j0s[5])/(panel_xs[5] - panel_xs[8]);
        // dx_j0[8] = wx_0/wx_2 * (panel_j0s[2] - panel_j0s[8])/(panel_xs[8] - panel_xs[2]) + wx_1/wx_2 * (panel_j0s[5] - panel_j0s[8])/(panel_xs[8] - panel_xs[5]);
        
        double hx = panel_xs[3]- panel_xs[0];
        // dx_j0[0] = (-3*panel_j0s[0] + 4* panel_j0s[3] - panel_j0s[6])/(2*hx);
        dx_j0[0] = (-3*(panel_j0s[0] - panel_j0s[3]) + (panel_j0s[3] - panel_j0s[6]))/(2*hx);
        dx_j0[3] = (panel_j0s[6] - panel_j0s[0])/(2*hx);
        // dx_j0[6] = (3*panel_j0s[6] - 4*panel_j0s[3] + panel_j0s[0])/(2*hx);
        dx_j0[6] = (3*(panel_j0s[6] - panel_j0s[3]) + (panel_j0s[0] - panel_j0s[3]))/(2*hx);

        // dx_j0[1] = (-3*panel_j0s[1] + 4* panel_j0s[4] - panel_j0s[7])/(2*hx);
        dx_j0[1] = (-3*(panel_j0s[1] - panel_j0s[4]) + (panel_j0s[4] - panel_j0s[7]))/(2*hx);
        dx_j0[4] = (panel_j0s[7] - panel_j0s[1])/(2*hx);
        // dx_j0[7] = (3*panel_j0s[7] - 4*panel_j0s[4] + panel_j0s[1])/(2*hx);
        dx_j0[7] = (3*(panel_j0s[7] - panel_j0s[4]) + (panel_j0s[1] - panel_j0s[4]))/(2*hx);

        dx_j0[2] = (-3*(panel_j0s[2] - panel_j0s[5]) + (panel_j0s[5] - panel_j0s[8]))/(2*hx);
        dx_j0[5] = (panel_j0s[8] - panel_j0s[2])/(2*hx);
        dx_j0[8] = (3*(panel_j0s[8] - panel_j0s[5]) + (panel_j0s[2] - panel_j0s[5]))/(2*hx);

        // cout << panel_j0s[0] << ", " << panel_j0s[3] << ", " << panel_j0s[6] << endl;
        // cout << panel_j0s[1] << ", " << panel_j0s[4] << ", " << panel_j0s[7] << endl;
        // cout << panel_j0s[2] << ", " << panel_j0s[5] << ", " << panel_j0s[8] << endl;
        // cout << (-3*(panel_j0s[0] - panel_j0s[3]) + (panel_j0s[3] - panel_j0s[6])) << endl;
        // cout << (panel_j0s[6] - panel_j0s[0]) << endl;
        // cout << (3*(panel_j0s[6] - panel_j0s[3]) + (panel_j0s[0] - panel_j0s[3])) << endl;

        // cout << (-3*(panel_j0s[1] - panel_j0s[4]) + (panel_j0s[4] - panel_j0s[7])) << endl;
        // cout << (panel_j0s[7] - panel_j0s[1]) << endl;
        // cout << (3*(panel_j0s[7] - panel_j0s[4]) + (panel_j0s[1] - panel_j0s[4])) << endl;

        // cout << (-3*(panel_j0s[2] - panel_j0s[5]) + (panel_j0s[5] - panel_j0s[8])) << endl;
        // cout << (panel_j0s[8] - panel_j0s[2]) << endl;
        // cout << (3*(panel_j0s[8] - panel_j0s[5]) + (panel_j0s[2] - panel_j0s[5])) << endl;

        // first column
        // dy_j0[0] = wy_1/wy_0 * (panel_j0s[1] - panel_j0s[0])/(panel_ys[0] - panel_ys[1]) + wy_2/wy_0 * (panel_j0s[2] - panel_j0s[0])/(panel_ys[0] - panel_ys[2]);
        // dy_j0[1] = wy_0/wy_1 * (panel_j0s[0] - panel_j0s[1])/(panel_ys[1] - panel_ys[0]) + wy_2/wy_1 * (panel_j0s[2] - panel_j0s[1])/(panel_ys[1] - panel_ys[2]);
        // dy_j0[2] = wy_0/wy_2 * (panel_j0s[0] - panel_j0s[2])/(panel_ys[2] - panel_ys[0]) + wy_1/wy_2 * (panel_j0s[1] - panel_j0s[2])/(panel_ys[2] - panel_ys[1]);
        
        // // second column
        // dy_j0[3] = wy_1/wy_0 * (panel_j0s[4] - panel_j0s[3])/(panel_ys[3] - panel_ys[4]) + wy_2/wy_0 * (panel_j0s[5] - panel_j0s[3])/(panel_ys[3] - panel_ys[5]);
        // dy_j0[4] = wy_0/wy_1 * (panel_j0s[3] - panel_j0s[4])/(panel_ys[4] - panel_ys[3]) + wy_2/wy_1 * (panel_j0s[5] - panel_j0s[4])/(panel_ys[4] - panel_ys[5]);
        // dy_j0[5] = wy_0/wy_2 * (panel_j0s[3] - panel_j0s[5])/(panel_ys[5] - panel_ys[3]) + wy_1/wy_2 * (panel_j0s[4] - panel_j0s[5])/(panel_ys[5] - panel_ys[4]);
        
        // // third column
        // dy_j0[6] = wy_1/wy_0 * (panel_j0s[7] - panel_j0s[6])/(panel_ys[6] - panel_ys[7]) + wy_2/wy_0 * (panel_j0s[8] - panel_j0s[6])/(panel_ys[6] - panel_ys[8]);
        // dy_j0[7] = wy_0/wy_1 * (panel_j0s[6] - panel_j0s[7])/(panel_ys[7] - panel_ys[6]) + wy_2/wy_1 * (panel_j0s[8] - panel_j0s[7])/(panel_ys[7] - panel_ys[8]);
        // dy_j0[8] = wy_0/wy_2 * (panel_j0s[6] - panel_j0s[8])/(panel_ys[8] - panel_ys[6]) + wy_1/wy_2 * (panel_j0s[7] - panel_j0s[8])/(panel_ys[8] - panel_ys[7]);
        

        double hy = panel_ys[1]- panel_ys[0];
        // dy_j0[0] = (-3*panel_j0s[0] + 4* panel_j0s[1] - panel_j0s[2])/(2*hy);
        dy_j0[0] = (-3*(panel_j0s[0] - panel_j0s[1]) + (panel_j0s[1] - panel_j0s[2]))/(2*hy);
        dy_j0[1] = (panel_j0s[2] - panel_j0s[0])/(2*hy);
        // dy_j0[2] = (3*panel_j0s[2] - 4*panel_j0s[1] + panel_j0s[0])/(2*hy);
        dy_j0[2] = (3*(panel_j0s[2] - panel_j0s[1]) + (panel_j0s[0] - panel_j0s[1]))/(2*hy);

        dy_j0[3] = (-3*(panel_j0s[3] - panel_j0s[4]) + (panel_j0s[4] - panel_j0s[5]))/(2*hy);
        dy_j0[4] = (panel_j0s[5] - panel_j0s[3])/(2*hy);
        dy_j0[5] = (3*(panel_j0s[5] - panel_j0s[4]) + (panel_j0s[3] - panel_j0s[4]))/(2*hy);

        dy_j0[6] = (-3*(panel_j0s[6] - panel_j0s[7]) + (panel_j0s[7] - panel_j0s[8]))/(2*hy);
        dy_j0[7] = (panel_j0s[8] - panel_j0s[6])/(2*hy);
        dy_j0[8] = (3*(panel_j0s[8] - panel_j0s[7]) + (panel_j0s[6] - panel_j0s[7]))/(2*hy);



        //////////// w0 ///////////////
        // dx_w0[0] = wx_1/wx_0 * (panel_w0s[3] - panel_w0s[0])/(panel_xs[0] - panel_xs[3]) + wx_2/wx_0 * (panel_w0s[6] - panel_w0s[0])/(panel_xs[0] - panel_xs[6]);
        // dx_w0[3] = wx_0/wx_1 * (panel_w0s[0] - panel_w0s[3])/(panel_xs[3] - panel_xs[0]) + wx_2/wx_1 * (panel_w0s[6] - panel_w0s[3])/(panel_xs[3] - panel_xs[6]);
        // dx_w0[6] = wx_0/wx_2 * (panel_w0s[0] - panel_w0s[6])/(panel_xs[6] - panel_xs[0]) + wx_1/wx_2 * (panel_w0s[3] - panel_w0s[6])/(panel_xs[6] - panel_xs[3]);
        // dx_w0[1] = wx_1/wx_0 * (panel_w0s[4] - panel_w0s[1])/(panel_xs[1] - panel_xs[4]) + wx_2/wx_0 * (panel_w0s[7] - panel_w0s[1])/(panel_xs[1] - panel_xs[7]);
        // dx_w0[4] = wx_0/wx_1 * (panel_w0s[1] - panel_w0s[4])/(panel_xs[4] - panel_xs[1]) + wx_2/wx_1 * (panel_w0s[7] - panel_w0s[4])/(panel_xs[4] - panel_xs[7]);
        // dx_w0[7] = wx_0/wx_2 * (panel_w0s[1] - panel_w0s[7])/(panel_xs[7] - panel_xs[1]) + wx_1/wx_2 * (panel_w0s[4] - panel_w0s[7])/(panel_xs[7] - panel_xs[4]);
        // dx_w0[2] = wx_1/wx_0 * (panel_w0s[5] - panel_w0s[2])/(panel_xs[2] - panel_xs[5]) + wx_2/wx_0 * (panel_w0s[8] - panel_w0s[2])/(panel_xs[2] - panel_xs[8]);
        // dx_w0[5] = wx_0/wx_1 * (panel_w0s[2] - panel_w0s[5])/(panel_xs[5] - panel_xs[2]) + wx_2/wx_1 * (panel_w0s[8] - panel_w0s[5])/(panel_xs[5] - panel_xs[8]);
        // dx_w0[8] = wx_0/wx_2 * (panel_w0s[2] - panel_w0s[8])/(panel_xs[8] - panel_xs[2]) + wx_1/wx_2 * (panel_w0s[5] - panel_w0s[8])/(panel_xs[8] - panel_xs[5]);

        // dy_w0[0] = wy_1/wy_0 * (panel_w0s[1] - panel_w0s[0])/(panel_ys[0] - panel_ys[1]) + wy_2/wy_0 * (panel_w0s[2] - panel_w0s[0])/(panel_ys[0] - panel_ys[2]);
        // dy_w0[1] = wy_0/wy_1 * (panel_w0s[0] - panel_w0s[1])/(panel_ys[1] - panel_ys[0]) + wy_2/wy_1 * (panel_w0s[2] - panel_w0s[1])/(panel_ys[1] - panel_ys[2]);
        // dy_w0[2] = wy_0/wy_2 * (panel_w0s[0] - panel_w0s[2])/(panel_ys[2] - panel_ys[0]) + wy_1/wy_2 * (panel_w0s[1] - panel_w0s[2])/(panel_ys[2] - panel_ys[1]);
        // dy_w0[3] = wy_1/wy_0 * (panel_w0s[4] - panel_w0s[3])/(panel_ys[3] - panel_ys[4]) + wy_2/wy_0 * (panel_w0s[5] - panel_w0s[3])/(panel_ys[3] - panel_ys[5]);
        // dy_w0[4] = wy_0/wy_1 * (panel_w0s[3] - panel_w0s[4])/(panel_ys[4] - panel_ys[3]) + wy_2/wy_1 * (panel_w0s[5] - panel_w0s[4])/(panel_ys[4] - panel_ys[5]);
        // dy_w0[5] = wy_0/wy_2 * (panel_w0s[3] - panel_w0s[5])/(panel_ys[5] - panel_ys[3]) + wy_1/wy_2 * (panel_w0s[4] - panel_w0s[5])/(panel_ys[5] - panel_ys[4]);
        // dy_w0[6] = wy_1/wy_0 * (panel_w0s[7] - panel_w0s[6])/(panel_ys[6] - panel_ys[7]) + wy_2/wy_0 * (panel_w0s[8] - panel_w0s[6])/(panel_ys[6] - panel_ys[8]);
        // dy_w0[7] = wy_0/wy_1 * (panel_w0s[6] - panel_w0s[7])/(panel_ys[7] - panel_ys[6]) + wy_2/wy_1 * (panel_w0s[8] - panel_w0s[7])/(panel_ys[7] - panel_ys[8]);
        // dy_w0[8] = wy_0/wy_2 * (panel_w0s[6] - panel_w0s[8])/(panel_ys[8] - panel_ys[6]) + wy_1/wy_2 * (panel_w0s[7] - panel_w0s[8])/(panel_ys[8] - panel_ys[7]);
        
        dx_w0[0] = (-3*panel_w0s[0] + 4* panel_w0s[3] - panel_w0s[6])/(2*hx);
        dx_w0[3] = (panel_w0s[6] - panel_w0s[0])/(2*hx);
        dx_w0[6] = (3*panel_w0s[6] - 4*panel_w0s[3] + panel_w0s[0])/(2*hx);

        dx_w0[1] = (-3*panel_w0s[1] + 4* panel_w0s[4] - panel_w0s[7])/(2*hx);
        dx_w0[4] = (panel_w0s[7] - panel_w0s[1])/(2*hx);
        dx_w0[7] = (3*panel_w0s[7] - 4*panel_w0s[4] + panel_w0s[1])/(2*hx);

        dx_w0[2] = (-3*panel_w0s[2] + 4* panel_w0s[5] - panel_w0s[8])/(2*hx);
        dx_w0[5] = (panel_w0s[8] - panel_w0s[2])/(2*hx);
        dx_w0[8] = (3*panel_w0s[8] - 4*panel_w0s[5] + panel_w0s[2])/(2*hx);
        
        
        dy_w0[0] = (-3*panel_w0s[0] + 4* panel_w0s[1] - panel_w0s[2])/(2*hy);
        dy_w0[1] = (panel_w0s[2] - panel_w0s[0])/(2*hy);
        dy_w0[2] = (3*panel_w0s[2] - 4*panel_w0s[1] + panel_w0s[0])/(2*hy);

        dy_w0[3] = (-3*panel_w0s[3] + 4* panel_w0s[4] - panel_w0s[5])/(2*hy);
        dy_w0[4] = (panel_w0s[5] - panel_w0s[3])/(2*hy);
        dy_w0[5] = (3*panel_w0s[5] - 4*panel_w0s[4] + panel_w0s[3])/(2*hy);

        dy_w0[6] = (-3*panel_w0s[6] + 4* panel_w0s[7] - panel_w0s[8])/(2*hy);
        dy_w0[7] = (panel_w0s[8] - panel_w0s[6])/(2*hy);
        dy_w0[8] = (3*panel_w0s[8] - 4*panel_w0s[7] + panel_w0s[6])/(2*hy);        
        
        
        
        
        
        //////////// u1s ///////////////
        // dx_u1s[0] = wx_1/wx_0 * (panel_u1s[3] - panel_u1s[0])/(panel_xs[0] - panel_xs[3]) + wx_2/wx_0 * (panel_u1s[6] - panel_u1s[0])/(panel_xs[0] - panel_xs[6]);
        // dx_u1s[3] = wx_0/wx_1 * (panel_u1s[0] - panel_u1s[3])/(panel_xs[3] - panel_xs[0]) + wx_2/wx_1 * (panel_u1s[6] - panel_u1s[3])/(panel_xs[3] - panel_xs[6]);
        // dx_u1s[6] = wx_0/wx_2 * (panel_u1s[0] - panel_u1s[6])/(panel_xs[6] - panel_xs[0]) + wx_1/wx_2 * (panel_u1s[3] - panel_u1s[6])/(panel_xs[6] - panel_xs[3]);
        // dx_u1s[1] = wx_1/wx_0 * (panel_u1s[4] - panel_u1s[1])/(panel_xs[1] - panel_xs[4]) + wx_2/wx_0 * (panel_u1s[7] - panel_u1s[1])/(panel_xs[1] - panel_xs[7]);
        // dx_u1s[4] = wx_0/wx_1 * (panel_u1s[1] - panel_u1s[4])/(panel_xs[4] - panel_xs[1]) + wx_2/wx_1 * (panel_u1s[7] - panel_u1s[4])/(panel_xs[4] - panel_xs[7]);
        // dx_u1s[7] = wx_0/wx_2 * (panel_u1s[1] - panel_u1s[7])/(panel_xs[7] - panel_xs[1]) + wx_1/wx_2 * (panel_u1s[4] - panel_u1s[7])/(panel_xs[7] - panel_xs[4]);
        // dx_u1s[2] = wx_1/wx_0 * (panel_u1s[5] - panel_u1s[2])/(panel_xs[2] - panel_xs[5]) + wx_2/wx_0 * (panel_u1s[8] - panel_u1s[2])/(panel_xs[2] - panel_xs[8]);
        // dx_u1s[5] = wx_0/wx_1 * (panel_u1s[2] - panel_u1s[5])/(panel_xs[5] - panel_xs[2]) + wx_2/wx_1 * (panel_u1s[8] - panel_u1s[5])/(panel_xs[5] - panel_xs[8]);
        // dx_u1s[8] = wx_0/wx_2 * (panel_u1s[2] - panel_u1s[8])/(panel_xs[8] - panel_xs[2]) + wx_1/wx_2 * (panel_u1s[5] - panel_u1s[8])/(panel_xs[8] - panel_xs[5]);

        // dy_u1s[0] = wy_1/wy_0 * (panel_u1s[1] - panel_u1s[0])/(panel_ys[0] - panel_ys[1]) + wy_2/wy_0 * (panel_u1s[2] - panel_u1s[0])/(panel_ys[0] - panel_ys[2]);
        // dy_u1s[1] = wy_0/wy_1 * (panel_u1s[0] - panel_u1s[1])/(panel_ys[1] - panel_ys[0]) + wy_2/wy_1 * (panel_u1s[2] - panel_u1s[1])/(panel_ys[1] - panel_ys[2]);
        // dy_u1s[2] = wy_0/wy_2 * (panel_u1s[0] - panel_u1s[2])/(panel_ys[2] - panel_ys[0]) + wy_1/wy_2 * (panel_u1s[1] - panel_u1s[2])/(panel_ys[2] - panel_ys[1]);
        // dy_u1s[3] = wy_1/wy_0 * (panel_u1s[4] - panel_u1s[3])/(panel_ys[3] - panel_ys[4]) + wy_2/wy_0 * (panel_u1s[5] - panel_u1s[3])/(panel_ys[3] - panel_ys[5]);
        // dy_u1s[4] = wy_0/wy_1 * (panel_u1s[3] - panel_u1s[4])/(panel_ys[4] - panel_ys[3]) + wy_2/wy_1 * (panel_u1s[5] - panel_u1s[4])/(panel_ys[4] - panel_ys[5]);
        // dy_u1s[5] = wy_0/wy_2 * (panel_u1s[3] - panel_u1s[5])/(panel_ys[5] - panel_ys[3]) + wy_1/wy_2 * (panel_u1s[4] - panel_u1s[5])/(panel_ys[5] - panel_ys[4]);
        // dy_u1s[6] = wy_1/wy_0 * (panel_u1s[7] - panel_u1s[6])/(panel_ys[6] - panel_ys[7]) + wy_2/wy_0 * (panel_u1s[8] - panel_u1s[6])/(panel_ys[6] - panel_ys[8]);
        // dy_u1s[7] = wy_0/wy_1 * (panel_u1s[6] - panel_u1s[7])/(panel_ys[7] - panel_ys[6]) + wy_2/wy_1 * (panel_u1s[8] - panel_u1s[7])/(panel_ys[7] - panel_ys[8]);
        // dy_u1s[8] = wy_0/wy_2 * (panel_u1s[6] - panel_u1s[8])/(panel_ys[8] - panel_ys[6]) + wy_1/wy_2 * (panel_u1s[7] - panel_u1s[8])/(panel_ys[8] - panel_ys[7]);
        
        
        dx_u1s[0] = (-3*panel_u1s[0] + 4* panel_u1s[3] - panel_u1s[6])/(2*hx);
        dx_u1s[3] = (panel_u1s[6] - panel_u1s[0])/(2*hx);
        dx_u1s[6] = (3*panel_u1s[6] - 4*panel_u1s[3] + panel_u1s[0])/(2*hx);

        dx_u1s[1] = (-3*panel_u1s[1] + 4* panel_u1s[4] - panel_u1s[7])/(2*hx);
        dx_u1s[4] = (panel_u1s[7] - panel_u1s[1])/(2*hx);
        dx_u1s[7] = (3*panel_u1s[7] - 4*panel_u1s[4] + panel_u1s[1])/(2*hx);

        dx_u1s[2] = (-3*panel_u1s[2] + 4* panel_u1s[5] - panel_u1s[8])/(2*hx);
        dx_u1s[5] = (panel_u1s[8] - panel_u1s[2])/(2*hx);
        dx_u1s[8] = (3*panel_u1s[8] - 4*panel_u1s[5] + panel_u1s[2])/(2*hx);
        
        
        dy_u1s[0] = (-3*panel_u1s[0] + 4* panel_u1s[1] - panel_u1s[2])/(2*hy);
        dy_u1s[1] = (panel_u1s[2] - panel_u1s[0])/(2*hy);
        dy_u1s[2] = (3*panel_u1s[2] - 4*panel_u1s[1] + panel_u1s[0])/(2*hy);

        dy_u1s[3] = (-3*panel_u1s[3] + 4* panel_u1s[4] - panel_u1s[5])/(2*hy);
        dy_u1s[4] = (panel_u1s[5] - panel_u1s[3])/(2*hy);
        dy_u1s[5] = (3*panel_u1s[5] - 4*panel_u1s[4] + panel_u1s[3])/(2*hy);

        dy_u1s[6] = (-3*panel_u1s[6] + 4* panel_u1s[7] - panel_u1s[8])/(2*hy);
        dy_u1s[7] = (panel_u1s[8] - panel_u1s[6])/(2*hy);
        dy_u1s[8] = (3*panel_u1s[8] - 4*panel_u1s[7] + panel_u1s[6])/(2*hy);        
        
        
        
        //////////// u2s ///////////////
        // dx_u2s[0] = wx_1/wx_0 * (panel_u2s[3] - panel_u2s[0])/(panel_xs[0] - panel_xs[3]) + wx_2/wx_0 * (panel_u2s[6] - panel_u2s[0])/(panel_xs[0] - panel_xs[6]);
        // dx_u2s[3] = wx_0/wx_1 * (panel_u2s[0] - panel_u2s[3])/(panel_xs[3] - panel_xs[0]) + wx_2/wx_1 * (panel_u2s[6] - panel_u2s[3])/(panel_xs[3] - panel_xs[6]);
        // dx_u2s[6] = wx_0/wx_2 * (panel_u2s[0] - panel_u2s[6])/(panel_xs[6] - panel_xs[0]) + wx_1/wx_2 * (panel_u2s[3] - panel_u2s[6])/(panel_xs[6] - panel_xs[3]);
        // dx_u2s[1] = wx_1/wx_0 * (panel_u2s[4] - panel_u2s[1])/(panel_xs[1] - panel_xs[4]) + wx_2/wx_0 * (panel_u2s[7] - panel_u2s[1])/(panel_xs[1] - panel_xs[7]);
        // dx_u2s[4] = wx_0/wx_1 * (panel_u2s[1] - panel_u2s[4])/(panel_xs[4] - panel_xs[1]) + wx_2/wx_1 * (panel_u2s[7] - panel_u2s[4])/(panel_xs[4] - panel_xs[7]);
        // dx_u2s[7] = wx_0/wx_2 * (panel_u2s[1] - panel_u2s[7])/(panel_xs[7] - panel_xs[1]) + wx_1/wx_2 * (panel_u2s[4] - panel_u2s[7])/(panel_xs[7] - panel_xs[4]);
        // dx_u2s[2] = wx_1/wx_0 * (panel_u2s[5] - panel_u2s[2])/(panel_xs[2] - panel_xs[5]) + wx_2/wx_0 * (panel_u2s[8] - panel_u2s[2])/(panel_xs[2] - panel_xs[8]);
        // dx_u2s[5] = wx_0/wx_1 * (panel_u2s[2] - panel_u2s[5])/(panel_xs[5] - panel_xs[2]) + wx_2/wx_1 * (panel_u2s[8] - panel_u2s[5])/(panel_xs[5] - panel_xs[8]);
        // dx_u2s[8] = wx_0/wx_2 * (panel_u2s[2] - panel_u2s[8])/(panel_xs[8] - panel_xs[2]) + wx_1/wx_2 * (panel_u2s[5] - panel_u2s[8])/(panel_xs[8] - panel_xs[5]);

        // dy_u2s[0] = wy_1/wy_0 * (panel_u2s[1] - panel_u2s[0])/(panel_ys[0] - panel_ys[1]) + wy_2/wy_0 * (panel_u2s[2] - panel_u2s[0])/(panel_ys[0] - panel_ys[2]);
        // dy_u2s[1] = wy_0/wy_1 * (panel_u2s[0] - panel_u2s[1])/(panel_ys[1] - panel_ys[0]) + wy_2/wy_1 * (panel_u2s[2] - panel_u2s[1])/(panel_ys[1] - panel_ys[2]);
        // dy_u2s[2] = wy_0/wy_2 * (panel_u2s[0] - panel_u2s[2])/(panel_ys[2] - panel_ys[0]) + wy_1/wy_2 * (panel_u2s[1] - panel_u2s[2])/(panel_ys[2] - panel_ys[1]);
        // dy_u2s[3] = wy_1/wy_0 * (panel_u2s[4] - panel_u2s[3])/(panel_ys[3] - panel_ys[4]) + wy_2/wy_0 * (panel_u2s[5] - panel_u2s[3])/(panel_ys[3] - panel_ys[5]);
        // dy_u2s[4] = wy_0/wy_1 * (panel_u2s[3] - panel_u2s[4])/(panel_ys[4] - panel_ys[3]) + wy_2/wy_1 * (panel_u2s[5] - panel_u2s[4])/(panel_ys[4] - panel_ys[5]);
        // dy_u2s[5] = wy_0/wy_2 * (panel_u2s[3] - panel_u2s[5])/(panel_ys[5] - panel_ys[3]) + wy_1/wy_2 * (panel_u2s[4] - panel_u2s[5])/(panel_ys[5] - panel_ys[4]);
        // dy_u2s[6] = wy_1/wy_0 * (panel_u2s[7] - panel_u2s[6])/(panel_ys[6] - panel_ys[7]) + wy_2/wy_0 * (panel_u2s[8] - panel_u2s[6])/(panel_ys[6] - panel_ys[8]);
        // dy_u2s[7] = wy_0/wy_1 * (panel_u2s[6] - panel_u2s[7])/(panel_ys[7] - panel_ys[6]) + wy_2/wy_1 * (panel_u2s[8] - panel_u2s[7])/(panel_ys[7] - panel_ys[8]);
        // dy_u2s[8] = wy_0/wy_2 * (panel_u2s[6] - panel_u2s[8])/(panel_ys[8] - panel_ys[6]) + wy_1/wy_2 * (panel_u2s[7] - panel_u2s[8])/(panel_ys[8] - panel_ys[7]);
        
        
        dx_u2s[0] = (-3*panel_u2s[0] + 4* panel_u2s[3] - panel_u2s[6])/(2*hx);
        dx_u2s[3] = (panel_u2s[6] - panel_u2s[0])/(2*hx);
        dx_u2s[6] = (3*panel_u2s[6] - 4*panel_u2s[3] + panel_u2s[0])/(2*hx);

        dx_u2s[1] = (-3*panel_u2s[1] + 4* panel_u2s[4] - panel_u2s[7])/(2*hx);
        dx_u2s[4] = (panel_u2s[7] - panel_u2s[1])/(2*hx);
        dx_u2s[7] = (3*panel_u2s[7] - 4*panel_u2s[4] + panel_u2s[1])/(2*hx);

        dx_u2s[2] = (-3*panel_u2s[2] + 4* panel_u2s[5] - panel_u2s[8])/(2*hx);
        dx_u2s[5] = (panel_u2s[8] - panel_u2s[2])/(2*hx);
        dx_u2s[8] = (3*panel_u2s[8] - 4*panel_u2s[5] + panel_u2s[2])/(2*hx);
        
        
        dy_u2s[0] = (-3*panel_u2s[0] + 4* panel_u2s[1] - panel_u2s[2])/(2*hy);
        dy_u2s[1] = (panel_u2s[2] - panel_u2s[0])/(2*hy);
        dy_u2s[2] = (3*panel_u2s[2] - 4*panel_u2s[1] + panel_u2s[0])/(2*hy);

        dy_u2s[3] = (-3*panel_u2s[3] + 4* panel_u2s[4] - panel_u2s[5])/(2*hy);
        dy_u2s[4] = (panel_u2s[5] - panel_u2s[3])/(2*hy);
        dy_u2s[5] = (3*panel_u2s[5] - 4*panel_u2s[4] + panel_u2s[3])/(2*hy);

        dy_u2s[6] = (-3*panel_u2s[6] + 4* panel_u2s[7] - panel_u2s[8])/(2*hy);
        dy_u2s[7] = (panel_u2s[8] - panel_u2s[6])/(2*hy);
        dy_u2s[8] = (3*panel_u2s[8] - 4*panel_u2s[7] + panel_u2s[6])/(2*hy);   
        
        
        //////////// b1s ///////////////
        // dx_b1s[0] = wx_1/wx_0 * (panel_b1s[3] - panel_b1s[0])/(panel_xs[0] - panel_xs[3]) + wx_2/wx_0 * (panel_b1s[6] - panel_b1s[0])/(panel_xs[0] - panel_xs[6]);
        // dx_b1s[3] = wx_0/wx_1 * (panel_b1s[0] - panel_b1s[3])/(panel_xs[3] - panel_xs[0]) + wx_2/wx_1 * (panel_b1s[6] - panel_b1s[3])/(panel_xs[3] - panel_xs[6]);
        // dx_b1s[6] = wx_0/wx_2 * (panel_b1s[0] - panel_b1s[6])/(panel_xs[6] - panel_xs[0]) + wx_1/wx_2 * (panel_b1s[3] - panel_b1s[6])/(panel_xs[6] - panel_xs[3]);
        // dx_b1s[1] = wx_1/wx_0 * (panel_b1s[4] - panel_b1s[1])/(panel_xs[1] - panel_xs[4]) + wx_2/wx_0 * (panel_b1s[7] - panel_b1s[1])/(panel_xs[1] - panel_xs[7]);
        // dx_b1s[4] = wx_0/wx_1 * (panel_b1s[1] - panel_b1s[4])/(panel_xs[4] - panel_xs[1]) + wx_2/wx_1 * (panel_b1s[7] - panel_b1s[4])/(panel_xs[4] - panel_xs[7]);
        // dx_b1s[7] = wx_0/wx_2 * (panel_b1s[1] - panel_b1s[7])/(panel_xs[7] - panel_xs[1]) + wx_1/wx_2 * (panel_b1s[4] - panel_b1s[7])/(panel_xs[7] - panel_xs[4]);
        // dx_b1s[2] = wx_1/wx_0 * (panel_b1s[5] - panel_b1s[2])/(panel_xs[2] - panel_xs[5]) + wx_2/wx_0 * (panel_b1s[8] - panel_b1s[2])/(panel_xs[2] - panel_xs[8]);
        // dx_b1s[5] = wx_0/wx_1 * (panel_b1s[2] - panel_b1s[5])/(panel_xs[5] - panel_xs[2]) + wx_2/wx_1 * (panel_b1s[8] - panel_b1s[5])/(panel_xs[5] - panel_xs[8]);
        // dx_b1s[8] = wx_0/wx_2 * (panel_b1s[2] - panel_b1s[8])/(panel_xs[8] - panel_xs[2]) + wx_1/wx_2 * (panel_b1s[5] - panel_b1s[8])/(panel_xs[8] - panel_xs[5]);

        // dy_b1s[0] = wy_1/wy_0 * (panel_b1s[1] - panel_b1s[0])/(panel_ys[0] - panel_ys[1]) + wy_2/wy_0 * (panel_b1s[2] - panel_b1s[0])/(panel_ys[0] - panel_ys[2]);
        // dy_b1s[1] = wy_0/wy_1 * (panel_b1s[0] - panel_b1s[1])/(panel_ys[1] - panel_ys[0]) + wy_2/wy_1 * (panel_b1s[2] - panel_b1s[1])/(panel_ys[1] - panel_ys[2]);
        // dy_b1s[2] = wy_0/wy_2 * (panel_b1s[0] - panel_b1s[2])/(panel_ys[2] - panel_ys[0]) + wy_1/wy_2 * (panel_b1s[1] - panel_b1s[2])/(panel_ys[2] - panel_ys[1]);
        // dy_b1s[3] = wy_1/wy_0 * (panel_b1s[4] - panel_b1s[3])/(panel_ys[3] - panel_ys[4]) + wy_2/wy_0 * (panel_b1s[5] - panel_b1s[3])/(panel_ys[3] - panel_ys[5]);
        // dy_b1s[4] = wy_0/wy_1 * (panel_b1s[3] - panel_b1s[4])/(panel_ys[4] - panel_ys[3]) + wy_2/wy_1 * (panel_b1s[5] - panel_b1s[4])/(panel_ys[4] - panel_ys[5]);
        // dy_b1s[5] = wy_0/wy_2 * (panel_b1s[3] - panel_b1s[5])/(panel_ys[5] - panel_ys[3]) + wy_1/wy_2 * (panel_b1s[4] - panel_b1s[5])/(panel_ys[5] - panel_ys[4]);
        // dy_b1s[6] = wy_1/wy_0 * (panel_b1s[7] - panel_b1s[6])/(panel_ys[6] - panel_ys[7]) + wy_2/wy_0 * (panel_b1s[8] - panel_b1s[6])/(panel_ys[6] - panel_ys[8]);
        // dy_b1s[7] = wy_0/wy_1 * (panel_b1s[6] - panel_b1s[7])/(panel_ys[7] - panel_ys[6]) + wy_2/wy_1 * (panel_b1s[8] - panel_b1s[7])/(panel_ys[7] - panel_ys[8]);
        // dy_b1s[8] = wy_0/wy_2 * (panel_b1s[6] - panel_b1s[8])/(panel_ys[8] - panel_ys[6]) + wy_1/wy_2 * (panel_b1s[7] - panel_b1s[8])/(panel_ys[8] - panel_ys[7]);
        

        dx_b1s[0] = (-3*panel_b1s[0] + 4* panel_b1s[3] - panel_b1s[6])/(2*hx);
        dx_b1s[3] = (panel_b1s[6] - panel_b1s[0])/(2*hx);
        dx_b1s[6] = (3*panel_b1s[6] - 4*panel_b1s[3] + panel_b1s[0])/(2*hx);

        dx_b1s[1] = (-3*panel_b1s[1] + 4* panel_b1s[4] - panel_b1s[7])/(2*hx);
        dx_b1s[4] = (panel_b1s[7] - panel_b1s[1])/(2*hx);
        dx_b1s[7] = (3*panel_b1s[7] - 4*panel_b1s[4] + panel_b1s[1])/(2*hx);

        dx_b1s[2] = (-3*panel_b1s[2] + 4* panel_b1s[5] - panel_b1s[8])/(2*hx);
        dx_b1s[5] = (panel_b1s[8] - panel_b1s[2])/(2*hx);
        dx_b1s[8] = (3*panel_b1s[8] - 4*panel_b1s[5] + panel_b1s[2])/(2*hx);
        
        
        dy_b1s[0] = (-3*panel_b1s[0] + 4* panel_b1s[1] - panel_b1s[2])/(2*hy);
        dy_b1s[1] = (panel_b1s[2] - panel_b1s[0])/(2*hy);
        dy_b1s[2] = (3*panel_b1s[2] - 4*panel_b1s[1] + panel_b1s[0])/(2*hy);

        dy_b1s[3] = (-3*panel_b1s[3] + 4* panel_b1s[4] - panel_b1s[5])/(2*hy);
        dy_b1s[4] = (panel_b1s[5] - panel_b1s[3])/(2*hy);
        dy_b1s[5] = (3*panel_b1s[5] - 4*panel_b1s[4] + panel_b1s[3])/(2*hy);

        dy_b1s[6] = (-3*panel_b1s[6] + 4* panel_b1s[7] - panel_b1s[8])/(2*hy);
        dy_b1s[7] = (panel_b1s[8] - panel_b1s[6])/(2*hy);
        dy_b1s[8] = (3*panel_b1s[8] - 4*panel_b1s[7] + panel_b1s[6])/(2*hy);   




        //////////// b2s ///////////////
        // dx_b2s[0] = wx_1/wx_0 * (panel_b2s[3] - panel_b2s[0])/(panel_xs[0] - panel_xs[3]) + wx_2/wx_0 * (panel_b2s[6] - panel_b2s[0])/(panel_xs[0] - panel_xs[6]);
        // dx_b2s[3] = wx_0/wx_1 * (panel_b2s[0] - panel_b2s[3])/(panel_xs[3] - panel_xs[0]) + wx_2/wx_1 * (panel_b2s[6] - panel_b2s[3])/(panel_xs[3] - panel_xs[6]);
        // dx_b2s[6] = wx_0/wx_2 * (panel_b2s[0] - panel_b2s[6])/(panel_xs[6] - panel_xs[0]) + wx_1/wx_2 * (panel_b2s[3] - panel_b2s[6])/(panel_xs[6] - panel_xs[3]);
        // dx_b2s[1] = wx_1/wx_0 * (panel_b2s[4] - panel_b2s[1])/(panel_xs[1] - panel_xs[4]) + wx_2/wx_0 * (panel_b2s[7] - panel_b2s[1])/(panel_xs[1] - panel_xs[7]);
        // dx_b2s[4] = wx_0/wx_1 * (panel_b2s[1] - panel_b2s[4])/(panel_xs[4] - panel_xs[1]) + wx_2/wx_1 * (panel_b2s[7] - panel_b2s[4])/(panel_xs[4] - panel_xs[7]);
        // dx_b2s[7] = wx_0/wx_2 * (panel_b2s[1] - panel_b2s[7])/(panel_xs[7] - panel_xs[1]) + wx_1/wx_2 * (panel_b2s[4] - panel_b2s[7])/(panel_xs[7] - panel_xs[4]);
        // dx_b2s[2] = wx_1/wx_0 * (panel_b2s[5] - panel_b2s[2])/(panel_xs[2] - panel_xs[5]) + wx_2/wx_0 * (panel_b2s[8] - panel_b2s[2])/(panel_xs[2] - panel_xs[8]);
        // dx_b2s[5] = wx_0/wx_1 * (panel_b2s[2] - panel_b2s[5])/(panel_xs[5] - panel_xs[2]) + wx_2/wx_1 * (panel_b2s[8] - panel_b2s[5])/(panel_xs[5] - panel_xs[8]);
        // dx_b2s[8] = wx_0/wx_2 * (panel_b2s[2] - panel_b2s[8])/(panel_xs[8] - panel_xs[2]) + wx_1/wx_2 * (panel_b2s[5] - panel_b2s[8])/(panel_xs[8] - panel_xs[5]);

        // dy_b2s[0] = wy_1/wy_0 * (panel_b2s[1] - panel_b2s[0])/(panel_ys[0] - panel_ys[1]) + wy_2/wy_0 * (panel_b2s[2] - panel_b2s[0])/(panel_ys[0] - panel_ys[2]);
        // dy_b2s[1] = wy_0/wy_1 * (panel_b2s[0] - panel_b2s[1])/(panel_ys[1] - panel_ys[0]) + wy_2/wy_1 * (panel_b2s[2] - panel_b2s[1])/(panel_ys[1] - panel_ys[2]);
        // dy_b2s[2] = wy_0/wy_2 * (panel_b2s[0] - panel_b2s[2])/(panel_ys[2] - panel_ys[0]) + wy_1/wy_2 * (panel_b2s[1] - panel_b2s[2])/(panel_ys[2] - panel_ys[1]);
        // dy_b2s[3] = wy_1/wy_0 * (panel_b2s[4] - panel_b2s[3])/(panel_ys[3] - panel_ys[4]) + wy_2/wy_0 * (panel_b2s[5] - panel_b2s[3])/(panel_ys[3] - panel_ys[5]);
        // dy_b2s[4] = wy_0/wy_1 * (panel_b2s[3] - panel_b2s[4])/(panel_ys[4] - panel_ys[3]) + wy_2/wy_1 * (panel_b2s[5] - panel_b2s[4])/(panel_ys[4] - panel_ys[5]);
        // dy_b2s[5] = wy_0/wy_2 * (panel_b2s[3] - panel_b2s[5])/(panel_ys[5] - panel_ys[3]) + wy_1/wy_2 * (panel_b2s[4] - panel_b2s[5])/(panel_ys[5] - panel_ys[4]);
        // dy_b2s[6] = wy_1/wy_0 * (panel_b2s[7] - panel_b2s[6])/(panel_ys[6] - panel_ys[7]) + wy_2/wy_0 * (panel_b2s[8] - panel_b2s[6])/(panel_ys[6] - panel_ys[8]);
        // dy_b2s[7] = wy_0/wy_1 * (panel_b2s[6] - panel_b2s[7])/(panel_ys[7] - panel_ys[6]) + wy_2/wy_1 * (panel_b2s[8] - panel_b2s[7])/(panel_ys[7] - panel_ys[8]);
        // dy_b2s[8] = wy_0/wy_2 * (panel_b2s[6] - panel_b2s[8])/(panel_ys[8] - panel_ys[6]) + wy_1/wy_2 * (panel_b2s[7] - panel_b2s[8])/(panel_ys[8] - panel_ys[7]);
        

        dx_b2s[0] = (-3*panel_b2s[0] + 4* panel_b2s[3] - panel_b2s[6])/(2*hx);
        dx_b2s[3] = (panel_b2s[6] - panel_b2s[0])/(2*hx);
        dx_b2s[6] = (3*panel_b2s[6] - 4*panel_b2s[3] + panel_b2s[0])/(2*hx);

        dx_b2s[1] = (-3*panel_b2s[1] + 4* panel_b2s[4] - panel_b2s[7])/(2*hx);
        dx_b2s[4] = (panel_b2s[7] - panel_b2s[1])/(2*hx);
        dx_b2s[7] = (3*panel_b2s[7] - 4*panel_b2s[4] + panel_b2s[1])/(2*hx);

        dx_b2s[2] = (-3*panel_b2s[2] + 4* panel_b2s[5] - panel_b2s[8])/(2*hx);
        dx_b2s[5] = (panel_b2s[8] - panel_b2s[2])/(2*hx);
        dx_b2s[8] = (3*panel_b2s[8] - 4*panel_b2s[5] + panel_b2s[2])/(2*hx);
        
        
        dy_b2s[0] = (-3*panel_b2s[0] + 4* panel_b2s[1] - panel_b2s[2])/(2*hy);
        dy_b2s[1] = (panel_b2s[2] - panel_b2s[0])/(2*hy);
        dy_b2s[2] = (3*panel_b2s[2] - 4*panel_b2s[1] + panel_b2s[0])/(2*hy);

        dy_b2s[3] = (-3*panel_b2s[3] + 4* panel_b2s[4] - panel_b2s[5])/(2*hy);
        dy_b2s[4] = (panel_b2s[5] - panel_b2s[3])/(2*hy);
        dy_b2s[5] = (3*panel_b2s[5] - 4*panel_b2s[4] + panel_b2s[3])/(2*hy);

        dy_b2s[6] = (-3*panel_b2s[6] + 4* panel_b2s[7] - panel_b2s[8])/(2*hy);
        dy_b2s[7] = (panel_b2s[8] - panel_b2s[6])/(2*hy);
        dy_b2s[8] = (3*panel_b2s[8] - 4*panel_b2s[7] + panel_b2s[6])/(2*hy);   




        for (int ii = 0; ii < 9; ++ii) {
            int pind = panel_point_inds[ii];
            vorticity_grad_x[pind] += dx_w0[ii];
            j_grad_x[pind]         += dx_j0[ii];
            u1s_grad_x[pind]       += dx_u1s[ii];
            u2s_grad_x[pind]       += dx_u2s[ii];
            b1s_grad_x[pind]       += dx_b1s[ii];
            b2s_grad_x[pind]       += dx_b2s[ii];

            vorticity_grad_y[pind] += dy_w0[ii];
            j_grad_y[pind]         += dy_j0[ii];
            u1s_grad_y[pind]       += dy_u1s[ii];
            u2s_grad_y[pind]       += dy_u2s[ii]; 
            b1s_grad_y[pind]       += dy_b1s[ii];
            b2s_grad_y[pind]       += dy_b2s[ii];

            grad_count[pind] += 1;
        }


        // print out the gradient for each panel 
        // cout << "leaf panel: " << panel->panel_ind <<endl;
        // for (int ii = 0; ii < 9; ++ii) {
        //     int pind = panel_point_inds[ii];
        //     cout << "i=" << ii
        //     << " x=" << xs[pind]
        //     << " y=" << ys[pind]
        //     << " b1=" << b1s[pind]
        //     << " b2=" << b2s[pind]
        //     << " j=" << j0s[pind]
        //     // << " w=" << w0s[pind]
        //     // << " w_dx=" << dx_w0[ii]
        //     // << " w_dy=" << dy_w0[ii]

        //     // << " u1_dx=" << dx_u1s[ii]
        //     // << " u1_dy=" << dy_u1s[ii]
        //     // << " u2_dx=" << dx_u2s[ii]
        //     // << " u2_dy=" << dy_u2s[ii]

        //     // << " b1_dx=" << dx_u1s[ii]
        //     // << " b1_dy=" << dy_u1s[ii]
        //     // << " b2_dx=" << dx_u2s[ii]
        //     // << " b2_dy=" << dy_u2s[ii]
        //     << " j_dx=" << dx_j0[ii]
        //     << " j_dy=" << dy_j0[ii]
        //     << "\n";
        // }

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
            // cout << " left bdry panel: " << panel->panel_ind << endl;
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




    // // Lax–Friedrichs, loop through panels 
    // std::vector<int> average_count(xs.size(), 0);

    // std::vector<int> xs_average_state(xs.size(), 0);
    // std::vector<int> w0s_average_state(w0s.size(), 0);

    // for (int panel_ind = 0; panel_ind < panels.size(); panel_ind++) {
    //     Panel* panel = &(panels[panel_ind]);
    //     if (panel->child_inds_start > -1) {
    //         continue;
    //     }
    //     const int* panel_point_inds = panel->point_inds;
    //     double panel_xs[9], panel_ys[9];
    //     double panel_w0s[9], panel_j0s[9];

    //     for (int ii = 0; ii < 9; ++ii) {
    //         int pind = panel_point_inds[ii];
    //         panel_xs[ii] = xs[pind];
    //         panel_ys[ii] = ys[pind];
    //         panel_w0s[ii] = w0s[pind];
    //         panel_j0s[ii] = j0s[pind];
    //     }

    //     // go through each index 0 to 8
    //     w0s_average_state[panel_point_inds[0]] += panel_w0s[1] + panel_w0s[3];
    //     average_count[panel_point_inds[0]] += 2;
    //     w0s_average_state[panel_point_inds[1]] += panel_w0s[0] + panel_w0s[2] + panel_w0s[4];
    //     average_count[panel_point_inds[1]] += 3;
    //     w0s_average_state[panel_point_inds[2]] += panel_w0s[1] + panel_w0s[5];
    //     average_count[panel_point_inds[2]] += 2;
    //     w0s_average_state[panel_point_inds[3]] += panel_w0s[0] + panel_w0s[4] + panel_w0s[6];
    //     average_count[panel_point_inds[3]] += 3;
    //     w0s_average_state[panel_point_inds[4]] += panel_w0s[1] + panel_w0s[3] + panel_w0s[5] + panel_w0s[7];
    //     average_count[panel_point_inds[4]] += 4;
    //     w0s_average_state[panel_point_inds[5]] += panel_w0s[2] + panel_w0s[4] + panel_w0s[8];
    //     average_count[panel_point_inds[5]] += 3;
    //     w0s_average_state[panel_point_inds[6]] += panel_w0s[3] + panel_w0s[7];
    //     average_count[panel_point_inds[6]] += 2;
    //     w0s_average_state[panel_point_inds[7]] += panel_w0s[4] + panel_w0s[6] + panel_w0s[8];
    //     average_count[panel_point_inds[7]] += 3;
    //     w0s_average_state[panel_point_inds[8]] += panel_w0s[5] + panel_w0s[7];
    //     average_count[panel_point_inds[8]] += 2;
    //     if (panel->is_left_bdry) {

    //     }

    //     if (panel->is_top_bdry) {
            
    //     }

    //     if (panel->is_bottom_bdry) {
            
    //     }

    // }

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
