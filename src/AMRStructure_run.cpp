#include "AMRStructure.hpp"
 
int AMRStructure::run() {
    while (iter_num < num_steps) {
        // only for kelvin helmholz test, set vorticity for each point at initial step:
        if (iter_num == 0){
            // find points of vortex sheet at middle
            int vortex_size = power(2, initial_height + 1) + 1;
            std::vector<double> vortex_sheet_x(vortex_size, 0.0);
            std::vector<double> vortex_sheet_y(vortex_size, 0.0);
            double middle = (y_min + y_max) / 2;
            int count = 0;
            for (int i = 0; i < ys.size(); i++) {
                if(ys[i] == middle) {
                    vortex_sheet_x[count] = xs[i];
                    vortex_sheet_y[count] = ys[i];
                    count++;
                }
                else {
                    // not at vortex sheet, set to be 0 
                    w0s[i] = 0.0;
                }
            }
            cout << "vortex_sheet size = " << count-- << endl;

            // set all other points' vorticity 
            double reg_delta = 0.5;
            const double pi = M_PI;
            for (int i = 0; i < xs.size(); i++) {
                for (int k = 0; k < vortex_sheet_x.size(); k++) {
                    if(xs[i] != vortex_sheet_x[k] && ys[i] != vortex_sheet_y[k]) {
                        double y_diff = ys[i] - vortex_sheet_y[k];
                        double x_diff = xs[i] - vortex_sheet_x[k];
                        double cst_period = 2*pi/Lx;
                        w0s[i] += pi * reg_delta * reg_delta / Lx /  vortex_sheet_x.size() * 
                                (cosh(cst_period * y_diff) + cos(cst_period * x_diff)) 
                                / (cosh(cst_period * y_diff) + cos(cst_period * x_diff) + reg_delta * reg_delta) / (cosh(cst_period * y_diff) + cos(cst_period * x_diff) + reg_delta * reg_delta);           
                    }
                }
            }
            // give a perturbation to the vortex sheet 


            // write to file
            write_to_file();
        }
       


        step();
    }
    return 0;
}

int AMRStructure::step() {
    iter_num += 1;
    std::cout << "step " << iter_num << std::endl;

    // rk4 step
    // rk4_step(false);

    // euler 
    euler();
    t += dt;

    // if remesh:
    if (iter_num % n_steps_remesh == 0) {
        remesh();
    }
    // // if not remesh: evaluate e (remeshing ends by calculating e on uniform grid)
    // else {
    //     evaluate_field(es, xs, q_ws, t);
    // }

    // if dump : write to file
    if (iter_num % n_steps_diag == 0) {
        write_to_file();
    }

    return 0;
}


int AMRStructure::euler() {
    cout << "enter euler" << endl;
    u1s.assign(u1s.size(), 0.0);
    u2s.assign(u2s.size(), 0.0);
    evaluate_u_field(u1s, u2s, xs, ys, u_weights, t);
    cout << "after u field evaluation" << endl;
    cout << "u1s/u2s first 5:" << endl;
    for (int i = 0; i < std::min<int>(5, (int)u1s.size()); ++i) {
        cout << i << " u1=" << u1s[i] << " u2=" << u2s[i] << endl;
    }
    b1s.assign(b1s.size(), 0.0);
    b2s.assign(b2s.size(), 0.0);
    evaluate_b_field(b1s, b2s, xs, ys, b_weights, t);
    cout << "after b field evaluation" << endl;
    for (int i = 0; i < xs.size(); i++) {
        xs[i] += dt * u1s[i];
    }

    for (int i = 0; i < ys.size(); i++) {
        ys[i] += dt * u2s[i];
    }

    return 0;
}