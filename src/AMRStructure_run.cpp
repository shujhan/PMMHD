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

    // rk4 step
    // rk4_step(false);

    // euler 
    euler();
    t += dt;

    // // if remesh:
    // if (iter_num % n_steps_remesh == 0) {
    //     remesh();
    // }
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
    // std::vector<double> u1s_local (xs.size());
    // std::vector<double> u2s_local (xs.size());
    // std::vector<double> b1s_local (xs.size());
    // std::vector<double> b2s_local (xs.size());

    evaluate_u_field(u1s, u2s, xs, ys, u_weights, t);
    evaluate_b_field(b1s, b2s, xs, ys, b_weights, t);
    for (int i = 0; i < xs.size(); i++) {
        xs[i] += dt * u1s[i];
    }

    for (int i = 0; i < ys.size(); i++) {
        ys[i] += dt * u2s[i];
    }

    return 0;
}