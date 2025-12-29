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