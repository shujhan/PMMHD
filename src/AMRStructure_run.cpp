#include "AMRStructure.hpp"
 
int AMRStructure::run() {
    while (iter_num < num_steps) {
        step();
    }
    return 0;
}

int AMRSimulation::step() {
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


int AMRSimulation::euler() {
    evaluate_u_field();
    evaluate_b_field();
    for (int i = 0; i < xs.size(); i++) {
        xs[i] += dt * u1_s[i];
    }

    for (int i = 0; i < ys.size(); i++) {
        ys[i] += dt * u2_s[i];
    }

    return 0;
}