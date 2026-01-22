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
    // rk4();

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

    #ifdef DEBUG
        cout << "u1s/u2s first 5:" << endl;
        for (int i = 0; i < std::min<int>(5, (int)u1s.size()); ++i) {
            cout << i << " u1=" << u1s[i] << " u2=" << u2s[i] << endl;
        }
    #endif

    // B evaluation
    b1s.assign(b1s.size(), 0.0);
    b2s.assign(b2s.size(), 0.0);
    evaluate_b_field(b1s, b2s, xs, ys, b_weights, t);

    // u1s_grad, u2s_grad, b1s_grad, b2s_grad evaluation

    u1s_grad_x.assign(xs.size(), 0.0);
    u1s_grad_y.assign(xs.size(), 0.0);

    evaluate_u1s_grad(u1s_grad_x, u1s_grad_y, xs, ys, u_weights, t);

    u2s_grad_x.assign(xs.size(), 0.0);
    u2s_grad_y.assign(xs.size(), 0.0);

    evaluate_u2s_grad(u2s_grad_x, u2s_grad_y, xs, ys, u_weights, t);

    b1s_grad_x.assign(xs.size(), 0.0);
    b1s_grad_y.assign(xs.size(), 0.0);

    evaluate_b1s_grad(b1s_grad_x, b1s_grad_y, xs, ys, b_weights, t);

    b2s_grad_x.assign(xs.size(), 0.0);
    b2s_grad_y.assign(xs.size(), 0.0);

    evaluate_b2s_grad(b2s_grad_x, b2s_grad_y, xs, ys, b_weights, t);

    // vortex_grad, j_grad evaluation
    vorticity_grad_x.assign(xs.size(), 0.0);
    vorticity_grad_y.assign(xs.size(), 0.0);

    evaluate_vorticity_grad(vorticity_grad_x, vorticity_grad_y, xs, ys, u_weights, t);

    j_grad_x.assign(xs.size(), 0.0);
    j_grad_y.assign(xs.size(), 0.0);

    evaluate_j_grad(j_grad_x, j_grad_y, xs, ys, b_weights, t);

    // vortex_laplacian, j_laplacian evaluation

    vorticity_laplacian.assign(xs.size(), 0.0);
    j_laplacian.assign(xs.size(), 0.0);

    std::vector<double> vorticity_none_local(xs.size(), 0.0);
    std::vector<double> j_none_local(xs.size(), 0.0);

    evaluate_vorticity_laplacian(vorticity_laplacian, vorticity_none_local, xs, ys, u_weights, t);
    evaluate_j_laplacian(j_laplacian, j_none_local, xs, ys, b_weights, t);


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
    }

    for (int i = 0; i < ys.size(); i++) {
        ys[i] += dt * u2s[i];
    }

    for (int i = 0; i < xs.size(); i++) {
        w0s[i] += dt * (nu * vorticity_laplacian[i] + B_dot_grad_j[i]);
    }

    for (int i = 0; i < xs.size(); i++) {
        j0s[i] += dt * (mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i]);
    }

    return 0;
}




int AMRStructure::rk4() {
    std::cout << "enter rk4" << std::endl;

    const int xs_size = (int)xs.size();

    // Stage storage: k1, k2, k3, k4 for x and y
    std::vector<double> k1x(xs_size, 0.0), k1y(xs_size, 0.0);
    std::vector<double> k2x(xs_size, 0.0), k2y(xs_size, 0.0);
    std::vector<double> k3x(xs_size, 0.0), k3y(xs_size, 0.0);
    std::vector<double> k4x(xs_size, 0.0), k4y(xs_size, 0.0);

    // Temporary positions
    std::vector<double> x2(xs_size), y2(xs_size);
    std::vector<double> x3(xs_size), y3(xs_size);
    std::vector<double> x4(xs_size), y4(xs_size);

    // -------------------------
    // k1 = u(x_n, t_n)
    // -------------------------

    evaluate_u_field(k1x, k1y, xs, ys, u_weights, t);
    for (int i = 0; i < xs_size; ++i) {
        // k1x[i] = u1s[i];
        // k1y[i] = u2s[i];
        x2[i]  = xs[i] + 0.5 * dt * k1x[i];
        y2[i]  = ys[i] + 0.5 * dt * k1y[i];
    }

    // -------------------------
    // k2 = u(x_n + dt/2*k1, t_n + dt/2)
    // -------------------------

    evaluate_u_field(k2x, k2y, x2, y2, u_weights, t + 0.5 * dt);
    for (int i = 0; i < xs_size; ++i) {
        x3[i]  = xs[i] + 0.5 * dt * k2x[i];
        y3[i]  = ys[i] + 0.5 * dt * k2y[i];
    }

    // -------------------------
    // k3 = u(x_n + dt/2*k2, t_n + dt/2)
    // -------------------------
    evaluate_u_field(k3x, k3y, x3, y3, u_weights, t + 0.5 * dt);
    for (int i = 0; i < xs_size; ++i) {
        x4[i]  = xs[i] + dt * k3x[i];
        y4[i]  = ys[i] + dt * k3y[i];
    }

    // -------------------------
    // k4 = u(x_n + dt*k3, t_n + dt)
    // -------------------------
    evaluate_u_field(k4x, k4y, x4, y4, u_weights, t + dt);

    // -------------------------
    // Combine to update positions
    // x_{n+1} = x_n + dt/6*(k1 + 2k2 + 2k3 + k4)
    // -------------------------
    for (int i = 0; i < xs_size; ++i) {
        xs[i] += (dt / 6.0) * (k1x[i] + 2.0 * k2x[i] + 2.0 * k3x[i] + k4x[i]);
        ys[i] += (dt / 6.0) * (k1y[i] + 2.0 * k2y[i] + 2.0 * k3y[i] + k4y[i]);
    }


    // Debug print
    std::cout << "after rk4 u field evaluation (combined) - first 5:" << std::endl;
    for (int i = 0; i < std::min(5, xs_size); ++i) {
        // std::cout << i
        //           << " k1=(" << k1x[i] << "," << k1y[i] << ")"
        //           << " k2=(" << k2x[i] << "," << k2y[i] << ")"
        //           << " k3=(" << k3x[i] << "," << k3y[i] << ")"
        //           << " k4=(" << k4x[i] << "," << k4y[i] << ")"
        //           << std::endl;
        std::cout << i
                  << " xs=(" << xs[i] << ")"
                  << " ys=(" << ys[i] << ")"
                  << std::endl;
    }

    return 0;
}
