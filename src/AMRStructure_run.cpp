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

    // rebuild u_weights amd b_weights : first divide by w/j before pushing
    // then multiply new w/j after pushing 
    // or just use weights[ii] * ws2, js2
    
    // // divide w/j before pushing 
    // for (int i = 0; i < u_weights.size(); i++) {
    //     u_weights[i] = u_weights[i] / w0s[i];
    //     b_weights[i] = b_weights[i] / j0s[i];
    // }

    // pushing 
    for (int i = 0; i < xs.size(); i++) {
        xs2[i] = xs[i] + 0.5 * dt * u1s[i];
    }
    for (int i = 0; i < ys.size(); i++) {
        ys2[i] = ys[i] + 0.5 * dt * u2s[i];
    }
    for (int i = 0; i < xs.size(); i++) {
        k1_dws[i] = nu * vorticity_laplacian[i] + B_dot_grad_j[i]; 
        ws2[i] = w0s[i] + 0.5 * dt * k1_dws[i];
    }
    for (int i = 0; i < xs.size(); i++) {
        k1_djs[i] = mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i];
        js2[i] = j0s[i] + 0.5 * dt * k1_djs[i];
    }

    // multiply new w/j after pushing 
    for (int i = 0; i < u_weights.size(); i++) {
        u_weights[i] = weights[i] * ws2[i];
        b_weights[i] = weights[i] * js2[i];
    }


    // -------------------------
    // Stage 2 : k2_u1s, k2_u2s, k2_dws, k2_djs
    // -------------------------

    evaluate_u_field(k2_u1s, k2_u2s, xs2, ys2, u_weights, t + 0.5 * dt);

    // B evaluation
    b1s.assign(b1s.size(), 0.0);
    b2s.assign(b2s.size(), 0.0);
    evaluate_b_field(b1s, b2s, xs2, ys2, b_weights, t+ 0.5 * dt);

    // u1s_grad, u2s_grad, b1s_grad, b2s_grad evaluation
    u1s_grad_x.assign(xs.size(), 0.0);
    u1s_grad_y.assign(xs.size(), 0.0);
    evaluate_u1s_grad(u1s_grad_x, u1s_grad_y, xs2, ys2, u_weights, t+ 0.5 * dt);


    u2s_grad_x.assign(xs.size(), 0.0);
    u2s_grad_y.assign(xs.size(), 0.0);

    evaluate_u2s_grad(u2s_grad_x, u2s_grad_y, xs2, ys2, u_weights, t+ 0.5 * dt);

    b1s_grad_x.assign(xs.size(), 0.0);
    b1s_grad_y.assign(xs.size(), 0.0);

    evaluate_b1s_grad(b1s_grad_x, b1s_grad_y, xs2, ys2, b_weights, t+ 0.5 * dt);

    b2s_grad_x.assign(xs.size(), 0.0);
    b2s_grad_y.assign(xs.size(), 0.0);

    evaluate_b2s_grad(b2s_grad_x, b2s_grad_y, xs2, ys2, b_weights, t+ 0.5 * dt);

    // vortex_grad, j_grad evaluation
    vorticity_grad_x.assign(xs.size(), 0.0);
    vorticity_grad_y.assign(xs.size(), 0.0);

    evaluate_vorticity_grad(vorticity_grad_x, vorticity_grad_y, xs2, ys2, u_weights, t+ 0.5 * dt);

    j_grad_x.assign(xs.size(), 0.0);
    j_grad_y.assign(xs.size(), 0.0);

    evaluate_j_grad(j_grad_x, j_grad_y, xs2, ys2, b_weights, t+ 0.5 * dt);

    // vortex_laplacian, j_laplacian evaluation

    vorticity_laplacian.assign(xs.size(), 0.0);
    j_laplacian.assign(xs.size(), 0.0);

    vorticity_none_local.assign(xs.size(), 0.0);
    j_none_local.assign(xs.size(), 0.0);

    evaluate_vorticity_laplacian(vorticity_laplacian, vorticity_none_local, xs2, ys2, u_weights, t+ 0.5 * dt);
    evaluate_j_laplacian(j_laplacian, j_none_local, xs2, ys2, b_weights, t+ 0.5 * dt);

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


    // pushing 
    for (int i = 0; i < xs.size(); i++) {
        xs3[i] = xs[i] + 0.5 * dt * k2_u1s[i];
    }
    for (int i = 0; i < ys.size(); i++) {
        ys3[i] = ys[i] + 0.5 * dt * k2_u2s[i];
    }
    for (int i = 0; i < xs.size(); i++) {
        k2_dws[i] = nu * vorticity_laplacian[i] + B_dot_grad_j[i]; 
        ws3[i] = w0s[i] + 0.5 * dt * k2_dws[i];
    }
    for (int i = 0; i < xs.size(); i++) {
        k2_djs[i] = mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i];
        js3[i] = j0s[i] + 0.5 * dt * k2_djs[i];
    }

    // multiply new w/j after pushing 
    for (int i = 0; i < u_weights.size(); i++) {
        u_weights[i] = weights[i] * ws3[i];
        b_weights[i] = weights[i] * js3[i];
    }




    // -------------------------
    // Stage 3: k3_u1s, k3_u2s, k3_dws, k3_djs
    // -------------------------

    evaluate_u_field(k3_u1s, k3_u2s, xs3, ys3, u_weights, t + 0.5 * dt);

    // B evaluation
    b1s.assign(b1s.size(), 0.0);
    b2s.assign(b2s.size(), 0.0);
    evaluate_b_field(b1s, b2s, xs3, ys3, b_weights, t+ 0.5 * dt);

    // u1s_grad, u2s_grad, b1s_grad, b2s_grad evaluation
    u1s_grad_x.assign(xs.size(), 0.0);
    u1s_grad_y.assign(xs.size(), 0.0);
    evaluate_u1s_grad(u1s_grad_x, u1s_grad_y, xs3, ys3, u_weights, t+ 0.5 * dt);


    u2s_grad_x.assign(xs.size(), 0.0);
    u2s_grad_y.assign(xs.size(), 0.0);

    evaluate_u2s_grad(u2s_grad_x, u2s_grad_y, xs3, ys3, u_weights, t+ 0.5 * dt);

    b1s_grad_x.assign(xs.size(), 0.0);
    b1s_grad_y.assign(xs.size(), 0.0);

    evaluate_b1s_grad(b1s_grad_x, b1s_grad_y, xs3, ys3, b_weights, t+ 0.5 * dt);

    b2s_grad_x.assign(xs.size(), 0.0);
    b2s_grad_y.assign(xs.size(), 0.0);

    evaluate_b2s_grad(b2s_grad_x, b2s_grad_y, xs3, ys3, b_weights, t+ 0.5 * dt);

    // vortex_grad, j_grad evaluation
    vorticity_grad_x.assign(xs.size(), 0.0);
    vorticity_grad_y.assign(xs.size(), 0.0);

    evaluate_vorticity_grad(vorticity_grad_x, vorticity_grad_y, xs3, ys3, u_weights, t+ 0.5 * dt);

    j_grad_x.assign(xs.size(), 0.0);
    j_grad_y.assign(xs.size(), 0.0);

    evaluate_j_grad(j_grad_x, j_grad_y, xs3, ys3, b_weights, t+ 0.5 * dt);

    // vortex_laplacian, j_laplacian evaluation

    vorticity_laplacian.assign(xs.size(), 0.0);
    j_laplacian.assign(xs.size(), 0.0);

    vorticity_none_local.assign(xs.size(), 0.0);
    j_none_local.assign(xs.size(), 0.0);

    evaluate_vorticity_laplacian(vorticity_laplacian, vorticity_none_local, xs3, ys3, u_weights, t+ 0.5 * dt);
    evaluate_j_laplacian(j_laplacian, j_none_local, xs3, ys3, b_weights, t+ 0.5 * dt);

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

    // pushing 
    for (int i = 0; i < xs.size(); i++) {
        xs4[i] = xs[i] + dt * k3_u1s[i];
    }
    for (int i = 0; i < ys.size(); i++) {
        ys4[i] = ys[i] + dt * k3_u2s[i];
    }
    for (int i = 0; i < xs.size(); i++) {
        k3_dws[i] = nu * vorticity_laplacian[i] + B_dot_grad_j[i]; 
        ws4[i] = w0s[i] + dt * k3_dws[i];
    }
    for (int i = 0; i < xs.size(); i++) {
        k3_djs[i] = mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i];
        js4[i] = j0s[i] + dt * k3_djs[i];
    }

    // multiply new w/j after pushing 
    for (int i = 0; i < u_weights.size(); i++) {
        u_weights[i] = weights[i] * ws4[i];
        b_weights[i] = weights[i] * js4[i];
    }




    // -------------------------
    // Stage 4: k4_u1s, k4_u2s, k4_dws, k4_djs
    // -------------------------
    // evaluate_u_field(k4x, k4y, x4, y4, u_weights, t + dt);

    evaluate_u_field(k4_u1s, k4_u2s, xs4, ys4, u_weights, t + dt);

    // B evaluation
    b1s.assign(b1s.size(), 0.0);
    b2s.assign(b2s.size(), 0.0);
    evaluate_b_field(b1s, b2s, xs4, ys4, b_weights, t + dt);

    // u1s_grad, u2s_grad, b1s_grad, b2s_grad evaluation
    u1s_grad_x.assign(xs.size(), 0.0);
    u1s_grad_y.assign(xs.size(), 0.0);
    evaluate_u1s_grad(u1s_grad_x, u1s_grad_y, xs4, ys4, u_weights, t+ dt);


    u2s_grad_x.assign(xs.size(), 0.0);
    u2s_grad_y.assign(xs.size(), 0.0);

    evaluate_u2s_grad(u2s_grad_x, u2s_grad_y, xs4, ys4, u_weights, t+ dt);

    b1s_grad_x.assign(xs.size(), 0.0);
    b1s_grad_y.assign(xs.size(), 0.0);

    evaluate_b1s_grad(b1s_grad_x, b1s_grad_y, xs4, ys4, b_weights, t+ dt);

    b2s_grad_x.assign(xs.size(), 0.0);
    b2s_grad_y.assign(xs.size(), 0.0);

    evaluate_b2s_grad(b2s_grad_x, b2s_grad_y, xs4, ys4, b_weights, t+ dt);

    // vortex_grad, j_grad evaluation
    vorticity_grad_x.assign(xs.size(), 0.0);
    vorticity_grad_y.assign(xs.size(), 0.0);

    evaluate_vorticity_grad(vorticity_grad_x, vorticity_grad_y, xs4, ys4, u_weights, t+ dt);

    j_grad_x.assign(xs.size(), 0.0);
    j_grad_y.assign(xs.size(), 0.0);

    evaluate_j_grad(j_grad_x, j_grad_y, xs4, ys4, b_weights, t+ dt);

    // vortex_laplacian, j_laplacian evaluation

    vorticity_laplacian.assign(xs.size(), 0.0);
    j_laplacian.assign(xs.size(), 0.0);

    vorticity_none_local.assign(xs.size(), 0.0);
    j_none_local.assign(xs.size(), 0.0);

    evaluate_vorticity_laplacian(vorticity_laplacian, vorticity_none_local, xs4, ys4, u_weights, t+ dt);
    evaluate_j_laplacian(j_laplacian, j_none_local, xs4, ys4, b_weights, t+ dt);

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

    // source terms: 
    for (int i = 0; i < xs.size(); i++) {
        k4_dws[i] = nu * vorticity_laplacian[i] + B_dot_grad_j[i]; 
    }
    for (int i = 0; i < xs.size(); i++) {
        k4_djs[i] = mu * j_laplacian[i] + B_dot_grad_vorticity[i] + 2 * B_grad_x_dot_u2_grad[i] - 2 * B_grad_y_dot_u1_grad[i];
    }


    // -------------------------
    // Combine to update positions
    // x_{n+1} = x_n + dt/6*(k1 + 2k2 + 2k3 + k4)
    // (u1s, u2s, k1_dws, k1_djs), (k2_u1s, k2_u2s, k2_dws, k2_djs), (k3_u1s, k3_u2s, k3_dws, k3_djs), (k4_u1s, k4_u2s, k4_dws, k4_djs)
    // -------------------------
    for (int i = 0; i < xs_size; ++i) {
        xs[i] += (dt / 6.0) * (u1s[i] + 2.0 * k2_u1s[i] + 2.0 * k3_u1s[i] + k4_u1s[i]);
        ys[i] += (dt / 6.0) * (u2s[i] + 2.0 * k2_u2s[i] + 2.0 * k3_u2s[i] + k4_u2s[i]);
        w0s[i] += (dt / 6.0) * (k1_dws[i] + 2.0 * k2_dws[i] + 2.0 * k3_dws[i] + k4_dws[i]);
        j0s[i] += (dt / 6.0) * (k1_djs[i] + 2.0 * k2_djs[i] + 2.0 * k3_djs[i] + k4_djs[i]);
    }


    #ifdef DEBUG
    std::cout << "after rk4 u field evaluation (combined) - first 5:" << std::endl;
    for (int i = 0; i < std::min(5, xs_size); ++i) {
        std::cout << i
                  << " k1=(" << k1x[i] << "," << k1y[i] << ")"
                  << " k2=(" << k2x[i] << "," << k2y[i] << ")"
                  << " k3=(" << k3x[i] << "," << k3y[i] << ")"
                  << " k4=(" << k4x[i] << "," << k4y[i] << ")"
                  << std::endl;
    }
    #endif

    return 0;
}
