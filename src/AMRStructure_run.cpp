#include "AMRStructure.hpp"
#include <iomanip>

int AMRStructure::run() {
    while (iter_num < num_steps) {
        step();
    }

    write_to_file();
    return 0;
}

int AMRStructure::step() {
    if (iter_num % n_steps_diag == 0) {
        write_to_file();
    }
    iter_num += 1;
    std::cout << "step " << iter_num << std::endl;

    if (method > 0) {
        euler();
    }
    else {
        rk4();
    }

    t += dt;

    remesh();
    init_fields();

    return 0;
}


int AMRStructure::init_fields() {
    // initial U/B evaluation on the structured mesh; safe to call only after the
    // periodizer has been set (for periodic bcs). Fields are recomputed every step.
    u1s.assign(xs.size(), 0.0);
    u2s.assign(xs.size(), 0.0);
    b1s.assign(xs.size(), 0.0);
    b2s.assign(xs.size(), 0.0);
    evaluate_u_field(u1s, u2s, xs, ys, u_weights, t);
    evaluate_b_field(b1s, b2s, xs, ys, b_weights, t);

    for (size_t i = 0; i < xs.size(); ++i) {
        b1s[i] += B0x;
        b2s[i] += B0y;
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
    const int N = (int)xs.size();

    // RHS at (X_n, t_n): fields + S on the current (structured) mesh
    std::vector<double> S(N);

    compute_source_S(xs, ys, w0s, j0s, t, S);

    // two copies start co-located on the mesh, advect with U∓B
    xs_plus = xs;  ys_plus = ys;  xs_minus = xs;  ys_minus = ys;
    for (int i = 0; i < N; ++i) {
        xs_plus[i]  += dt * (u1s[i] - b1s[i]);
        ys_plus[i]  += dt * (u2s[i] - b2s[i]);
        q_plus[i]   += dt * S[i];

        xs_minus[i] += dt * (u1s[i] + b1s[i]);
        ys_minus[i] += dt * (u2s[i] + b2s[i]);
        q_minus[i]  -= dt * S[i];
    }
    return 0;
}




int AMRStructure::rk4() {
    // NOT IMPLEMENTED for the Elsasser (q+/q-) scheme yet.
    // The body below is the old in-place omega/j RK4: it calls compute_rhs_state
    // (now replaced by compute_source_S) and never fills xs_plus/q_plus/..., so it is
    // inconsistent with euler() and remesh(). Disabled for now -- use method>0 (euler).
    std::cerr << "rk4() is not implemented for the q+/q- scheme yet; "
                 "run with method>0 to use euler." << std::endl;
    return 1;

#if 0
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
#endif
}



// Nodal x/y derivatives of `field` on the 9 points of leaf `panel`.
// Uses a centered cross-panel difference when the neighbor on that side is a
// same-level leaf (so on a uniform mesh this reproduces the old result exactly),
// and a second-order one-sided in-panel difference otherwise -- i.e. at any
// coarse-fine interface or missing/refined neighbor, where reading a neighbor's
// node with this panel's spacing is geometrically wrong.
void AMRStructure::leaf_field_gradient(const std::vector<double>& field,
                                       const Panel* panel,
                                       double dx[9], double dy[9])
{
    const int* P = panel->point_inds;
    double f[9];
    for (int ii = 0; ii < 9; ++ii) { f[ii] = field[P[ii]]; }

    double hx = xs[P[3]] - xs[P[0]];   // column spacing (points 0,3,6 increase in x)
    double hy = ys[P[1]] - ys[P[0]];   // row spacing    (points 0,1,2 increase in y)

    auto same_level_leaf = [&](int nbr) -> bool {
        if (nbr < 0) { return false; }            // -1 (coarser) / -2 (domain boundary)
        const Panel& q = panels[nbr];
        return (q.child_inds_start < 0) && (q.level == panel->level);
    };

    // ---- x-derivatives ----
    dx[3] = (f[6] - f[0]) / (2*hx);               // middle column: in-panel centered
    dx[4] = (f[7] - f[1]) / (2*hx);
    dx[5] = (f[8] - f[2]) / (2*hx);

    if (same_level_leaf(panel->left_nbr_ind)) {
        const int* L = panels[panel->left_nbr_ind].point_inds;
        dx[0] = (f[3] - field[L[3]]) / (2*hx);
        dx[1] = (f[4] - field[L[4]]) / (2*hx);
        dx[2] = (f[5] - field[L[5]]) / (2*hx);
    } else {                                      // forward one-sided, in-panel
        dx[0] = (-3*f[0] + 4*f[3] - f[6]) / (2*hx);
        dx[1] = (-3*f[1] + 4*f[4] - f[7]) / (2*hx);
        dx[2] = (-3*f[2] + 4*f[5] - f[8]) / (2*hx);
    }

    if (same_level_leaf(panel->right_nbr_ind)) {
        const int* R = panels[panel->right_nbr_ind].point_inds;
        dx[6] = (field[R[3]] - f[3]) / (2*hx);
        dx[7] = (field[R[4]] - f[4]) / (2*hx);
        dx[8] = (field[R[5]] - f[5]) / (2*hx);
    } else {                                      // backward one-sided, in-panel
        dx[6] = (3*f[6] - 4*f[3] + f[0]) / (2*hx);
        dx[7] = (3*f[7] - 4*f[4] + f[1]) / (2*hx);
        dx[8] = (3*f[8] - 4*f[5] + f[2]) / (2*hx);
    }

    // ---- y-derivatives ----
    dy[1] = (f[2] - f[0]) / (2*hy);               // middle row: in-panel centered
    dy[4] = (f[5] - f[3]) / (2*hy);
    dy[7] = (f[8] - f[6]) / (2*hy);

    if (same_level_leaf(panel->bottom_nbr_ind)) {
        const int* B = panels[panel->bottom_nbr_ind].point_inds;
        dy[0] = (f[1] - field[B[1]]) / (2*hy);
        dy[3] = (f[4] - field[B[4]]) / (2*hy);
        dy[6] = (f[7] - field[B[7]]) / (2*hy);
    } else {                                      // forward one-sided, in-panel
        dy[0] = (-3*f[0] + 4*f[1] - f[2]) / (2*hy);
        dy[3] = (-3*f[3] + 4*f[4] - f[5]) / (2*hy);
        dy[6] = (-3*f[6] + 4*f[7] - f[8]) / (2*hy);
    }

    if (same_level_leaf(panel->top_nbr_ind)) {
        const int* T = panels[panel->top_nbr_ind].point_inds;
        dy[2] = (field[T[1]] - f[1]) / (2*hy);
        dy[5] = (field[T[4]] - f[4]) / (2*hy);
        dy[8] = (field[T[7]] - f[7]) / (2*hy);
    } else {                                      // backward one-sided, in-panel
        dy[2] = (3*f[2] - 4*f[1] + f[0]) / (2*hy);
        dy[5] = (3*f[5] - 4*f[4] + f[3]) / (2*hy);
        dy[8] = (3*f[8] - 4*f[7] + f[6]) / (2*hy);
    }
}

int AMRStructure::compute_source_S( std::vector<double>& xs_in, std::vector<double>& ys_in,
                            std::vector<double>& w0s_in, std::vector<double>& j0s_in,
                            double t_in, std::vector<double>& S_out
) {
    cout << "enter compute_source" << endl;

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
    // For each leaf panel, fill nodal gradients of U and B (and w0/j0) using an
    // AMR-aware stencil: centered across a same-level leaf neighbor (identical to
    // the uniform-mesh result), one-sided in-panel at coarse-fine interfaces where
    // a uniform-spacing cross-panel difference is geometrically invalid and was
    // injecting spurious O(1) gradients into S at every refinement boundary.
    for (int panel_ind = 0; panel_ind < panels.size(); panel_ind++) {
        Panel* panel = &(panels[panel_ind]);
        if (panel->child_inds_start > -1) { continue; } // leaves only

        const int* P = panel->point_inds;
        double dxv[9], dyv[9];

        leaf_field_gradient(w0s, panel, dxv, dyv);
        for (int ii = 0; ii < 9; ++ii) { vorticity_grad_x[P[ii]] = dxv[ii]; vorticity_grad_y[P[ii]] = dyv[ii]; }

        leaf_field_gradient(j0s, panel, dxv, dyv);
        for (int ii = 0; ii < 9; ++ii) { j_grad_x[P[ii]] = dxv[ii]; j_grad_y[P[ii]] = dyv[ii]; }

        leaf_field_gradient(u1s, panel, dxv, dyv);
        for (int ii = 0; ii < 9; ++ii) { u1s_grad_x[P[ii]] = dxv[ii]; u1s_grad_y[P[ii]] = dyv[ii]; }

        leaf_field_gradient(u2s, panel, dxv, dyv);
        for (int ii = 0; ii < 9; ++ii) { u2s_grad_x[P[ii]] = dxv[ii]; u2s_grad_y[P[ii]] = dyv[ii]; }

        leaf_field_gradient(b1s, panel, dxv, dyv);
        for (int ii = 0; ii < 9; ++ii) { b1s_grad_x[P[ii]] = dxv[ii]; b1s_grad_y[P[ii]] = dyv[ii]; }

        leaf_field_gradient(b2s, panel, dxv, dyv);
        for (int ii = 0; ii < 9; ++ii) { b2s_grad_x[P[ii]] = dxv[ii]; b2s_grad_y[P[ii]] = dyv[ii]; }
    }

    // S term (same as two-panel): built directly from the U and B gradients.
    // Note: the old tendency RHS (dw0s_dt / dj0s_dt with nu*Laplacian, B.grad terms)
    // is intentionally gone -- in the Elsasser split, transport carries only the ideal
    // coupling through S, and dissipation is applied separately after remesh.
    S_out.assign(xs.size(), 0.0);
    for (int i = 0; i < xs.size(); ++i) {
        S_out[i] = 2.0 * (b1s_grad_x[i]*u2s_grad_x[i] + b2s_grad_x[i]*u2s_grad_y[i])
                 - 2.0 * (b1s_grad_y[i]*u1s_grad_x[i] + b2s_grad_y[i]*u1s_grad_y[i]);
    }

    return 0;
}