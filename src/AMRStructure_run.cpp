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
    // Classic RK4 for the Elsasser (q+/q-) transport, built as a 4-stage
    // generalisation of euler(), same idea as multi_species' rk4_step. Each
    // stage evaluates exactly the slopes euler() uses:
    //     q+ copy advects with (U - B), source +S
    //     q- copy advects with (U + B), source -S
    // multi_species can sample its field at scattered displaced points, so it
    // never remeshes between stages. Here S needs U/B gradients from grid FD,
    // so each stage state is formed like a full euler step: push the two copies
    // off the base grid, remesh them back onto the (fixed) grid via
    // interpolate_to_initial_xys, recover w/j and weights, refresh U/B fields.
    //
    // The grid (xs, ys, panels) is held FIXED across the 4 stages so the four
    // slope sets live at the same points and combine. AMR is applied once per
    // step at the step()-level remesh() (the grid we enter with was already
    // refined there); re-refining between stages would change the point count
    // and make the k1..k4 combine ill-defined. With AMR off the fixed grid is
    // just the uniform grid, so the same path serves both cases.
    const int N = (int)xs.size();
    const int nx_points = 2 * npanels_x + 1;
    const int ny_points = 2 * npanels_y + 1;
    xs_plus = xs;  ys_plus = ys;  xs_minus = xs;  ys_minus = ys;

    // base q+/q- at the fixed grid (held fixed; every stage state branches off these)
    std::vector<double> qp0 = q_plus;
    std::vector<double> qm0 = q_minus;

    // per-stage slopes at the fixed grid: dx/dt, dy/dt, dq/dt for each copy
    std::vector<double> kx_p(N), ky_p(N), kq_p(N);   // q+ copy
    std::vector<double> kx_m(N), ky_m(N), kq_m(N);   // q- copy

    // RK4 weighted-sum accumulators (the dt/6 is applied at the final combine)
    std::vector<double> sx_p(N, 0.0), sy_p(N, 0.0), sq_p(N, 0.0);
    std::vector<double> sx_m(N, 0.0), sy_m(N, 0.0), sq_m(N, 0.0);

    std::vector<double> S(N);

    // stage substep fraction (offset from base) and combine weight
    const double stage_c[4] = {0.0, 0.5, 0.5, 1.0};
    const double stage_b[4] = {1.0, 2.0, 2.0, 1.0};

    for (int stage = 0; stage < 4; ++stage) {
        // --- slopes at the current stage state ---
        // Fields (u1s/u2s/b1s/b2s) are already current on the fixed grid: from
        // the previous step's init_fields() at stage 0, and from this routine's
        // own init_fields() (below) at later stages.
        compute_source_S(xs, ys, w0s, j0s, t + stage_c[stage] * dt, S);
        for (int i = 0; i < N; ++i) {
            kx_p[i] = u1s[i] - b1s[i];  ky_p[i] = u2s[i] - b2s[i];  kq_p[i] =  S[i];
            kx_m[i] = u1s[i] + b1s[i];  ky_m[i] = u2s[i] + b2s[i];  kq_m[i] = -S[i];
        }

        // accumulate this stage into the RK4 sum
        const double b = stage_b[stage];
        for (int i = 0; i < N; ++i) {
            sx_p[i] += b * kx_p[i];  sy_p[i] += b * ky_p[i];  sq_p[i] += b * kq_p[i];
            sx_m[i] += b * kx_m[i];  sy_m[i] += b * ky_m[i];  sq_m[i] += b * kq_m[i];
        }

        // --- form the next stage state (push copies, remesh back, fields) ---
        if (stage < 3) {
            const double c = stage_c[stage + 1];   // 0.5, 0.5, 1.0
            for (int i = 0; i < N; ++i) {
                xs_plus[i]  = xs[i] + c * dt * kx_p[i];
                ys_plus[i]  = ys[i] + c * dt * ky_p[i];
                q_plus[i]   = qp0[i] + c * dt * kq_p[i];

                xs_minus[i] = xs[i] + c * dt * kx_m[i];
                ys_minus[i] = ys[i] + c * dt * ky_m[i];
                q_minus[i]  = qm0[i] + c * dt * kq_m[i];
            }

            old_panels = panels;
            old_xs = xs_plus;  old_ys = ys_plus;  old_q0s = q_plus;
            interpolate_to_initial_xys(q_plus, xs, ys, nx_points, ny_points);
            old_xs = xs_minus; old_ys = ys_minus; old_q0s = q_minus;
            interpolate_to_initial_xys(q_minus, xs, ys, nx_points, ny_points);

            for (int i = 0; i < N; ++i) {
                w0s[i] = 0.5 * (q_plus[i] + q_minus[i]);
                j0s[i] = 0.5 * (q_plus[i] - q_minus[i]);
                u_weights[i] = weights[i] * w0s[i];
                b_weights[i] = weights[i] * j0s[i];
            }
            init_fields(); 
        }
    }

    // --- final RK4 combine ---
    // Push the two copies off the base grid by dt/6 * (k1 + 2k2 + 2k3 + k4).
    // step() then calls remesh() (which AMR-refines and recovers w/j on the new
    // grid) followed by init_fields(), exactly as it does after euler().
    for (int i = 0; i < N; ++i) {
        xs_plus[i]  = xs[i] + (dt / 6.0) * sx_p[i];
        ys_plus[i]  = ys[i] + (dt / 6.0) * sy_p[i];
        q_plus[i]   = qp0[i] + (dt / 6.0) * sq_p[i];

        xs_minus[i] = xs[i] + (dt / 6.0) * sx_m[i];
        ys_minus[i] = ys[i] + (dt / 6.0) * sy_m[i];
        q_minus[i]  = qm0[i] + (dt / 6.0) * sq_m[i];
    }

    return 0;
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