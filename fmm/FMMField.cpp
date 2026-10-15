#include <iostream>
#include <cstdlib>
#include <cmath>

#include <Kokkos_Core.hpp>

#include "barytreek-config.h"
#include "structs.hpp"
#include "tree_construction.hpp"
#include "interaction_list.hpp"
#include "fmm_traversal_2d.hpp"
#include "upward_pass.hpp"
#include "fmm_kernels_2d.hpp"
#include "fmm_downward_2d.hpp"

#include "FMMField.hpp"

using namespace std;

void fmm_initialize(int& argc, char** argv) {
    if (!Kokkos::is_initialized()) {
        Kokkos::initialize(argc, argv);
    }
}

void fmm_finalize() {
    if (Kokkos::is_initialized() && !Kokkos::is_finalized()) {
        Kokkos::finalize();
    }
}

U_FMM::U_FMM() {}
U_FMM::U_FMM(double x_min, double L, double epsilon, double mac, int degree, int max_source)
    : x_min(x_min), L(L), epsilon(epsilon), mac(mac), degree(degree), max_source(max_source), mode(periodic_xy) {}
U_FMM::~U_FMM() = default;

void* U_FMM::operator new(size_t size) { return ::operator new(size); }
void U_FMM::operator delete(void* p) { ::operator delete(p); }

void U_FMM::print_field_obj() {
    cout << "[U_FMM]" << endl;
    cout << "  x_min=" << x_min << " L=" << L << " epsilon=" << epsilon << " theta=" << mac << " degree=" << degree
         << " cluster_size=" << max_source << endl;
}

void U_FMM::set_mode(KernelMode m) {
    mode = m;
}

void U_FMM::operator()(double* e1s, double* e2s,
                       double* x_vals, int nx,
                       double* y_vals, double* q_ws, int ny)
{
    if (mode != periodic_xy && mode != periodic_xy_potentials &&
        mode != original && mode != original_potentials) {
        cout << "U_FMM: only periodic_xy(_potentials) and original(_potentials) are supported" << endl;
        std::abort();
    }
    // channel kernel (periodic in x, open in y) for mode original(_potentials)
    const int periodic_x = (mode == original || mode == original_potentials) ? 1 : 0;
    const bool velocity = (mode == periodic_xy || mode == original);

    // run parameters (deck values, no namelist)
    RunConfig run_config;
    run_config.dim = 2;
    run_config.interp_degree = (degree > max_degree) ? max_degree : degree;
    run_config.interp_point_count = (run_config.interp_degree + 1) * (run_config.interp_degree + 1);
    run_config.fmm_theta = mac;
    run_config.fmm_cluster_thresh = max_source;
    run_config.ker_eps = epsilon;
    run_config.mpi_id = 0;
    run_config.mpi_p = 1;

    // targets: first nx points, sources: all ny points
    view_real_host x_co_t ("target x coordinates", nx);
    view_real_host y_co_t ("target y coordinates", nx);
    view_real_host x_co_s ("source x coordinates", ny);
    view_real_host y_co_s ("source y coordinates", ny);
    view_real_host charges ("point charges", ny);
    view_real_host sols_x ("solution x", nx);
    view_real_host sols_y ("solution y", nx);

    // channel kernel: wrap x into [x_min, x_min + L), as U_Treecode does with fmod
    auto wrap_x = [&](double x) {
        if (!periodic_x) return x;
        double xw = x_min + std::fmod(x - x_min, L);
        if (xw < x_min) xw += L;
        if (xw >= x_min + L) xw -= L;
        return xw;
    };
    for (int i = 0; i < nx; i++) {
        x_co_t(i) = wrap_x(x_vals[i]);
        y_co_t(i) = y_vals[i];
    }
    for (int i = 0; i < ny; i++) {
        x_co_s(i) = wrap_x(x_vals[i]);
        y_co_s(i) = y_vals[i];
        charges(i) = q_ws[i];
    }

    // trees
    TreeInfo tree_info_target, tree_info_source;
    view_panel_2d_host blfmm_panels_target ("target blfmm tree panels", 1);
    view_panel_2d_host blfmm_panels_source ("source blfmm tree panels", 1);
    view_int_host point_leaf_panel_target ("target point leaf panel indices", nx);
    view_intt_host panel_points_inside_target ("target leaf panels contained points", 1, 1);
    view_int_host point_leaf_panel_source ("source point leaf panel indices", ny);
    view_intt_host panel_points_inside_source ("source leaf panels contained points", 1, 1);

    blfmm_tree_construction_2d(run_config, tree_info_target, x_co_t, y_co_t, blfmm_panels_target, point_leaf_panel_target, panel_points_inside_target);
    blfmm_tree_construction_2d(run_config, tree_info_source, x_co_s, y_co_s, blfmm_panels_source, point_leaf_panel_source, panel_points_inside_source);

    // interaction lists
    view_interact_host interaction_list ("blfmm interactions", 1);
    fmm_dual_tree_traversal_2d(run_config, blfmm_panels_target, blfmm_panels_source, interaction_list, periodic_x, L);

    view_interact_host pp_interactions ("pp interactions", run_config.fmm_pp_count);
    view_interact_host pc_interactions ("pc interactions", run_config.fmm_pc_count);
    view_interact_host cp_interactions ("cp interactions", run_config.fmm_cp_count);
    view_interact_host cc_interactions ("cc interactions", run_config.fmm_cc_count);
    split_interactions(interaction_list, pp_interactions, pc_interactions, cp_interactions, cc_interactions);

    cout << "FMM (" << Kokkos::DefaultExecutionSpace::name() << "): targets " << nx << ", sources " << ny
         << ", target panels " << tree_info_target.panel_count << ", source panels " << tree_info_source.panel_count
         << ", interactions pp/pc/cp/cc " << run_config.fmm_pp_count << "/" << run_config.fmm_pc_count
         << "/" << run_config.fmm_cp_count << "/" << run_config.fmm_cc_count << endl;

    // device copies
    view_real d_x_co_t ("device target x coordinates", nx);
    view_real d_y_co_t ("device target y coordinates", nx);
    view_real d_x_co_s ("device source x coordinates", ny);
    view_real d_y_co_s ("device source y coordinates", ny);
    view_real d_charges ("device charges", ny);
    view_real d_sol_x ("device solution x", nx);
    view_real d_sol_y ("device solution y", nx);
    view_interact d_pp_ints ("device pp interactions", run_config.fmm_pp_count);
    view_interact d_pc_ints ("device pc interactions", run_config.fmm_pc_count);
    view_interact d_cp_ints ("device cp interactions", run_config.fmm_cp_count);
    view_interact d_cc_ints ("device cc interactions", run_config.fmm_cc_count);
    view_panel_2d d_blfmm_panels_target ("device target blfmm panels", tree_info_target.panel_count);
    view_panel_2d d_blfmm_panels_source ("device source blfmm panels", tree_info_source.panel_count);
    view_int d_point_leaf_panel_source ("device source point leaf panel indices", ny);
    view_intt d_panel_points_inside_target ("device target leaf panels contained points", tree_info_target.panel_count, run_config.fmm_cluster_thresh);
    view_intt d_panel_points_inside_source ("device source leaf panels contained points", tree_info_source.panel_count, run_config.fmm_cluster_thresh);

    Kokkos::deep_copy(d_x_co_t, x_co_t);
    Kokkos::deep_copy(d_y_co_t, y_co_t);
    Kokkos::deep_copy(d_x_co_s, x_co_s);
    Kokkos::deep_copy(d_y_co_s, y_co_s);
    Kokkos::deep_copy(d_charges, charges);
    Kokkos::deep_copy(d_blfmm_panels_target, blfmm_panels_target);
    Kokkos::deep_copy(d_blfmm_panels_source, blfmm_panels_source);
    Kokkos::deep_copy(d_point_leaf_panel_source, point_leaf_panel_source);
    Kokkos::deep_copy(d_panel_points_inside_target, panel_points_inside_target);
    Kokkos::deep_copy(d_panel_points_inside_source, panel_points_inside_source);
    Kokkos::deep_copy(d_pp_ints, pp_interactions);
    Kokkos::deep_copy(d_pc_ints, pc_interactions);
    Kokkos::deep_copy(d_cp_ints, cp_interactions);
    Kokkos::deep_copy(d_cc_ints, cc_interactions);

    view_reall proxy_source_weights ("proxy source weights", tree_info_source.panel_count, run_config.interp_point_count);
    view_reall proxy_target_weights_x ("proxy target weights x", tree_info_target.panel_count, run_config.interp_point_count);
    view_reall proxy_target_weights_y ("proxy target weights y", tree_info_target.panel_count, run_config.interp_point_count);
    Kokkos::deep_copy(proxy_source_weights, 0.0);
    Kokkos::deep_copy(proxy_target_weights_x, 0.0);
    Kokkos::deep_copy(proxy_target_weights_y, 0.0);
    Kokkos::deep_copy(d_sol_x, 0.0);
    Kokkos::deep_copy(d_sol_y, 0.0);

    // upward pass: charges -> source proxy weights
    upward_pass_2d(run_config, tree_info_source, d_x_co_s, d_y_co_s, d_charges, d_blfmm_panels_source, proxy_source_weights, d_point_leaf_panel_source);

    if (velocity) {
        // Biot-Savart (free-space or channel kernel, see fmm_kernels_2d.hpp)
        velocity_fmm_interactions_2d(run_config, d_x_co_t, d_y_co_t, d_x_co_s, d_y_co_s, d_charges, d_sol_x, d_sol_y,
                                        d_panel_points_inside_source, d_panel_points_inside_target, proxy_source_weights,
                                        proxy_target_weights_x, proxy_target_weights_y,
                                        d_pp_ints, d_pc_ints, d_cp_ints, d_cc_ints, d_blfmm_panels_source, d_blfmm_panels_target, periodic_x, L);
        // downward pass once per component
        fmm_downward_pass_2d(run_config, tree_info_target, d_x_co_t, d_y_co_t, d_sol_x, d_blfmm_panels_target,
                         proxy_target_weights_x, d_panel_points_inside_target);
        fmm_downward_pass_2d(run_config, tree_info_target, d_x_co_t, d_y_co_t, d_sol_y, d_blfmm_panels_target,
                         proxy_target_weights_y, d_panel_points_inside_target);
    }
    else {
        // potential (free-space or channel kernel, see fmm_kernels_2d.hpp)
        potential_fmm_interactions_2d(run_config, d_x_co_t, d_y_co_t, d_x_co_s, d_y_co_s, d_charges, d_sol_x,
                                    d_panel_points_inside_source, d_panel_points_inside_target, proxy_source_weights,
                                    proxy_target_weights_x, d_pp_ints, d_pc_ints, d_cp_ints, d_cc_ints,
                                    d_blfmm_panels_source, d_blfmm_panels_target, periodic_x, L);
        fmm_downward_pass_2d(run_config, tree_info_target, d_x_co_t, d_y_co_t, d_sol_x, d_blfmm_panels_target,
                         proxy_target_weights_x, d_panel_points_inside_target);
    }
    Kokkos::fence();

    Kokkos::deep_copy(sols_x, d_sol_x);
    Kokkos::deep_copy(sols_y, d_sol_y);
    for (int i = 0; i < nx; i++) {
        e1s[i] = sols_x(i);
        e2s[i] = sols_y(i);   // zero in potentials mode
    }

    free(tree_info_target.level_start);
    free(tree_info_source.level_start);
}