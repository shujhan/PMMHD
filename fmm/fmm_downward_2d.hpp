#ifndef PMMHD_FMM_DOWNWARD_2D_HPP
#define PMMHD_FMM_DOWNWARD_2D_HPP

// BarytreeK's 2d downward pass (downward_pass.cpp / downward_pass_impl.hpp),
// copied unchanged except that the unreachable `throw` in the child switch is
// Kokkos::abort, because device code cannot throw. This lets the FMM run on
// the CUDA backend without modifying BarytreeK.

#include <Kokkos_Core.hpp>
#include "barytreek-config.h"
#include "structs.hpp"
#include "bli_impl.hpp"

struct fmm_parent_to_child_2d {
	view_panel_2d blfmm_panels;
	view_reall proxy_target_weights;
	int interp_deg;

	fmm_parent_to_child_2d(view_panel_2d& blfmm_panels_, view_reall& proxy_target_weights_, int interp_deg_) : 
					blfmm_panels(blfmm_panels_), proxy_target_weights(proxy_target_weights_), interp_deg(interp_deg_) {}

	KOKKOS_INLINE_FUNCTION
	void operator()(const int i, const int j) const {
		// interpolates proxy target weights of panel i to child panels
		int child;
		switch (j) {
			case 0:
				child = blfmm_panels(i).child1;
				break;
			case 1:
				child = blfmm_panels(i).child2;
				break;
			case 2:
				child = blfmm_panels(i).child3;
				break;
			case 3:
				child = blfmm_panels(i).child4;
				break;
			default:
				Kokkos::abort("j not between 0 and 3, parent to child");
		}
		real min_x_p = blfmm_panels(i).min_x;
		real max_x_p = blfmm_panels(i).max_x;
		real min_y_p = blfmm_panels(i).min_y;
		real max_y_p = blfmm_panels(i).max_y;
		real min_x_c = blfmm_panels(child).min_x;
		real max_x_c = blfmm_panels(child).max_x;
		real min_y_c = blfmm_panels(child).min_y;
		real max_y_c = blfmm_panels(child).max_y;
		real cheb_x[max_degree+1], cheb_y[max_degree+1], x_basis_vals[max_degree+1], y_basis_vals[max_degree+1];
		bli_points_shift(cheb_x, min_x_c, max_x_c, interp_deg);
		bli_points_shift(cheb_y, min_y_c, max_y_c, interp_deg);
		real x, y;
		int index1, index2;
		index1 = 0;
		for (int j = 0; j < interp_deg+1; j++) { // child x loop
			x = cheb_x[j];
			for (int k = 0; k < interp_deg+1; k++) { // child y loop
				y = cheb_y[k];
				interp_vals_bli(x_basis_vals, x, min_x_p, max_x_p, interp_deg);
				interp_vals_bli(y_basis_vals, y, min_y_p, max_y_p, interp_deg);
				index2 = 0;
				for (int j2 = 0; j2 < interp_deg+1; j2++) {
					for (int k2 = 0; k2 < interp_deg+1; k2++) {
						Kokkos::atomic_add(&proxy_target_weights(child,index1), x_basis_vals[j2]*y_basis_vals[k2]*proxy_target_weights(i,index2));
						index2 += 1;
					}
				}
				index1 += 1;
			}
		}
	}
};

struct fmm_leaf_to_point_2d {
	view_real xcos;
	view_real ycos;
	view_real soln;
	view_reall proxy_target_weights;
	view_panel_2d blfmm_panels;
	view_intt panel_points_inside;
	int interp_deg;

	fmm_leaf_to_point_2d(view_real& xcos_, view_real& ycos_, view_real& soln_, view_reall& proxy_target_weights_, view_panel_2d& blfmm_panels_, view_intt& panel_points_inside_, int interp_deg_) :
					xcos(xcos_), ycos(ycos_), soln(soln_), proxy_target_weights(proxy_target_weights_), blfmm_panels(blfmm_panels_), panel_points_inside(panel_points_inside_), interp_deg(interp_deg_) {}

	KOKKOS_INLINE_FUNCTION
	void operator()(const int i) const {
		int target_count = blfmm_panels(i).point_count;
		real min_x = blfmm_panels(i).min_x;
		real max_x = blfmm_panels(i).max_x;
		real min_y = blfmm_panels(i).min_y;
		real max_y = blfmm_panels(i).max_y;
		int i_t;
		real x, y, x_basis_vals[max_degree+1], y_basis_vals[max_degree+1];
		int index;
		for (int j = 0; j < target_count; j++) {
			i_t = panel_points_inside(i,j);
			x = xcos(i_t);
			y = ycos(i_t);
			interp_vals_bli(x_basis_vals, x, min_x, max_x, interp_deg);
			interp_vals_bli(y_basis_vals, y, min_y, max_y, interp_deg);
			index = 0;
			for (int k = 0; k < interp_deg+1; k++) {
				for (int l = 0; l < interp_deg+1; l++) {
					Kokkos::atomic_add(&soln(i_t), proxy_target_weights(i,index) * x_basis_vals[k] * y_basis_vals[l]);
					index += 1;
				}
			}
		}
	}
};

inline void fmm_downward_pass_2d(const RunConfig& run_config, const TreeInfo& tree_info, view_real& xcos, view_real& ycos, view_real& soln, view_panel_2d& blfmm_panels, view_reall& proxy_target_weights, view_intt& panel_points_inside) {
	int lb, ub;
	for (int i = 0; i < tree_info.levels-1; i++) {
		lb = tree_info.level_start[i];
		ub = tree_info.level_start[i+1];
		Kokkos::parallel_for(Kokkos::MDRangePolicy({lb, 0},{ub, 4}), fmm_parent_to_child_2d(blfmm_panels, proxy_target_weights, run_config.interp_degree));
	}
	lb = tree_info.level_start[tree_info.levels-1];
	ub = tree_info.level_start[tree_info.levels];
	Kokkos::parallel_for(Kokkos::RangePolicy(lb, ub), fmm_leaf_to_point_2d(xcos, ycos, soln, proxy_target_weights, blfmm_panels, panel_points_inside, run_config.interp_degree));
}

#endif
