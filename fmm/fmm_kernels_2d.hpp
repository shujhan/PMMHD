#ifndef PMMHD_FMM_KERNELS_2D_HPP
#define PMMHD_FMM_KERNELS_2D_HPP

// 2d interaction kernels (pp, pc, cp, cc) for BarytreeK, kept on the PMMHD side
// so BarytreeK itself is used unmodified (only its tree, traversal, upward and
// downward passes are compiled).
//
//   potential:   1/(4 pi) log(r^2 + eps^2) q           (PMMHD periodic_xy_potentials)
//                copied from BarytreeK poisson_fmm_impl.hpp (2d), renamed
//   Biot-Savart: same structure, two outputs           (PMMHD periodic_xy)

#include <Kokkos_Core.hpp>
#include "barytreek-config.h"
#include "structs.hpp"
#include <numbers>
#include "bli_impl.hpp"

// ---------------- potential ----------------

struct potential_pp_interaction_2d {
	view_real xcos_t;
	view_real ycos_t;
	view_real xcos_s;
	view_real ycos_s;
	view_real charges;
	view_real soln;
	view_intt panel_points_inside_target;
	view_intt panel_points_inside_source;
	view_interact interaction_list;
	view_panel_2d blfmm_panels_target;
	view_panel_2d blfmm_panels_source;
	real eps;

	potential_pp_interaction_2d(view_real& xcos_t_, view_real& ycos_t_, view_real& xcos_s_, view_real& ycos_s_, view_real& charges_, 
							view_real& soln_, view_intt& panel_points_inside_target_, view_intt& panel_points_inside_source_, 
							view_interact& interactions_, view_panel_2d& blfmm_panels_target_, view_panel_2d& blfmm_panels_source_, real eps_) :
							xcos_t(xcos_t_), ycos_t(ycos_t_), xcos_s(xcos_s_), ycos_s(ycos_s_), charges(charges_), soln(soln_), 
							panel_points_inside_target(panel_points_inside_target_), panel_points_inside_source(panel_points_inside_source_), 
							interaction_list(interactions_), blfmm_panels_target(blfmm_panels_target_), blfmm_panels_source(blfmm_panels_source_), eps(eps_) {}

	KOKKOS_INLINE_FUNCTION
	void operator()(const int i) const {
		int target_panel = interaction_list(i).target_panel;
		int source_panel = interaction_list(i).source_panel;
		int target_count = blfmm_panels_target(target_panel).point_count;
		int source_count = blfmm_panels_source(source_panel).point_count;
		real tx, ty, sx, sy, gfv;
		int i_t, i_s;
		real gfc = 1.0/(4.0*std::numbers::pi_v<real>);
		for (int j = 0; j < target_count; j++) {
			i_t = panel_points_inside_target(target_panel,j);
			tx = xcos_t(i_t);
			ty = ycos_t(i_t);
			for (int k = 0; k < source_count; k++) {
				i_s = panel_points_inside_source(source_panel,k);
				sx = xcos_s(i_s);
				sy = ycos_s(i_s);
				gfv = gfc*Kokkos::log((tx-sx)*(tx-sx)+(ty-sy)*(ty-sy)+eps*eps);
				Kokkos::atomic_add(&soln(i_t), gfv*charges(i_s));
			}
		}
	}
};

struct potential_pc_interaction_2d {
	view_real xcos_t;
	view_real ycos_t;
	view_real soln;
	view_intt panel_points_inside_target;
	view_reall proxy_source_weights;
	view_interact interaction_list;
	view_panel_2d blfmm_panels_target;
	view_panel_2d blfmm_panels_source;
	real eps;
	int interp_deg;

	potential_pc_interaction_2d(view_real& xcos_t_, view_real& ycos_t_, view_real& soln_, view_intt& panel_points_inside_target_, 
							view_reall& proxy_source_weights_, view_interact& interaction_list_, view_panel_2d& blfmm_panels_target_, 
							view_panel_2d& blfmm_panels_source_, real eps_, int interp_deg_) : xcos_t(xcos_t_), ycos_t(ycos_t_), soln(soln_), 
							panel_points_inside_target(panel_points_inside_target_), proxy_source_weights(proxy_source_weights_), 
							interaction_list(interaction_list_), blfmm_panels_target(blfmm_panels_target_), blfmm_panels_source(blfmm_panels_source_), 
							eps(eps_), interp_deg(interp_deg_) {}

	KOKKOS_INLINE_FUNCTION
	void operator()(const int i) const {
		int target_panel = interaction_list(i).target_panel;
		int source_panel = interaction_list(i).source_panel;
		int target_count = blfmm_panels_target(target_panel).point_count;
		real gfc = 1.0/(4.0*std::numbers::pi_v<real>);
		real min_x, max_x, x, min_y, max_y, y, cheb_x[max_degree+1], cheb_y[max_degree+1], tx, ty, sx, sy, gfv;
		min_x = blfmm_panels_source(source_panel).min_x;
		max_x = blfmm_panels_source(source_panel).max_x;
		min_y = blfmm_panels_source(source_panel).min_y;
		max_y = blfmm_panels_source(source_panel).max_y;
		bli_points_shift(cheb_x, min_x, max_x, interp_deg);
		bli_points_shift(cheb_y, min_y, max_y, interp_deg);
		int i_t, index;
		for (int l = 0; l < target_count; l++) {
			i_t = panel_points_inside_target(target_panel,l);
			tx = xcos_t(i_t);
			ty = ycos_t(i_t);
			index = 0;
			for (int j = 0; j < interp_deg+1; j++) { // x loop
				for (int k = 0; k < interp_deg+1; k++) { // y loop
					sx = cheb_x[j];
					sy = cheb_y[k];
					gfv = gfc*Kokkos::log((tx-sx)*(tx-sx)+(ty-sy)*(ty-sy)+eps*eps);
					Kokkos::atomic_add(&soln(i_t), gfv*proxy_source_weights(source_panel,index));
					index += 1;
				}
			}
		}
	}
};

struct potential_cp_interaction_2d {
	view_real xcos_s;
	view_real ycos_s;
	view_real charges;
	view_intt panel_points_inside_source;
	view_reall proxy_target_weights;
	view_interact interaction_list;
	view_panel_2d blfmm_panels_target;
	view_panel_2d blfmm_panels_source;
	real eps;
	int interp_deg;

	potential_cp_interaction_2d(view_real& xcos_s_, view_real& ycos_s_, view_real& charges_, view_intt& panel_points_inside_source_, 
							view_reall& proxy_target_weights_, view_interact& interaction_list_, view_panel_2d& blfmm_panels_target_, 
							view_panel_2d& blfmm_panels_source_, real eps_, int interp_deg_) : xcos_s(xcos_s_), ycos_s(ycos_s_), charges(charges_), 
							panel_points_inside_source(panel_points_inside_source_), proxy_target_weights(proxy_target_weights_), 
							interaction_list(interaction_list_), blfmm_panels_target(blfmm_panels_target_), blfmm_panels_source(blfmm_panels_source_), 
							eps(eps_), interp_deg(interp_deg_) {}

	KOKKOS_INLINE_FUNCTION
	void operator()(const int i) const {
		int target_panel = interaction_list(i).target_panel;
		int source_panel = interaction_list(i).source_panel;
		int source_count = blfmm_panels_source(source_panel).point_count;
		real gfc = 1.0/(4.0*std::numbers::pi_v<real>);
		real min_x, max_x, x, min_y, max_y, y, cheb_x[max_degree+1], cheb_y[max_degree+1], tx, ty, sx, sy, gfv;
		min_x = blfmm_panels_target(target_panel).min_x;
		max_x = blfmm_panels_target(target_panel).max_x;
		min_y = blfmm_panels_target(target_panel).min_y;
		max_y = blfmm_panels_target(target_panel).max_y;
		bli_points_shift(cheb_x, min_x, max_x, interp_deg);
		bli_points_shift(cheb_y, min_y, max_y, interp_deg);
		int i_s, index = 0;
		for (int j = 0; j < interp_deg+1; j++) { // x loop
			for (int k = 0; k < interp_deg+1; k++) { // y loop
				tx = cheb_x[j];
				ty = cheb_y[k];
				for (int l = 0; l < source_count; l++) {
					i_s = panel_points_inside_source(source_panel,l);
					sx = xcos_s(i_s);
					sy = ycos_s(i_s);
					gfv = gfc*Kokkos::log((tx-sx)*(tx-sx)+(ty-sy)*(ty-sy)+eps*eps);
					Kokkos::atomic_add(&proxy_target_weights(target_panel,index), gfv*charges(i_s));
				}
				index += 1;
			}
		}
	}
};

struct potential_cc_interaction_2d {
	view_reall proxy_target_weights;
	view_reall proxy_source_weights;
	view_interact interaction_list;
	view_panel_2d blfmm_panels_target;
	view_panel_2d blfmm_panels_source;
	real eps;
	int interp_deg;

	potential_cc_interaction_2d(view_reall& proxy_target_weights_, view_reall& proxy_source_weights_, view_interact& interaction_list_, 
							view_panel_2d& blfmm_panels_target_, view_panel_2d& blfmm_panels_source_, real eps_, int interp_deg_) :
							proxy_target_weights(proxy_target_weights_), proxy_source_weights(proxy_source_weights_), 
							interaction_list(interaction_list_), blfmm_panels_target(blfmm_panels_target_), 
							blfmm_panels_source(blfmm_panels_source_), eps(eps_), interp_deg(interp_deg_) {}

	KOKKOS_INLINE_FUNCTION
	void operator()(const int i) const {
		int target_panel = interaction_list(i).target_panel;
		int source_panel = interaction_list(i).source_panel;
		real gfc = 1.0/(4.0*std::numbers::pi_v<real>);
		real min_x_t = blfmm_panels_target(target_panel).min_x;
		real max_x_t = blfmm_panels_target(target_panel).max_x;
		real min_y_t = blfmm_panels_target(target_panel).min_y;
		real max_y_t = blfmm_panels_target(target_panel).max_y;
		real min_x_s = blfmm_panels_source(source_panel).min_x;
		real max_x_s = blfmm_panels_source(source_panel).max_x;
		real min_y_s = blfmm_panels_source(source_panel).min_y;
		real max_y_s = blfmm_panels_source(source_panel).max_y;
		real tx, ty, sx, sy, gfv;
		real cheb_x_t[max_degree+1], cheb_y_t[max_degree+1], cheb_x_s[max_degree+1], cheb_y_s[max_degree+1];
		bli_points_shift(cheb_x_t, min_x_t, max_x_t, interp_deg);
		bli_points_shift(cheb_y_t, min_y_t, max_y_t, interp_deg);
		bli_points_shift(cheb_x_s, min_x_s, max_x_s, interp_deg);
		bli_points_shift(cheb_y_s, min_y_s, max_y_s, interp_deg);
		int index_t, index_s;
		index_t = 0;
		for (int j1 = 0; j1 < interp_deg+1; j1++) { // target x loop
			for (int k1 = 0; k1 < interp_deg+1; k1++) { // target y loop
				tx = cheb_x_t[j1];
				ty = cheb_y_t[k1];
				index_s = 0;
				for (int j2 = 0; j2 < interp_deg+1; j2++) { // source x loop
					for (int k2 = 0; k2 < interp_deg+1; k2++) {	// source y loop
						sx = cheb_x_s[j2];
						sy = cheb_y_s[k2];
						gfv = gfc*Kokkos::log((tx-sx)*(tx-sx)+(ty-sy)*(ty-sy)+eps*eps);
						Kokkos::atomic_add(&proxy_target_weights(target_panel,index_t), gfv*proxy_source_weights(source_panel,index_s));
						index_s += 1;
					}
				}
				index_t += 1;
			}
		}
	}
};

inline void potential_fmm_interactions_2d(const RunConfig& run_config, view_real& xcos_t, view_real& ycos_t, view_real& xcos_s, view_real& ycos_s, view_real& charges, view_real& soln, view_intt& panel_points_inside_source, view_intt& panel_points_inside_target, 
								view_reall& proxy_source_weights, view_reall& proxy_target_weights, view_interact& pp_ints, view_interact& pc_ints, view_interact& cp_ints, view_interact& cc_ints, 
								view_panel_2d& blfmm_panels_source, view_panel_2d& blfmm_panels_target) {
	Kokkos::parallel_for(run_config.fmm_pp_count, potential_pp_interaction_2d(xcos_t, ycos_t, xcos_s, ycos_s, charges, soln, panel_points_inside_target, panel_points_inside_source, pp_ints, blfmm_panels_target, blfmm_panels_source, run_config.ker_eps));
	Kokkos::parallel_for(run_config.fmm_pc_count, potential_pc_interaction_2d(xcos_t, ycos_t, soln, panel_points_inside_target, proxy_source_weights, pc_ints, blfmm_panels_target, blfmm_panels_source, run_config.ker_eps, run_config.interp_degree));
	Kokkos::parallel_for(run_config.fmm_cp_count, potential_cp_interaction_2d(xcos_s, ycos_s, charges, panel_points_inside_source, proxy_target_weights, cp_ints, blfmm_panels_target, blfmm_panels_source, run_config.ker_eps, run_config.interp_degree));
	Kokkos::parallel_for(run_config.fmm_cc_count, potential_cc_interaction_2d(proxy_target_weights, proxy_source_weights, cc_ints, blfmm_panels_target, blfmm_panels_source, run_config.ker_eps, run_config.interp_degree));
	Kokkos::fence();
}

// ---------------- Biot-Savart velocity ----------------
// scalar charge q, two outputs (same kernel as PMMHD periodic_xy)
// vel_x += -1/(2 pi) * (ty-sy) / (r^2 + eps^2) * q
// vel_y +=  1/(2 pi) * (tx-sx) / (r^2 + eps^2) * q

struct velocity_pp_interaction_2d {
	view_real xcos_t;
	view_real ycos_t;
	view_real xcos_s;
	view_real ycos_s;
	view_real charges;
	view_real vel_x;
	view_real vel_y;
	view_intt panel_points_inside_target;
	view_intt panel_points_inside_source;
	view_interact interaction_list;
	view_panel_2d blfmm_panels_target;
	view_panel_2d blfmm_panels_source;
	real eps;

	velocity_pp_interaction_2d(view_real& xcos_t_, view_real& ycos_t_, view_real& xcos_s_, view_real& ycos_s_, view_real& charges_,
							view_real& vel_x_, view_real& vel_y_, view_intt& panel_points_inside_target_, view_intt& panel_points_inside_source_,
							view_interact& interactions_, view_panel_2d& blfmm_panels_target_, view_panel_2d& blfmm_panels_source_, real eps_) :
							xcos_t(xcos_t_), ycos_t(ycos_t_), xcos_s(xcos_s_), ycos_s(ycos_s_), charges(charges_), vel_x(vel_x_), vel_y(vel_y_),
							panel_points_inside_target(panel_points_inside_target_), panel_points_inside_source(panel_points_inside_source_),
							interaction_list(interactions_), blfmm_panels_target(blfmm_panels_target_), blfmm_panels_source(blfmm_panels_source_), eps(eps_) {}

	KOKKOS_INLINE_FUNCTION
	void operator()(const int i) const {
		int target_panel = interaction_list(i).target_panel;
		int source_panel = interaction_list(i).source_panel;
		int target_count = blfmm_panels_target(target_panel).point_count;
		int source_count = blfmm_panels_source(source_panel).point_count;
		real tx, ty, sx, sy, gfv;
		int i_t, i_s;
		real gfc = 1.0/(2.0*std::numbers::pi_v<real>);
		for (int j = 0; j < target_count; j++) {
			i_t = panel_points_inside_target(target_panel,j);
			tx = xcos_t(i_t);
			ty = ycos_t(i_t);
			real ux = 0;
			real uy = 0;
			for (int k = 0; k < source_count; k++) {
				i_s = panel_points_inside_source(source_panel,k);
				sx = xcos_s(i_s);
				sy = ycos_s(i_s);
				gfv = gfc*charges(i_s)/((tx-sx)*(tx-sx)+(ty-sy)*(ty-sy)+eps*eps);
				ux -= (ty-sy)*gfv;
				uy += (tx-sx)*gfv;
			}
			Kokkos::atomic_add(&vel_x(i_t), ux);
			Kokkos::atomic_add(&vel_y(i_t), uy);
		}
	}
};

struct velocity_pc_interaction_2d {
	view_real xcos_t;
	view_real ycos_t;
	view_real vel_x;
	view_real vel_y;
	view_intt panel_points_inside_target;
	view_reall proxy_source_weights;
	view_interact interaction_list;
	view_panel_2d blfmm_panels_target;
	view_panel_2d blfmm_panels_source;
	real eps;
	int interp_deg;

	velocity_pc_interaction_2d(view_real& xcos_t_, view_real& ycos_t_, view_real& vel_x_, view_real& vel_y_, view_intt& panel_points_inside_target_,
							view_reall& proxy_source_weights_, view_interact& interaction_list_, view_panel_2d& blfmm_panels_target_,
							view_panel_2d& blfmm_panels_source_, real eps_, int interp_deg_) : xcos_t(xcos_t_), ycos_t(ycos_t_), vel_x(vel_x_), vel_y(vel_y_),
							panel_points_inside_target(panel_points_inside_target_), proxy_source_weights(proxy_source_weights_),
							interaction_list(interaction_list_), blfmm_panels_target(blfmm_panels_target_), blfmm_panels_source(blfmm_panels_source_),
							eps(eps_), interp_deg(interp_deg_) {}

	KOKKOS_INLINE_FUNCTION
	void operator()(const int i) const {
		int target_panel = interaction_list(i).target_panel;
		int source_panel = interaction_list(i).source_panel;
		int target_count = blfmm_panels_target(target_panel).point_count;
		real gfc = 1.0/(2.0*std::numbers::pi_v<real>);
		real cheb_x[max_degree+1], cheb_y[max_degree+1], tx, ty, sx, sy, gfv;
		bli_points_shift(cheb_x, blfmm_panels_source(source_panel).min_x, blfmm_panels_source(source_panel).max_x, interp_deg);
		bli_points_shift(cheb_y, blfmm_panels_source(source_panel).min_y, blfmm_panels_source(source_panel).max_y, interp_deg);
		int i_t, index;
		for (int l = 0; l < target_count; l++) {
			i_t = panel_points_inside_target(target_panel,l);
			tx = xcos_t(i_t);
			ty = ycos_t(i_t);
			real ux = 0;
			real uy = 0;
			index = 0;
			for (int j = 0; j < interp_deg+1; j++) { // x loop
				for (int k = 0; k < interp_deg+1; k++) { // y loop
					sx = cheb_x[j];
					sy = cheb_y[k];
					gfv = gfc*proxy_source_weights(source_panel,index)/((tx-sx)*(tx-sx)+(ty-sy)*(ty-sy)+eps*eps);
					ux -= (ty-sy)*gfv;
					uy += (tx-sx)*gfv;
					index += 1;
				}
			}
			Kokkos::atomic_add(&vel_x(i_t), ux);
			Kokkos::atomic_add(&vel_y(i_t), uy);
		}
	}
};

struct velocity_cp_interaction_2d {
	view_real xcos_s;
	view_real ycos_s;
	view_real charges;
	view_intt panel_points_inside_source;
	view_reall proxy_target_weights_x;
	view_reall proxy_target_weights_y;
	view_interact interaction_list;
	view_panel_2d blfmm_panels_target;
	view_panel_2d blfmm_panels_source;
	real eps;
	int interp_deg;

	velocity_cp_interaction_2d(view_real& xcos_s_, view_real& ycos_s_, view_real& charges_, view_intt& panel_points_inside_source_,
							view_reall& proxy_target_weights_x_, view_reall& proxy_target_weights_y_, view_interact& interaction_list_,
							view_panel_2d& blfmm_panels_target_, view_panel_2d& blfmm_panels_source_, real eps_, int interp_deg_) :
							xcos_s(xcos_s_), ycos_s(ycos_s_), charges(charges_), panel_points_inside_source(panel_points_inside_source_),
							proxy_target_weights_x(proxy_target_weights_x_), proxy_target_weights_y(proxy_target_weights_y_),
							interaction_list(interaction_list_), blfmm_panels_target(blfmm_panels_target_), blfmm_panels_source(blfmm_panels_source_),
							eps(eps_), interp_deg(interp_deg_) {}

	KOKKOS_INLINE_FUNCTION
	void operator()(const int i) const {
		int target_panel = interaction_list(i).target_panel;
		int source_panel = interaction_list(i).source_panel;
		int source_count = blfmm_panels_source(source_panel).point_count;
		real gfc = 1.0/(2.0*std::numbers::pi_v<real>);
		real cheb_x[max_degree+1], cheb_y[max_degree+1], tx, ty, sx, sy, gfv;
		bli_points_shift(cheb_x, blfmm_panels_target(target_panel).min_x, blfmm_panels_target(target_panel).max_x, interp_deg);
		bli_points_shift(cheb_y, blfmm_panels_target(target_panel).min_y, blfmm_panels_target(target_panel).max_y, interp_deg);
		int i_s, index = 0;
		for (int j = 0; j < interp_deg+1; j++) { // x loop
			for (int k = 0; k < interp_deg+1; k++) { // y loop
				tx = cheb_x[j];
				ty = cheb_y[k];
				real ux = 0;
				real uy = 0;
				for (int l = 0; l < source_count; l++) {
					i_s = panel_points_inside_source(source_panel,l);
					sx = xcos_s(i_s);
					sy = ycos_s(i_s);
					gfv = gfc*charges(i_s)/((tx-sx)*(tx-sx)+(ty-sy)*(ty-sy)+eps*eps);
					ux -= (ty-sy)*gfv;
					uy += (tx-sx)*gfv;
				}
				Kokkos::atomic_add(&proxy_target_weights_x(target_panel,index), ux);
				Kokkos::atomic_add(&proxy_target_weights_y(target_panel,index), uy);
				index += 1;
			}
		}
	}
};

struct velocity_cc_interaction_2d {
	view_reall proxy_target_weights_x;
	view_reall proxy_target_weights_y;
	view_reall proxy_source_weights;
	view_interact interaction_list;
	view_panel_2d blfmm_panels_target;
	view_panel_2d blfmm_panels_source;
	real eps;
	int interp_deg;

	velocity_cc_interaction_2d(view_reall& proxy_target_weights_x_, view_reall& proxy_target_weights_y_, view_reall& proxy_source_weights_,
							view_interact& interaction_list_, view_panel_2d& blfmm_panels_target_, view_panel_2d& blfmm_panels_source_, real eps_, int interp_deg_) :
							proxy_target_weights_x(proxy_target_weights_x_), proxy_target_weights_y(proxy_target_weights_y_),
							proxy_source_weights(proxy_source_weights_), interaction_list(interaction_list_), blfmm_panels_target(blfmm_panels_target_),
							blfmm_panels_source(blfmm_panels_source_), eps(eps_), interp_deg(interp_deg_) {}

	KOKKOS_INLINE_FUNCTION
	void operator()(const int i) const {
		int target_panel = interaction_list(i).target_panel;
		int source_panel = interaction_list(i).source_panel;
		real gfc = 1.0/(2.0*std::numbers::pi_v<real>);
		real tx, ty, sx, sy, gfv;
		real cheb_x_t[max_degree+1], cheb_y_t[max_degree+1], cheb_x_s[max_degree+1], cheb_y_s[max_degree+1];
		bli_points_shift(cheb_x_t, blfmm_panels_target(target_panel).min_x, blfmm_panels_target(target_panel).max_x, interp_deg);
		bli_points_shift(cheb_y_t, blfmm_panels_target(target_panel).min_y, blfmm_panels_target(target_panel).max_y, interp_deg);
		bli_points_shift(cheb_x_s, blfmm_panels_source(source_panel).min_x, blfmm_panels_source(source_panel).max_x, interp_deg);
		bli_points_shift(cheb_y_s, blfmm_panels_source(source_panel).min_y, blfmm_panels_source(source_panel).max_y, interp_deg);
		int index_t, index_s;
		index_t = 0;
		for (int j1 = 0; j1 < interp_deg+1; j1++) { // target x loop
			for (int k1 = 0; k1 < interp_deg+1; k1++) { // target y loop
				tx = cheb_x_t[j1];
				ty = cheb_y_t[k1];
				real ux = 0;
				real uy = 0;
				index_s = 0;
				for (int j2 = 0; j2 < interp_deg+1; j2++) { // source x loop
					for (int k2 = 0; k2 < interp_deg+1; k2++) {	// source y loop
						sx = cheb_x_s[j2];
						sy = cheb_y_s[k2];
						gfv = gfc*proxy_source_weights(source_panel,index_s)/((tx-sx)*(tx-sx)+(ty-sy)*(ty-sy)+eps*eps);
						ux -= (ty-sy)*gfv;
						uy += (tx-sx)*gfv;
						index_s += 1;
					}
				}
				Kokkos::atomic_add(&proxy_target_weights_x(target_panel,index_t), ux);
				Kokkos::atomic_add(&proxy_target_weights_y(target_panel,index_t), uy);
				index_t += 1;
			}
		}
	}
};

inline void velocity_fmm_interactions_2d(const RunConfig& run_config, view_real& xcos_t, view_real& ycos_t, view_real& xcos_s, view_real& ycos_s, view_real& charges, view_real& vel_x, view_real& vel_y, 
								view_intt& panel_points_inside_source, view_intt& panel_points_inside_target, view_reall& proxy_source_weights, view_reall& proxy_target_weights_x, view_reall& proxy_target_weights_y, 
								view_interact& pp_ints, view_interact& pc_ints, view_interact& cp_ints, view_interact& cc_ints, view_panel_2d& blfmm_panels_source, view_panel_2d& blfmm_panels_target) {
	Kokkos::parallel_for(run_config.fmm_pp_count, velocity_pp_interaction_2d(xcos_t, ycos_t, xcos_s, ycos_s, charges, vel_x, vel_y, panel_points_inside_target, panel_points_inside_source, pp_ints, blfmm_panels_target, blfmm_panels_source, run_config.ker_eps));
	Kokkos::parallel_for(run_config.fmm_pc_count, velocity_pc_interaction_2d(xcos_t, ycos_t, vel_x, vel_y, panel_points_inside_target, proxy_source_weights, pc_ints, blfmm_panels_target, blfmm_panels_source, run_config.ker_eps, run_config.interp_degree));
	Kokkos::parallel_for(run_config.fmm_cp_count, velocity_cp_interaction_2d(xcos_s, ycos_s, charges, panel_points_inside_source, proxy_target_weights_x, proxy_target_weights_y, cp_ints, blfmm_panels_target, blfmm_panels_source, run_config.ker_eps, run_config.interp_degree));
	Kokkos::parallel_for(run_config.fmm_cc_count, velocity_cc_interaction_2d(proxy_target_weights_x, proxy_target_weights_y, proxy_source_weights, cc_ints, blfmm_panels_target, blfmm_panels_source, run_config.ker_eps, run_config.interp_degree));
	Kokkos::fence();
}

#endif
