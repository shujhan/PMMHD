#ifndef AMRSTRUCTURE_HPP
#define AMRSTRUCTURE_HPP

#include <algorithm>            // std::sort, std::copy, std::find
#include <assert.h>             // assert
#include <chrono>               // high_resolution _clock, duration_cast, microseconds
using namespace std::chrono;
#include <functional>           // std::logical_and
#include <fstream>
#include <iostream>             // std::cout, std::endl
#include <iterator>             // std::ostream_iterator for vector printing
#include <math.h>
#include <numeric>              // std::iota, std::accumulate
#include <omp.h>
#include <set>                  // std::set, set.find
#include <stdexcept>            // exceptions
#include <stdio.h>              // printf
#include <string> 
#include <vector>               // std::vector

// external libraries to include
#include <Eigen/Dense>
using namespace Eigen;

// amr in-project dependencies
#include "initial_distributions.hpp"
#include "Panel.hpp"
#include "FieldStructure.hpp"
#include "Periodizer.hpp"

// named values for the integer `bcs` and `quad` members
// (bcs == 0 is doubly-periodic, bcs == 1 is open-in-y; quad 0 = trap, 1 = simpsons)
enum BCType   { periodic_bcs = 0, open_bcs = 1 };
enum QuadType { trap = 0, simpsons = 1 };



struct AMRStructure {
    std::string sim_dir;
    // domain parameters
    double Lx, Ly;
    double x_min, x_max;
    double y_min, y_max;
    int bcs;

    // initial condition
    distribution* w0;
    distribution* j0;

    // mesh parameters
    int initial_height;
    int y_height;
    int height;
    int max_height;
    int npanels_x, npanels_y;
    double initial_dx, initial_dy;


    bool is_initial_mesh_set;
    bool need_further_refinement;
    int minimum_unrefined_index;

    // for vorticity and current density 
    bool do_adaptively_refine_vorticity; 
    bool do_adaptively_refine_j;
    double amr_epsilons_vorticity;
    double amr_epsilons_j;


    std::vector <Panel> panels;
    std::vector <int> leaf_inds;

    //source terms
    int uniform_size;
    double uniform_dx, uniform_dy;
    std::vector<double> uniform_xs, uniform_ys;
    std::vector<double> xs, ys, w0s, j0s, weights;


    // transient Elsässer state used only inside step()
    std::vector<double> q_plus, q_minus;          // q± sampled at mesh points
    std::vector<double> xs_plus, ys_plus;         // advected positions of the q+ copy
    std::vector<double> xs_minus, ys_minus;       // advected positions of the q− copy
    std::vector<double> source_S;                 // the S term at mesh points

    // the two deformed Lagrangian copies saved for remeshing
    // (both reuse old_panels connectivity — only coordinates differ)
    std::vector <Panel> old_panels;
    std::vector<double> old_xs, old_ys, old_w0s, old_j0s, old_q0s;
    std::vector<double> old_xs_plus, old_ys_plus, old_q_plus;
    std::vector<double> old_xs_minus, old_ys_minus, old_q_minus;
    double B0x = 0.0, B0y = 0.0;        // optional uniform guide field
    Periodizer* periodizer = nullptr;   



    std::vector<double> u_weights, b_weights;
    std::vector<double> u1s, u2s, b1s, b2s; // u1 velocity in x, u2 velocity in y; b1, b2: magnetic field in x,y
    std::vector<double> u1s_grad_x, u1s_grad_y;
    std::vector<double> u2s_grad_x, u2s_grad_y;
    std::vector<double> b1s_grad_x, b1s_grad_y;
    std::vector<double> b2s_grad_x, b2s_grad_y;
    std::vector<double> vorticity_grad_x, vorticity_grad_y;
    std::vector<double> j_grad_x, j_grad_y;
    std::vector<double> vorticity_laplacian;
    std::vector<double> j_laplacian;


    //source terms calculation
    double nu,mu; // nu: fluid viscosity, mu: resistivity
    std::vector<double> B_dot_grad_j;
    std::vector<double> B_dot_grad_vorticity;

    std::vector<double> B_grad_x_dot_u2_grad;
    std::vector<double> B_grad_y_dot_u1_grad;

    bool allow_boundary_extrapolation = false;
    double w0_beyond_boundary = 0;
    double j0_beyond_boundary = 0;
    double q0_beyond_boundary = 0;
    bool do_unshear = false;
    bool sqrt_f = false;

    // time stepping parameters
    int iter_num;
    int num_steps;
    double dt;
    double t;
    int n_steps_remesh;
    int n_steps_diag;
    int method;

    //field parameters

    int quad; // 0 for trap; 1 for simpsons

    Field* calculate_e;
    double greens_epsilon;

    //profile parameters
    bool do_profile;
    std::vector<duration<double>> time_operations;
    std::vector<int> num_operations;

    // private functions
    int create_prerefined_mesh(bool is_initial_step);
    void refine_panels(std::function<double (double,double)> f, bool do_adaptive_refine, bool is_initial_step);
    void refine_panels_refine_v(std::function<double (double,double)> f, bool do_adaptive_refine,  bool is_initial_step);
    void test_panel(int panel_ind, bool verbose);

    int write_particles_to_file();
    int write_panels_to_file();
    int write_particles_to_file(bool pre_remesh);
    int write_panels_to_file(bool pre_remesh);

    public:
        AMRStructure();
        AMRStructure(std::string sim_dir, distribution* w0, distribution* j0, double nu, double mu,
                int initial_height, int y_height, int max_height, 
                double x_min, double x_max, double y_min, double y_max, 
                int bcs, Field* calculate_e,
                int quad, int num_steps, double dt, int method,
                int n_steps_remesh,
                int n_steps_diag,
                bool do_adaptively_refine_vorticity, double amr_epsilons_vorticity,
                bool do_adaptively_refine_j, double amr_epsilons_j, double greens_epsilon);

        // end constructors
        // destructor
        ~AMRStructure();
        // getters
        // std::vector<double> get_e();
        std::string get_sim_dir() const;
        // end getters

        // amr
        // void refine_panels_refine_v(std::function<double (double,double)> w0, std::function<double (double,double)> j0,
        //                 bool do_adaptively_refine_vorticity, bool do_adaptively_refine_j, bool is_initial_step);
        void set_leaves_weights();
        void recursively_set_leaves_weights(int panel_ind);

        // // remesh
        // void copy_to_old();
        // void reset_mesh();
        void remesh();

        // // interpolation functions
        bool use_limiter;
        double limit_val;
        void shift_xs(std::vector<double>& shifted_xs, const std::vector<double>& xs, const std::vector<double>& ys);
        int find_leaf_containing_xy_recursively(double &x, double &y, bool& beyond_boundary, int panel_ind);
        int find_leaf_containing_point_from_neighbor(double& tx, double& ty, bool& beyond_boundary, int leaf_ind, std::set<int>& history);
        // // int find_leaf_containing();
        void interpolate_to_initial_xys(std::vector<double>& q0s, std::vector<double>& xs, std::vector<double>& ys, int nx, int ny);
        // double interpolate_from_mesh(double xs, double vs, bool verbose);
        // void interpolate_from_mesh(std::vector<double> &values, std::vector<double>& x, std::vector<double>& v, bool verbose);
        // void interpolate_from_mesh_slow(std::vector<double> &values, std::vector<double>& x, std::vector<double>& v, bool verbose);
        // double interpolate_from_panel(double x, double v, int panel_ind, bool use_limiter, bool verbose);
        void interpolate_from_panel_to_points(std::vector<double>& values_q0, std::vector<double>& xs, std::vector<double>& ys,
                                                std::vector<int>& point_inds, int panel_ind, bool use_limiter, double limit_val);

        // field functions
        int init_fields();

        
        // u1 and u2, use u_weights
        int evaluate_u_field(std::vector<double>& u1s_local, std::vector<double>& u2s_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        // b1 and b2, use b_weights 
        int evaluate_b_field(std::vector<double>& b1s_local, std::vector<double>& b2s_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        
        // // gradients evaluation
        // int evaluate_u1s_grad(std::vector<double>& u1s_grad_x_local, std::vector<double>& u1s_grad_y_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        // int evaluate_u2s_grad(std::vector<double>& u2s_grad_x_local, std::vector<double>& u2s_grad_y_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        // int evaluate_b1s_grad(std::vector<double>& b1s_grad_x_local, std::vector<double>& b1s_grad_y_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        // int evaluate_b2s_grad(std::vector<double>& b2s_grad_x_local, std::vector<double>& b2s_grad_y_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        
        // // vorticity_gradient and j_gradient
        // int evaluate_vorticity_grad(std::vector<double>& vorticity_grad_x_local, std::vector<double>& vorticity_grad_y_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        // int evaluate_j_grad(std::vector<double>& j_grad_x_local, std::vector<double>& j_grad_y_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        
        // // vorticity_laplacian and j_laplacian
        // int evaluate_vorticity_laplacian(std::vector<double>& vorticity_laplacian_local, std::vector<double>& vorticity_none_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        // int evaluate_j_laplacian(std::vector<double>& j_laplacian_local, std::vector<double>& j_none_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        



        // run step functions, all in run.cpp
        int run();
        int step();
        int euler();
        // void step(bool get_4th_e);
        int rk4();
        // int compute_rhs_state(std::vector<double>& xs_in, std::vector<double>& ys_in, std::vector<double>& w0s_in,
        //         std::vector<double>& j0s_in, double t_in, std::vector<double>& dxs_dt, std::vector<double>& dys_dt, std::vector<double>& dw0s_dt, std::vector<double>& dj0s_dt);


        int  compute_source_S(std::vector<double>& xs_in, std::vector<double>& ys_in,
                            std::vector<double>& w0s_in, std::vector<double>& j0s_in,
                            double t_in, std::vector<double>& S_out,
                            std::vector<double>& u1_out, std::vector<double>& u2_out,
                            std::vector<double>& b1_out, std::vector<double>& b2_out);

        void set_periodizer(Periodizer* p) { periodizer = p; }

        // io
        friend std::ostream& operator<<(std::ostream& os, const AMRStructure& amr);
        void print_amr();
        int write_to_file();
        int write_to_file(bool pre_remesh);
        // void print_panel_points();

        // profiling
        // void add_time(ProfileTypes prof_type, duration<double> op_time);
        // void print_times();
};


#endif /* AMRSTRUCTURE_HPP */
