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
    std::vector<double> xs, ys, w0s, j0s, weights;
    std::vector<double> u_weights, b_weights;
    std::vector<double> u1s, u2s, b1s, b2s; // u1 velocity in x, u2 velocity in y; b1, b2: magnetic field in x,y

    std::vector <Panel> old_panels;
    std::vector<double> old_xs, old_ys, old_w0s, old_j0s;

    bool allow_boundary_extrapolation = false;
    double f_beyond_boundary = 0;
    bool do_unshear = false;
    bool sqrt_f = false;

    // time stepping parameters
    int iter_num;
    int num_steps;
    double dt;
    double t;
    int n_steps_remesh;
    int n_steps_diag;

    //field parameters

    int quad; // 0 for trap; 1 for simpsons

    Field* calculate_e;

    //profile parameters
    bool do_profile;
    std::vector<duration<double>> time_operations;
    std::vector<int> num_operations;

    // private functions
    int create_prerefined_mesh();
//     int create_prerefined_mesh_v_refinement();
    void refine_panels(std::function<double (double,double)> f, bool do_adaptive_refine);
//     void refine_panels_refine_v(std::function<double (double,double)> f, bool do_adaptive_refine);
    void test_panel(int panel_ind, bool verbose);

    int write_particles_to_file();
    int write_panels_to_file();
    int write_particles_to_file(bool pre_remesh);
    int write_panels_to_file(bool pre_remesh);

    public:
        AMRStructure();
        AMRStructure(std::string sim_dir, distribution* w0, distribution* j0,  
                int initial_height, int y_height, int max_height, 
                double x_min, double x_max, double y_min, double y_max, 
                int bcs, Field* calculate_e,
                int quad, int num_steps, double dt, 
                int n_steps_remesh,
                int n_steps_diag,
                bool do_adaptively_refine_vorticity, double amr_epsilons_vorticity,
                bool do_adaptively_refine_j, double amr_epsilons_j);


        // AMRStructure(std::string sim_dir, distribution* f0, 
        //         int initial_height, 
        //         double x_min, double x_max, double v_min, double v_max, 
        //         ElectricField* calculate_e, int num_steps, double dt);

        // AMRStructure(std::string sim_dir, distribution* f0, 
        //         double q, double m, 
        //         int initial_height, int v_height, int max_height, 
        //         double x_min, double x_max, double v_min, double v_max, 
        //         BoundaryConditions bcs,
        //         ElectricField* calculate_e, Quadrature quad, int num_steps, double dt, 
        //         bool do_adaptively_refine, std::vector<double>& amr_epsilons);

        // end constructors
        // destructor
        ~AMRStructure();
        // getters
        // std::vector<double> get_e();
        std::string get_sim_dir() const;
        // end getters

        // amr
        void generate_mesh(std::function<double (double,double)> w0, std::function<double (double,double)> j0,
                        bool do_adaptively_refine_vorticity, bool do_adaptively_refine_j, bool is_initial_step);
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
        int find_leaf_containing_xy_recursively(double &x, const double &y, bool& beyond_boundary, int panel_ind, bool verbose);
        int find_leaf_containing_point_from_neighbor(double& tx, double& ty, bool& beyond_boundary, int leaf_ind, std::set<int>& history, bool verbose);
        // // int find_leaf_containing();
        void interpolate_to_initial_xys(std::vector<double>& w0s, std::vector<double>& j0s, std::vector<double>& xs, std::vector<double>& ys, int nx, int ny,bool verbose);
        // double interpolate_from_mesh(double xs, double vs, bool verbose);
        // void interpolate_from_mesh(std::vector<double> &values, std::vector<double>& x, std::vector<double>& v, bool verbose);
        // void interpolate_from_mesh_slow(std::vector<double> &values, std::vector<double>& x, std::vector<double>& v, bool verbose);
        // double interpolate_from_panel(double x, double v, int panel_ind, bool use_limiter, bool verbose);
        void interpolate_from_panel_to_points(std::vector<double>& values, std::vector<double>& xs, std::vector<double>& ys,
                                                std::vector<int>& point_inds, int panel_ind, bool use_limiter, double limit_val);

        // field functions
        // void init_e();

        
        // u1 and u2, use u_weights
        int evaluate_u_field(std::vector<double>& u1s_local, std::vector<double>& u2s_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        // int evaluate_u2_field(std::vector<double>& u2s_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        // b1 and b2, use b_weights 
        int evaluate_b_field(std::vector<double>& b1s_local, std::vector<double>& b2s_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);
        // int evaluate_b2_field(std::vector<double>& b2s_local, std::vector<double>& xs_local,std::vector<double>& ys_local,std::vector<double>& ws_local,double t);



        // run step functions, all in run.cpp
        int run();
        int step();
        int euler();
        // void step(bool get_4th_e);
        // int rk4();

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
