#include "AMRStructure.hpp"

AMRStructure::AMRStructure() {}
AMRStructure::AMRStructure(std::string sim_dir, distribution* w0, distribution* j0, 
                            int initial_height, int y_height, int max_height, 
                            double x_min, double x_max, double y_min, double y_max, 
                            int bcs, Field* calculate_e,
                            int quad, int num_steps, double dt, 
                            int n_steps_remesh,
                            int n_steps_diag,
                            bool do_adaptively_refine_vorticity, double amr_epsilons_vorticity,
                            bool do_adaptively_refine_j, double amr_epsilons_j)
                           : w0(w0), j0(j0), 
                           initial_height(initial_height), y_height(y_height),
                           height(initial_height), max_height(max_height), 
                           x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), bcs(bcs),
                           calculate_e(calculate_e),
                           iter_num(0), num_steps(num_steps), dt(dt), t(0), n_steps_remesh(n_steps_remesh),
                           n_steps_diag(n_steps_diag), quad(quad),
                           is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
                           do_adaptively_refine_vorticity(do_adaptively_refine_vorticity),
                           amr_epsilons_vorticity(amr_epsilons_vorticity),
                           amr_epsilons_j(amr_epsilons_j),
                           do_adaptively_refine_j(do_adaptively_refine_j),
                           use_limiter(false), limit_val(0.0)
{
    // time_operations = std::vector<duration<double>>(last_time);
    // num_operations = std::vector<int> (last_time);

    this->sim_dir = sim_dir;
    Lx = x_max - x_min;
    Ly = y_max - y_min;
    npanels_x = int(pow(2, initial_height));
    npanels_y = int(pow(2, initial_height + y_height));
    initial_dx = Lx / npanels_x;
    initial_dy = Ly / npanels_y;

    bool is_initial_step = true;
    generate_mesh([&](double x, double y) { return (*w0)(x,y); }, [&](double x, double y) { return (*j0)(x,y); },do_adaptively_refine_vorticity, do_adaptively_refine_j, is_initial_step);
    // f_beyond_boundary = *std::min_element(fs.begin(), fs.end() );
    // cout << "extrapolating value is " << f_beyond_boundary << endl;
    u1s.assign(xs.size(), 0.0);
    u2s.assign(xs.size(), 0.0);
    b1s.assign(xs.size(), 0.0);
    b2s.assign(xs.size(), 0.0);
}




// AMRStructure::AMRStructure(std::string sim_dir, distribution* f0, //std::function<double (double,double)> f0, 
//                             int initial_height, 
//                             double x_min, double x_max, double v_min, double v_max, 
//                             ElectricField* calculate_e, int num_steps, double dt)
//                            : f0(f0), q(-1.0), qm(-1.0), 
//                            initial_height(initial_height) , v_height(0),
//                            height(initial_height), max_height(initial_height),
//                            x_min(x_min), x_max(x_max),
//                            v_min(v_min), v_max(v_max), bcs(periodic_bcs),
//                            iter_num(0), num_steps(num_steps), dt(dt),
//                            calculate_e(calculate_e), quad(trap),
//                            is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
//                            do_adaptively_refine(false),
//                            use_limiter(false), limit_val(0.0)
// {
//     time_operations = std::vector<duration<double>>(last_time);
//     num_operations = std::vector<int> (last_time);

//     this->sim_dir = sim_dir;
//     Lx = x_max - x_min;
//     Lv = v_max - v_min;
//     npanels_x = int(pow(2, initial_height));
//     npanels_v = int(pow(2, initial_height + v_height));
//     // create_prerefined_mesh();
//     bool is_initial_step = true;
    
//     generate_mesh([&](double x, double v) { return (*f0)(x,v); }, do_adaptively_refine, is_initial_step);
//     f_beyond_boundary = *std::min_element(fs.begin(), fs.end() );
//     cout << "extrapolating value is " << f_beyond_boundary << endl;
// }

// AMRStructure::AMRStructure(std::string sim_dir, distribution* f0, //std::function<double (double,double)> f0,
//                             int initial_height, int max_height, 
//                             double x_min, double x_max, double v_min, double v_max, 
//                             BoundaryConditions bcs,
//                             ElectricField* calculate_e, int num_steps, double dt, 
//                             bool do_adaptively_refine, std::vector<double>& amr_epsilons)
//                            : f0(f0), q(-1.0), qm(-1.0), 
//                            initial_height(initial_height), v_height(0),
//                            height(initial_height), max_height(max_height), 
//                            x_min(x_min), x_max(x_max), v_min(v_min), v_max(v_max), bcs(bcs),
//                            iter_num(0), num_steps(num_steps), dt(dt),
//                            calculate_e(calculate_e),
//                            is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
//                            do_adaptively_refine(do_adaptively_refine),
//                            use_limiter(false), limit_val(0.0)
// {
//     time_operations = std::vector<duration<double>>(last_time);
//     num_operations = std::vector<int> (last_time);

//     this->sim_dir = sim_dir;
//     Lx = x_max - x_min;
//     Lv = v_max - v_min;
//     npanels_x = int(pow(2, initial_height));
//     npanels_v = int(pow(2, initial_height + v_height));
//     initial_dx = Lx / npanels_x;
//     initial_dv = Lv / npanels_v;
//     this->amr_epsilons = amr_epsilons;

//     // create_prerefined_mesh();
//     bool is_initial_step = true;
//     generate_mesh([&](double x, double v) { return (*f0)(x,v); }, do_adaptively_refine, is_initial_step);
//     // calculate_e = new E_MQ_DirectSum(Lx, greens_epsilon);
//     f_beyond_boundary = *std::min_element(fs.begin(), fs.end() );
//     cout << "extrapolating value is " << f_beyond_boundary << endl;
// }

//end constructors

//destructor
AMRStructure::~AMRStructure() = default;

// getters
// std::vector<double> AMRStructure::get_e() { return es; };
std::string AMRStructure::get_sim_dir() const { return sim_dir; }
// end getters

// setters

// end setters

// void AMRStructure::add_time(ProfileTypes prof_type, duration<double> op_time) {
//     num_operations[prof_type] ++;
//     time_operations[prof_type] += op_time;
// }

