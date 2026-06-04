#include "AMRStructure.hpp"

AMRStructure::AMRStructure() {}
AMRStructure::AMRStructure(std::string sim_dir, distribution* w0, distribution* j0, double nu, double mu,
                            int initial_height, int y_height, int max_height, 
                            double x_min, double x_max, double y_min, double y_max, 
                            int bcs, Field* calculate_e,
                            int quad, int num_steps, double dt, int method,
                            int n_steps_remesh,
                            int n_steps_diag,
                            bool do_adaptively_refine_vorticity, double amr_epsilons_vorticity,
                            bool do_adaptively_refine_j, double amr_epsilons_j, double greens_epsilon)
                           : w0(w0), j0(j0), nu(nu),mu(mu),
                           initial_height(initial_height), y_height(y_height),
                           height(initial_height), max_height(max_height), 
                           x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), bcs(bcs),
                           calculate_e(calculate_e),
                           iter_num(0), num_steps(num_steps), dt(dt), method(method), t(0), n_steps_remesh(n_steps_remesh),
                           n_steps_diag(n_steps_diag), quad(quad),
                           is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
                           do_adaptively_refine_vorticity(do_adaptively_refine_vorticity),
                           amr_epsilons_vorticity(amr_epsilons_vorticity),
                           amr_epsilons_j(amr_epsilons_j),
                           do_adaptively_refine_j(do_adaptively_refine_j), greens_epsilon(greens_epsilon),
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

    // generate all uniform_xs, uniform_ys to make sure all xs, ys are the same for different panels
    uniform_size = int(pow(2, max_height+1)) + 1;
    uniform_xs.assign(uniform_size, 0.0);
    uniform_ys.assign(uniform_size, 0.0);
    uniform_dx = Lx / (uniform_size - 1);
    uniform_dy = Ly / (uniform_size - 1);
    for (int i = 0; i < uniform_size; i++) {
        uniform_xs[i] = x_min + i * uniform_dx;
        uniform_ys[i] = y_min + i * uniform_dy;
    }

    
    bool is_initial_step = true;
    generate_mesh([&](double x, double y) { return (*w0)(x,y); }, [&](double x, double y) { return (*j0)(x,y); },do_adaptively_refine_vorticity, do_adaptively_refine_j, is_initial_step);
    w0s_beyond_boundary = *std::min_element(w0s.begin(), w0s.end() );
    cout << "w0 extrapolating value is " << w0s_beyond_boundary << endl;
    j0s_beyond_boundary = *std::min_element(j0s.begin(), j0s.end() );
    cout << "j0 extrapolating value is " << j0s_beyond_boundary << endl;
    u1s.assign(xs.size(), 0.0);
    u2s.assign(xs.size(), 0.0);
    b1s.assign(xs.size(), 0.0);
    b2s.assign(xs.size(), 0.0);
    // calculate fileds
    evaluate_u_field(u1s, u2s, xs, ys, u_weights, t);
    evaluate_b_field(b1s, b2s, xs, ys, b_weights, t);
}


//destructor
AMRStructure::~AMRStructure() = default;


std::string AMRStructure::get_sim_dir() const { return sim_dir; }


// void AMRStructure::add_time(ProfileTypes prof_type, duration<double> op_time) {
//     num_operations[prof_type] ++;
//     time_operations[prof_type] += op_time;
// }

