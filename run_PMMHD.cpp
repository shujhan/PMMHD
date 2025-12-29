#include <algorithm> // std::copy
#include <exception> // std::exception
#include <fstream>   // std::ifstream
#include <iterator> // std::ostream_iterator
#include <math.h> // M_PI, exp, cos, sin
#include <numeric> // std::inner_product
#include <stdio.h> // printf
#include <string>
#include <string.h> // atof
#include <iostream> // cout, endl
using std::cout;
using std::endl;
#include <stdexcept> // invalid_argument exception
#include <thread> // std::thread::hardware_concurrency



#include "boost/property_tree/ptree.hpp"        //property_tree::ptree, property_tree::read_json
#include "boost/property_tree/json_parser.hpp"

// #include "/sw/pkgs/arc/stacks/gcc/10.3.0/boost/1.78.0/include/boost/property_tree/ptree.hpp"
// #include "/sw/pkgs/arc/stacks/gcc/10.3.0/boost/1.78.0/include/boost/property_tree/json_parser.hpp"
//#include "/sw/pkgs/arc/stacks/gcc/10.3.0/boost/1.78.0/include/boost/property_tree/ptree_fwd.hpp"

namespace pt = boost::property_tree;

// #include "Panel.hpp"
#include "AMRStructure.hpp"
// #include "initial_distributions.hpp"

// #define DEBUG // for debugging purposes
// #define DEBUG2


int main(int argc, char** argv) {

    std::string sim_dir, input_deck;
    if (argc > 1) {
        sim_dir = argv[1];
    } else {
        sim_dir = "";
    }
    if (argc > 3) {
        input_deck = argv[2];
        input_deck = input_deck + argv[3];
    } else {
        if (argc > 2) {
            input_deck = sim_dir + argv[2];
        } else {
            input_deck = sim_dir + "deck.json";
        }
    }

    pt::ptree deck;
    try {
        pt::read_json(input_deck, deck);
    } catch(std::exception& err) {
        cout << "unable to open input deck" << endl;
        return 1;
    }

    // Get simulation box parameters 
    std::string project_name = deck.get<std::string>("project_name", "no_name_found");
    double x_min = deck.get<double>("xmin", 0.0), x_max = deck.get<double>("xmax", 1.0);
    double y_min = deck.get<double>("ymin", 0.0), y_max = deck.get<double>("ymax", 1.0);
    int bcs = deck.get<int>("bcs",0); // 0 for periodic bc
    double Lx = x_max - x_min;
    double Ly = y_max - y_min;
    int quad = deck.get<int>("quadrature",0); // 0 for trap rule 

    // get amr structure parameters
    int initial_height = deck.get<int>("initial_height", 2);
    int y_height = deck.get<int>("y_height",0);
    int max_height = deck.get<int>("max_height", initial_height);

    // get treecode parameters
    double greens_epsilon = deck.get<double>("greens_epsilon", 0.5);
    int use_treecode = deck.get<int>("use_treecode", 0); 
    double mac = deck.get<double>("mac", -1.0);
    int degree = deck.get<int>("degree", -1); 
    int max_source = deck.get<int>("max_source", 200); 
    int max_target = deck.get<int>("max_target", 200); 

    // get simulation parameters 
    int num_steps = deck.get<int>("num_steps", 10);//atoi(argv[19]);//120;
    int n_steps_remesh = deck.get<int>("remesh_period", 1); //atoi(argv[20]);
    int n_steps_diag = deck.get<int>("diag_period", 1); //atoi(argv[21]);
    double dt = deck.get<double> ("dt", 0.5); //atof(argv[22]);//0.5;

    // get vorticity and current density distribution parameters 

    pt::ptree &initial_list_deck = deck.get_child("initial_list");
    auto it = initial_list_deck.begin();

    // get vorticity
    pt::ptree &vorticity_dk = it->second;
    double kx_vorticity = 2.0 * M_PI / Lx * vorticity_dk.get<double>("normalized_wavenumber",1.0);
    double amp_vorticity = vorticity_dk.get<double>("amp", 0.0);
    int ics_type_vorticity = vorticity_dk.get<int>("ics_type", 1);
    bool do_adaptively_refine_vorticity = vorticity_dk.get<bool> ("adaptively_refine", false);
    double amr_epsilons_vorticity = vorticity_dk.get<double>("amr_epsilons",0.1);
    // get current density 
    ++it; 
    pt::ptree &current_density_dk = it->second;
    double kx_j = 2.0 * M_PI / Lx * current_density_dk.get<double>("normalized_wavenumber",1.0);
    double amp_j = current_density_dk.get<double>("amp", 0.0);
    int ics_type_j = current_density_dk.get<int>("ics_type", 1);
    bool do_adaptively_refine_j = current_density_dk.get<bool> ("adaptively_refine", false);
    double amr_epsilons_j = current_density_dk.get<double>("amr_epsilons",0.1);



    // create distribution for vorticity 
    distribution* w0;
    switch (ics_type_vorticity)
    {
        case 1: // for vorticity 
            w0 = new w0_uniform();
            break;
        default:
            cout << "Using default initial conditions, all 1s" << endl;
            w0 = new w0_uniform();
            break;
    }
    
    // create distribution for current density 
    distribution* j0;
    switch (ics_type_j)
    {
        case 1: // for current density
            j0 = new j0_uniform();
            break;

        case 2: 
            j0 = new j0_non_uniform();
            break;

        default:
            cout << "Using default initial conditions, all 1s" << endl;
            j0 = new j0_uniform();
            break;
    }


    // create field solver 
    Field* calculate_field;
    if (use_treecode > 0) {
        // calculate_e = new E_MQ_Treecode(Lx, greens_epsilon, mac, degree, max_source, max_target);
        calculate_field = new U_DirectSum(Lx, greens_epsilon);
    }
    else {
        calculate_field = new U_DirectSum(Lx, greens_epsilon);
    }

    cout << "============================" << endl;
    cout << "Running a FARRSIGHT simulation" << endl;
    cout << "sim dir: " << sim_dir << endl;
    cout << "deck found in: " << input_deck << endl;
    cout << x_min << " <= x <= " << x_max << endl;
    cout << y_min << " <= y <= " << y_max << endl;
    switch (bcs) {
        case 1 : cout << "Using open boundary conditions" << endl;
            break;
        default : // periodic
            cout << "Using periodic boundary conditions" << endl;
            break;
    }

    cout << "vorticity ics type: " << ics_type_vorticity << endl;
    cout << "current density ics type: " << ics_type_j << endl;

    cout << "height " << initial_height << ", y height " << y_height << ", max height " << max_height << endl;
    switch (quad) {
        case 1 : cout << "Using Simpson's rule" << endl;
            break;
        default : cout << "Using trap rule" << endl;
            break;
    }
    cout << "green's epsilon = " << greens_epsilon << endl;
    cout << "Taking " << num_steps << " steps with dt = " << dt << endl;
    cout << "Remesh every " << n_steps_remesh << " step(s), diagnostic dump every " << n_steps_diag << " step(s)" << endl;

    if (use_treecode > 0) { 
        cout << "Using treecode with mac " << mac << " and degree " << degree << endl;
    } else {
        cout << "using direct sum" << endl;
    }

    cout << "k_vorticity = " << kx_vorticity << ", amp_vorticity = " << amp_vorticity <<  endl;
    cout << "k_current_density = " << kx_j << ", amp_current_density = " << amp_j <<  endl;

    if (do_adaptively_refine_vorticity) {
        cout << "Adaptively refining for vorticity, to height at most " << max_height << endl;
        cout << "Refinement epsilon(s) : " << amr_epsilons_vorticity << endl;
    } else {
        cout <<"Not adaptively refining for vorticity." << endl;
    }
    if (do_adaptively_refine_j) {
        cout << "Adaptively refining for current density, to height at most " << max_height << endl;
        cout << "Refinement epsilon(s) : " << amr_epsilons_j << endl;
    } else {
        cout <<"Not adaptively refining for current density." << endl;
    }
    cout << "============================" << endl;

    auto sim_start = high_resolution_clock::now();

    AMRStructure amr{sim_dir, w0, j0,
                initial_height, y_height, max_height,
                x_min, x_max, y_min, y_max, bcs,
                calculate_field, quad, num_steps, dt,
                n_steps_remesh, n_steps_diag,
                do_adaptively_refine_vorticity, amr_epsilons_vorticity,
                do_adaptively_refine_j, amr_epsilons_j};
                

    // amr.init_e();
    amr.write_to_file();

    // print AMR structure info 
    cout << amr << endl;

    // run simulation
    amr.run();


    delete w0;
    delete j0;
    delete calculate_field;
    return 0;
}

