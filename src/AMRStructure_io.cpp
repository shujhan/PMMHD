#include "AMRStructure.hpp"

int AMRStructure::write_particles_to_file() {
    bool pre_remesh = false;
    return write_particles_to_file(pre_remesh);
}

int AMRStructure::write_particles_to_file(bool pre_remesh) {
    std::ofstream xs_file;
    std::ofstream ys_file;
    std::ofstream w0s_file;
    std::ofstream j0s_file;
    std::ofstream uweights_file;
    std::ofstream u1s_file;
    std::ofstream u2s_file;
    std::ofstream b1s_file;
    std::ofstream b2s_file;


    std::string remesh_str = "";
    if (pre_remesh) {
        remesh_str = "preremesh_";
    }
    xs_file.open(sim_dir + "simulation_output/xs/xs_" + remesh_str + std::to_string(iter_num), std::ios::out | std::ios::binary); 
    ys_file.open(sim_dir + "simulation_output/ys/ys_" + remesh_str  + std::to_string(iter_num), std::ios::out | std::ios::binary); 
    w0s_file.open(sim_dir + "simulation_output/w0s/w0s_" + remesh_str  + std::to_string(iter_num), std::ios::out | std::ios::binary); 
    j0s_file.open(sim_dir + "simulation_output/j0s/j0s_" + remesh_str  + std::to_string(iter_num), std::ios::out | std::ios::binary);
    uweights_file.open(sim_dir + "simulation_output/uweights/uweights_" + remesh_str  + std::to_string(iter_num), std::ios::out | std::ios::binary); 
    u1s_file.open(sim_dir + "simulation_output/u1s/u1s_" + remesh_str  + std::to_string(iter_num), std::ios::out | std::ios::binary); 
    u2s_file.open(sim_dir + "simulation_output/u2s/u2s_" + remesh_str  + std::to_string(iter_num), std::ios::out | std::ios::binary); 
    b1s_file.open(sim_dir + "simulation_output/b1s/b1s_" + remesh_str  + std::to_string(iter_num), std::ios::out | std::ios::binary); 
    b2s_file.open(sim_dir + "simulation_output/b2s/b2s_" + remesh_str  + std::to_string(iter_num), std::ios::out | std::ios::binary); 

    std::cout << "#xs " << xs.size() << std::endl;
    std::cout << "#ys " << ys.size() << std::endl;
    std::cout << "#w0s " << w0s.size() << std::endl;
    std::cout << "#j0s " << j0s.size() << std::endl;
    std::cout << "#uweights " << u_weights.size() << std::endl;
    std::cout << "#u1s " << u1s.size() << std::endl;
    std::cout << "#u2s " << u2s.size() << std::endl;
    std::cout << "#b1s " << b1s.size() << std::endl;
    std::cout << "#b2s " << b2s.size() << std::endl;

    if (!xs_file | !ys_file | !w0s_file | !j0s_file | !uweights_file | !u1s_file | !u2s_file | !b1s_file | !b2s_file) {
        cout << "Unable to open step " << iter_num << " particle data files" << endl;
        return 1;
    }

    // assert(xs.size() == vs.size() && vs.size() == fs.size() && fs.size() == q_ws.size() && q_ws.size() == es.size());
    for (int ii = 0; ii < xs.size(); ++ii) {
        double x = xs[ii];
        double y = ys[ii];
        double w0 = w0s[ii];
        double j0 = j0s[ii];
        double uweight = u_weights[ii];
        double u1 = u1s[ii];
        double u2 = u2s[ii];
        double b1 = b1s[ii];
        double b2 = b2s[ii];
        xs_file.write((char *) &x, sizeof(double));
        ys_file.write((char *) &y, sizeof(double));
        w0s_file.write((char *) &w0, sizeof(double));
        j0s_file.write((char *) &j0, sizeof(double));
        uweights_file.write((char *) &uweight, sizeof(double));
        u1s_file.write((char *) &u1, sizeof(double));
        u2s_file.write((char *) &u2, sizeof(double));
        b1s_file.write((char *) &b1, sizeof(double));
        b2s_file.write((char *) &b2, sizeof(double));
    }

    if (!xs_file.good() | !ys_file.good() | !w0s_file.good() | !j0s_file.good() | !uweights_file.good() | !u1s_file.good()
        | !u2s_file.good() | !b1s_file.good() | !b2s_file.good()) {
        cout << "Error occurred writing step " << iter_num << " particle data files." << endl;
        return 1;
    }
    xs_file.close();
    ys_file.close();
    w0s_file.close();
    j0s_file.close();
    uweights_file.close();
    u1s_file.close();
    u2s_file.close();
    b1s_file.close();
    b2s_file.close();
    // cout << "Successfully wrote step " << iter_num << " particle data files" << endl;

    return 0;
}

int AMRStructure::write_panels_to_file() {
    bool pre_remesh = false;
    return write_panels_to_file(pre_remesh);
}
int AMRStructure::write_panels_to_file(bool pre_remesh) {
    std::ofstream panel_file;
    std::string remesh_str = "";
    if (pre_remesh) {
        remesh_str = "preremesh_";
    }

    panel_file.open(sim_dir + "simulation_output/panels/leaf_point_inds_" + remesh_str+ std::to_string(iter_num), std::ios::out | std::ios::binary);
    
    if (!panel_file) {
        cout << "Unable to open step " << iter_num << " panel data files" << endl;
        return 1;
    }

    for (int ii=0; ii < leaf_inds.size(); ++ii ) {
        const int *inds = panels[leaf_inds[ii]].point_inds;
        panel_file.write( (char*) inds, 9*sizeof(int));
    }

    panel_file.close();
    if (!panel_file.good() ) {
        cout << "Error writing step " << iter_num << " panel data files" << endl;
        return 1;
    }

    // cout << "Successfully wrote step " << iter_num << " panel data files" << endl;

    return 0;
}

int AMRStructure::write_to_file() { 
    bool pre_remesh=false;
    return write_to_file(pre_remesh);
}
int AMRStructure::write_to_file(bool pre_remesh) {
    write_particles_to_file(pre_remesh);
    write_panels_to_file(pre_remesh);
    return 0;
}



void AMRStructure::print_amr() {
    printf("This mesh structure has %lu panels.\n", panels.size() );
    for (const Panel& panel : panels) {
        panel.print_panel();
    }
}
// AMRStructure_io.cpp
std::ostream& operator<<(std::ostream& os, const AMRStructure& amr) {
    os << "=================" << endl;
    os << "This is an AMR structure" << endl;
    os << "=================" << endl;
    os << "Computational domain: (x,y) in [" << amr.x_min << ", " << amr.x_max << "]x[" << amr.y_min << ", " << amr.y_max << "]" << endl; 
    os << "-----------------" << endl;
    os << "Initial height: " << amr.initial_height << ", y height: " << amr.y_height << ", height: " << amr.height << ", max height: " << amr.max_height << endl;
    if (amr.do_adaptively_refine_vorticity) { os << "This structure is adaptively refined for vorticity" << endl;}
    else { os << "This structure is NOT adaptively refined for vorticity" << endl;}

    if (amr.do_adaptively_refine_j) { os << "This structure is adaptively refined for current density" << endl;}
    else { os << "This structure is NOT adaptively refined for current density" << endl;}

    if (amr.need_further_refinement) { os << "Structure is not done being refined" << endl;}
    else { os << "Structure has reached full refinement" << endl; }


    os << "---------" << endl;
    os << "num steps: " <<  amr.iter_num << ", dt: " << amr.dt << endl;
    os << "---------" << endl;


    // os << "Panels: " << endl << "==============" << endl;
    // std::copy(amr.panels.begin(), amr.panels.end(), std::ostream_iterator<Panel>(os));
    // os <<  "==============" << endl;

    // os << "Point data" << endl << "==============" << endl;
    // os << "xs, size = " << amr.xs.size() << endl;
    // std::copy(amr.xs.begin(), amr.xs.end(), std::ostream_iterator<double>(os, " "));
    // os << endl << "vs, size = " << amr.vs.size() << endl;
    // std::copy(amr.vs.begin(), amr.vs.end(), std::ostream_iterator<double>(os, " "));
    // os << endl << "fs, size = " << amr.fs.size() << endl;
    // std::copy(amr.fs.begin(), amr.fs.end(), std::ostream_iterator<double>(os, " "));
    // os << endl << "qw, size = " << amr.q_ws.size() << endl;
    // std::copy(amr.q_ws.begin(), amr.q_ws.end(), std::ostream_iterator<double>(os, " "));
    // os << endl << "es, size = " << amr. es.size() << endl;
    // std::copy(amr.es.begin(), amr.es.end(), std::ostream_iterator<double>(os, " "));

    return os;
}

// void AMRStructure::print_panel_points() {
//     auto panel_it = panels.begin();

//     for (int ii = 0; ii < panels.size(); ii++) {
//         cout << "==============" << endl;
//         cout << "Panel " << ii << " points" << endl;
//         cout << "xs: ";
//         for (int jj = 0; jj < 9; jj++) {
//             cout << xs[panel_it->point_inds[jj]] << ", ";
//         }
//         cout << endl << "vs: ";
//         for (int jj = 0; jj < 9; jj++) {
//             cout << vs[panel_it->point_inds[jj]] << ", ";
//         }
//         cout << endl << "fs: ";
//         for (int jj = 0; jj < 9; jj++) {
//             cout << fs[panel_it->point_inds[jj]] << ", ";
//         }
//         cout << endl << "qw: ";
//         for (int jj = 0; jj < 9; jj++) {
//             cout << q_ws[panel_it->point_inds[jj]] << ", ";
//         }
//         cout << endl;
//         panel_it++;
//     }
// }

// void AMRStructure::print_times() {
//     cout << "========= timing ==========" << endl;
//     cout << "Total sim time " << time_operations[sim_time].count() << " seconds" << endl;
//     cout << ". Total step time " << time_operations[step_time].count() << " seconds" << endl;
//     cout << "- " << num_operations[step_time] << " steps taken at ";
//     cout << time_operations[step_time].count() / num_operations[step_time] << " seconds per step" << endl;
//     cout << ".. Total field time " << time_operations[field_time].count() << " seconds" << endl;
//     cout << "- " << num_operations[field_time] << " field evaluations made at ";
//     cout << time_operations[field_time].count() / num_operations[field_time] << " seconds per step" << endl;
//     cout << ". Total remeshing time " << time_operations[remesh_time].count() << " seconds" << endl;
//     cout << "- " << num_operations[remesh_time] << " remeshings made at ";
//     cout << time_operations[remesh_time].count() / num_operations[remesh_time] << " seconds per remeshing" << endl;
//     cout << ".. Initial tree building  + refinement time " << time_operations[tree_build_time].count() << " seconds (this includes the interpolations for mesh refinement)" << endl;
//     cout << "- " << num_operations[tree_build_time] << " not sure what these ops all are, at ";
//     cout << time_operations[tree_build_time].count() / num_operations[tree_build_time] << " seconds per op" << endl;
//     if (do_adaptively_refine && v_height + initial_height < max_height) {
//         cout <<"... Total panel-testing-in-amr time " << time_operations[amr_test_time].count() << " seconds" << endl;
//         cout << "- " << num_operations[amr_test_time] << " test sessions made at ";
//         cout << time_operations[amr_test_time].count() / num_operations[amr_test_time] << " seconds per test session" << endl;
//         cout << "... Total amr panel refinement time " << time_operations[amr_refine_time].count() << " seconds per refinement session" << endl;
//         cout << "- " << num_operations[amr_refine_time] << " amr sessions made at ";
//         cout << time_operations[amr_refine_time].count() / num_operations[amr_refine_time] << " seconds per refinement session" << endl;
//     }
//     cout << ".. Total interpolation time " << time_operations[interp_time].count() << " seconds" << endl;
//     cout << "- " << num_operations[interp_time] << " interpolations made at ";
//     cout << time_operations[interp_time].count() / num_operations[interp_time] << " seconds per interpolation" << endl;
//     cout << "... Total panel search time " << time_operations[search_time].count() << " seconds" << endl;
//     cout << "- " << num_operations[search_time] << " searches made at ";
//     cout << time_operations[search_time].count() / num_operations[search_time] << " seconds per search" << endl;
//     cout << "... Total interpolant evaluation time " << time_operations[eval_time].count() << " seconds" << endl;
//     cout << "- " << num_operations[eval_time] << " evaluations made at ";
//     cout <<  time_operations[eval_time].count() / num_operations[eval_time] << " seconds per evaluation " << endl;
//     cout <<". Total file writing time " << time_operations[file_time].count() << " seconds" << endl;
// }

