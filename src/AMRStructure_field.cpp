#include "AMRStructure.hpp"

int AMRStructure::evaluate_u_field(std::vector<double>& u1s_local, std::vector<double>& u2s_local, 
                     std::vector<double>& xs_local,std::vector<double>& ys_local,
                     std::vector<double>& ws_local, double t) {
    std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
    std::vector<double> ws_temp_cpy(ws_local);
    KernelMode m = original;
    calculate_e->set_mode(m);
    (*calculate_e)(u1s_local.data(), u2s_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
                    ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
    return 0;
}



int AMRStructure::evaluate_b_field(std::vector<double>& b1s_local, std::vector<double>& b2s_local, 
                    std::vector<double>& xs_local,std::vector<double>& ys_local,
                        std::vector<double>& ws_local, double t) {
    std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
    std::vector<double> ws_temp_cpy(ws_local);
    KernelMode m = original;
    calculate_e->set_mode(m);
    (*calculate_e)(b1s_local.data(), b2s_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
                    ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
    return 0;
}


int AMRStructure::evaluate_u1s_grad(std::vector<double>& u1s_grad_x_local, std::vector<double>& u1s_grad_y_local, 
                    std::vector<double>& xs_local,std::vector<double>& ys_local,
                        std::vector<double>& ws_local, double t) {
    std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
    std::vector<double> ws_temp_cpy(ws_local);
    KernelMode m = u1_grad;
    calculate_e->set_mode(m);
    (*calculate_e)(u1s_grad_x_local.data(), u1s_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
                    ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
    return 0;
}

int AMRStructure::evaluate_u2s_grad(std::vector<double>& u2s_grad_x_local, std::vector<double>& u2s_grad_y_local, 
                    std::vector<double>& xs_local,std::vector<double>& ys_local,
                        std::vector<double>& ws_local, double t) {
    std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
    std::vector<double> ws_temp_cpy(ws_local);
    KernelMode m = u2_grad;
    calculate_e->set_mode(m);
    (*calculate_e)(u2s_grad_x_local.data(), u2s_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
                    ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
    return 0;
}

int AMRStructure::evaluate_b1s_grad(std::vector<double>& b1s_grad_x_local, std::vector<double>& b1s_grad_y_local, 
                    std::vector<double>& xs_local,std::vector<double>& ys_local,
                        std::vector<double>& ws_local, double t) {
    std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
    std::vector<double> ws_temp_cpy(ws_local);
    KernelMode m = u1_grad;
    calculate_e->set_mode(m);
    (*calculate_e)(b1s_grad_x_local.data(), b1s_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
                    ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
    return 0;
}

int AMRStructure::evaluate_b2s_grad(std::vector<double>& b2s_grad_x_local, std::vector<double>& b2s_grad_y_local, 
                    std::vector<double>& xs_local,std::vector<double>& ys_local,
                        std::vector<double>& ws_local, double t) {
    std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
    std::vector<double> ws_temp_cpy(ws_local);
    KernelMode m = u2_grad;
    calculate_e->set_mode(m);
    (*calculate_e)(b2s_grad_x_local.data(), b2s_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
                    ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
    return 0;
}


int AMRStructure::evaluate_vorticity_grad(std::vector<double>& vorticity_grad_x_local, std::vector<double>& vorticity_grad_y_local, 
                    std::vector<double>& xs_local,std::vector<double>& ys_local,
                        std::vector<double>& ws_local, double t) {
    std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
    std::vector<double> ws_temp_cpy(ws_local);
    KernelMode m = vorticity_grad;
    calculate_e->set_mode(m);
    (*calculate_e)(vorticity_grad_x_local.data(), vorticity_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
                    ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
    return 0;
}


int AMRStructure::evaluate_j_grad(std::vector<double>& j_grad_x_local, std::vector<double>& j_grad_y_local, 
                    std::vector<double>& xs_local,std::vector<double>& ys_local,
                        std::vector<double>& ws_local, double t) {
    std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
    std::vector<double> ws_temp_cpy(ws_local);
    KernelMode m = vorticity_grad;
    calculate_e->set_mode(m);
    (*calculate_e)(j_grad_x_local.data(), j_grad_y_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
                    ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
    return 0;
}


int AMRStructure::evaluate_vorticity_laplacian(std::vector<double>& vorticity_laplacian_local, std::vector<double>& vorticity_none_local, 
                    std::vector<double>& xs_local,std::vector<double>& ys_local,
                        std::vector<double>& ws_local, double t) {
    std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
    std::vector<double> ws_temp_cpy(ws_local);
    KernelMode m = laplacian;
    calculate_e->set_mode(m);
    (*calculate_e)(vorticity_laplacian_local.data(), vorticity_none_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
                    ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
    return 0;
}

int AMRStructure::evaluate_j_laplacian(std::vector<double>& j_laplacian_local, std::vector<double>& j_none_local, 
                    std::vector<double>& xs_local,std::vector<double>& ys_local,
                        std::vector<double>& ws_local, double t) {
    std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
    std::vector<double> ws_temp_cpy(ws_local);
    KernelMode m = laplacian;
    calculate_e->set_mode(m);
    (*calculate_e)(j_laplacian_local.data(),j_none_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
                    ytemp_cpy.data(), ws_temp_cpy.data(), ytemp_cpy.size());
    return 0;
}

