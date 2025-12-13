#include "AMRStructure.hpp"

int evaluate_u_field(std::vector<double>& u1s_local, std::vector<double>& u2s_local, 
                     std::vector<double>& xs_local,std::vector<double>& ys_local,
                     std::vector<double>& ws_local, double t) {
    std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
    (*calculate_e)(u1s_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
                    ytemp_cpy.data(), ws_local_local.data(), ytemp_cpy.size());
    return 0;
}


int evaluate_b_field(std::vector<double>& b1s_local, std::vector<double>& b2s_local, 
                    std::vector<double>& xs_local,std::vector<double>& ys_local,
                        std::vector<double>& ws_local, double t) {
    std::vector<double> xtemp_cpy(xs_local), ytemp_cpy(ys_local);
    (*calculate_e)(u1s_local.data(), xtemp_cpy.data(), xtemp_cpy.size(),
                    ytemp_cpy.data(), ws_local_local.data(), ytemp_cpy.size());
    return 0;
}