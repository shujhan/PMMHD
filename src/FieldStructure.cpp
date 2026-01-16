#include "FieldStructure.hpp"
#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <cfloat> // dbl_min
#include <cstddef>
using namespace std;
#if OPENACC_ENABLED
#include <accelmath.h>
#endif

Field::~Field() = default;


U_DirectSum::U_DirectSum() {}
U_DirectSum::U_DirectSum(double L, double epsilon) : L(L), epsilon(epsilon) {}
U_DirectSum::~U_DirectSum() = default;

void U_DirectSum::operator() (double* u1s, double* u2s, double* x_vals, int nx, 
                        double* y_vals, double* q_ws, int ny)
{    
    const double pi = std::atan(1.0) * 4.0;
    for (int i = 0; i < nx; i++) {
        for(int k = 0; k < ny; k++) {
            double denom = cosh(2* pi / L * (y_vals[i] - y_vals[k])) - cos(2* pi / L * (x_vals[i] - x_vals[k])) + epsilon * epsilon;
            u1s[i] -= 0.5/L * sinh(2 * pi / L * (y_vals[i] - y_vals[k])) / denom * q_ws[k];
            u2s[i] += 0.5/L * sin(2 * pi / L * (x_vals[i] - x_vals[k])) / denom * q_ws[k];
        }
    }
}

void U_DirectSum::print_field_obj() {
    cout << "-------------" << endl;
    cout << "Field object: " << endl;
}



// =========================
// treecode
// =========================
U_Treecode::U_Treecode() {}
U_Treecode::U_Treecode(double L_, double epsilon_,
                                   double mac_, int degree_, int max_source_,
                                   int max_target_, int verbosity_)
    : L(L_), epsilon(epsilon_), mac(mac_), degree(degree_),
      max_source(max_source_), max_target(max_target_), verbosity(verbosity_), 
      lambda(nullptr), particles_x(nullptr), particles_y(nullptr), velo_tc_reord_x(nullptr), velo_tc_reord_y(nullptr), 
      velo_tc_noreord_x(nullptr), velo_tc_noreord_y(nullptr), 
      tree_members{nullptr, nullptr}, leaf_members{nullptr, nullptr}, iList(nullptr), cList(nullptr)
      {
        #ifdef OPENACC_ENABLED
        #pragma acc enter data copyin(this)
        #endif
        P = degree_;
        PP = P + 1;
        Pflat = PP * PP;
}

U_Treecode::U_Treecode()
    {
    #ifdef OPENACC_ENABLED
        #pragma acc exit data delete(this)
    #endif
    };

void U_Treecode::print_field_obj() {
    std::cout << "[U_Treecode]\n";
    std::cout << "  L=" << L << " epsilon=" << epsilon << " beta=" << beta << "\n";
    std::cout << "  P=" << P << " PP=" << PP << " Pflat=" << Pflat << "\n";
    std::cout << "  N0=" << N0 << " theta=" << std::sqrt(sq_theta) << "\n";
    std::cout << "  delta=" << delta << "\n";
}

void U_Treecode::operator()(double* e1s, double* e2s,
                                 double* x_vals, int nx,
                                 double* y_vals, double* q_ws, int ny)
{
    numpars_s = (size_t)nx;
    // Make sure no false allocations
    cleanup();
    // Allocate internal arrays
    lambda      = new double[numpars_s];
    particles_x = new double[numpars_s];
    particles_y = new double[numpars_s];

    for (size_t i = 0; i < numpars_s; i++) {
        // particles_x[i] = x_vals[i];
        particles_x[i] = fmod(x_vals[i],L);
        particles_y[i] = y_vals[i];
        lambda[i]      = q_ws[i];
    }

#if OPENACC_ENABLED
std::cout << "Running with OpenACC" << std::endl;
#pragma acc enter data copyin(lambda[0:numpars_s])
#pragma acc enter data copyin(particles_x[0:numpars_s], particles_y[0:numpars_s])
#else
std::cout << "Running without OpenACC" << std::endl;
#endif


    //============  BLTC ============
    velo_tc_reord_x   = new double[numpars_s];
    velo_tc_reord_y   = new double[numpars_s];
    velo_tc_noreord_x = new double[numpars_s];
    velo_tc_noreord_y = new double[numpars_s];

#if OPENACC_ENABLED
#pragma acc enter data create(velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s])
#pragma acc enter data create(velo_tc_noreord_x[0:numpars_s], velo_tc_noreord_y[0:numpars_s])
#endif

    // Run BLTC treecode
    compute_RHS_BLTC();

    // copy result back to output arrays
    for (size_t i = 0; i < numpars_s; i++) {
        e1s[i] = velo_tc_noreord_x[i];
        e2s[i] = velo_tc_noreord_y[i];
    }

#if OPENACC_ENABLED
#pragma acc exit data delete(lambda[0:numpars_s])
#pragma acc exit data delete(particles_x[0:numpars_s], particles_y[0:numpars_s])
#pragma acc exit data delete(velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s])
#pragma acc exit data delete(velo_tc_noreord_x[0:numpars_s], velo_tc_noreord_y[0:numpars_s])
#endif

    cleanup();
}


void U_Treecode::cleanup() {
    // interaction + cluster lists (host frees, OpenACC frees inside functions)
    if (cList) free_cluster_list();
    if (iList) free_interaction_list();

    tree.clear();
    leaf.clear();

    node_count = 0;
    leaf_count = 0;

    if (tree_members[0]) { delete[] tree_members[0]; tree_members[0] = nullptr; }
    if (tree_members[1]) { delete[] tree_members[1]; tree_members[1] = nullptr; }
    if (leaf_members[0]) { delete[] leaf_members[0]; leaf_members[0] = nullptr; }
    if (leaf_members[1]) { delete[] leaf_members[1]; leaf_members[1] = nullptr; }

    if (lambda)      { delete[] lambda;      lambda = nullptr; }
    if (particles_x) { delete[] particles_x; particles_x = nullptr; }
    if (particles_y) { delete[] particles_y; particles_y = nullptr; }

    if (velo_tc_reord_x)   { delete[] velo_tc_reord_x;   velo_tc_reord_x = nullptr; }
    if (velo_tc_reord_y)   { delete[] velo_tc_reord_y;   velo_tc_reord_y = nullptr; }
    if (velo_tc_noreord_x) { delete[] velo_tc_noreord_x; velo_tc_noreord_x = nullptr; }
    if (velo_tc_noreord_y) { delete[] velo_tc_noreord_y; velo_tc_noreord_y = nullptr; }
}


// ===============================
// compute_RHS_BLTC
// ===============================
void U_Treecode::compute_RHS_BLTC()
{
#if OPENACC_ENABLED
#pragma acc update host(particles_x[0:numpars_s], particles_y[0:numpars_s])
#endif



}