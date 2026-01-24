#include "FieldStructure.hpp"
#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <cfloat> // dbl_min
#include <cstddef>
#include <sys/times.h>
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
    switch (mode)
    {
        case original:
            // cout << "orginal kernel: " <<endl;
            #ifdef OPENACC_ENABLED
            #pragma acc parallel loop independent
            #else
            #pragma omp parallel for
            #endif
                for (int i = 0; i < nx; i++) {
                    double u1 = 0.0;
                    double u2 = 0.0;
                #ifdef OPENACC_ENABLED
                #pragma acc loop independent reduction(+:u1, u2)
                #endif
                    for(int k = 0; k < ny; k++) {
                        double x_diff = 2*pi/ L * (x_vals[i] - x_vals[k]);
                        double y_diff = 2*pi/ L * (y_vals[i] - y_vals[k]);
                        double denom = cosh(y_diff) - cos(x_diff) + epsilon * epsilon;
                        u1 -= 0.5/L * sinh(y_diff) / denom * q_ws[k];
                        u2 += 0.5/L * sin(x_diff) / denom * q_ws[k];
                    }
                    u1s[i] = u1;
                    u2s[i] = u2;
                }
            break;

        case u1_grad: // for u1s_grad_x, u1s_grad_y, b1s_grad_x, b1s_grad_y
            // cout << "u1_grad kernel: " <<endl;
            #ifdef OPENACC_ENABLED
            #pragma acc parallel loop independent
            #else
            #pragma omp parallel for
            #endif
                for (int i = 0; i < nx; i++) {
                    double u1 = 0.0;
                    double u2 = 0.0;
                #ifdef OPENACC_ENABLED
                #pragma acc loop independent reduction(+:u1, u2)
                #endif
                    for(int k = 0; k < ny; k++) {
                        double x_diff = 2*pi/ L * (x_vals[i] - x_vals[k]);
                        double y_diff = 2*pi/ L * (y_vals[i] - y_vals[k]);
                        double denom_sqr = (cosh(y_diff) - cos(x_diff) + epsilon * epsilon) * (cosh(y_diff) - cos(x_diff) + epsilon * epsilon);
                        double constant_c = pi/L/L;
                        u1 += constant_c * sinh(y_diff) * sin(x_diff) / denom_sqr * q_ws[k];
                        u2 += constant_c * (cos(x_diff) * cosh(y_diff) - 1 - epsilon * epsilon * cosh(y_diff)) / denom_sqr * q_ws[k];
                    }
                    u1s[i] = u1;
                    u2s[i] = u2;
                }
            break;

        case u2_grad: // for u2s_grad_x, u2s_grad_y, b2s_grad_x, b2s_grad_y
            // cout << "u2_grad kernel: " <<endl;
            #ifdef OPENACC_ENABLED
            #pragma acc parallel loop independent
            #else
            #pragma omp parallel for
            #endif
                for (int i = 0; i < nx; i++) {
                    double u1 = 0.0;
                    double u2 = 0.0;
                #ifdef OPENACC_ENABLED
                #pragma acc loop independent reduction(+:u1, u2)
                #endif
                    for(int k = 0; k < ny; k++) {
                        double x_diff = 2*pi/ L * (x_vals[i] - x_vals[k]);
                        double y_diff = 2*pi/ L * (y_vals[i] - y_vals[k]);
                        double denom_sqr = (cosh(y_diff) - cos(x_diff) + epsilon * epsilon) * (cosh(y_diff) - cos(x_diff) + epsilon * epsilon);
                        double constant_c = pi/L/L;
                        u1 += constant_c * (cos(x_diff) * cosh(y_diff) - 1 + epsilon * epsilon * cos(x_diff)) / denom_sqr * q_ws[k];
                        u2 += -1 * constant_c * sin(x_diff) * sinh(y_diff) / denom_sqr * q_ws[k];
                    }
                    u1s[i] = u1;
                    u2s[i] = u2;
                }
            break;

        case vorticity_grad: // for vorticity_grad_x, vorticity_grad_x, j_grad_x, j_grad_y
            // cout << "vorticity_grad kernel: " <<endl;
            #ifdef OPENACC_ENABLED
            #pragma acc parallel loop independent
            #else
            #pragma omp parallel for
            #endif
                for (int i = 0; i < nx; i++) {
                    double u1 = 0.0;
                    double u2 = 0.0;
                #ifdef OPENACC_ENABLED
                #pragma acc loop independent reduction(+:u1, u2)
                #endif
                    for(int k = 0; k < ny; k++) {
                        double x_diff = 2*pi/ L * (x_vals[i] - x_vals[k]);
                        double y_diff = 2*pi/ L * (y_vals[i] - y_vals[k]);
                        double denom_cube = (cosh(y_diff) - cos(x_diff) + epsilon * epsilon) * (cosh(y_diff) - cos(x_diff) + epsilon * epsilon) * (cosh(y_diff) - cos(x_diff) + epsilon * epsilon);
                        double constant_c = 2 * pi * pi * epsilon * epsilon /L/L/L; 
                        u1 += constant_c * sin(x_diff) *(-3*cosh(y_diff) - cos(x_diff) - epsilon * epsilon) / denom_cube * q_ws[k];
                        u2 += constant_c * sinh(y_diff) *(-cosh(y_diff) - 3 * cos(x_diff) + epsilon * epsilon) / denom_cube * q_ws[k];
                    }
                    u1s[i] = u1;
                    u2s[i] = u2;
                }
            break;

        case laplacian: // for vorticity_laplacian, j_laplacian;
            // cout << "laplacian kernel: " <<endl;
            #ifdef OPENACC_ENABLED
            #pragma acc parallel loop independent
            #else
            #pragma omp parallel for
            #endif
                for (int i = 0; i < nx; i++) {
                    double u1 = 0.0;
                    double u2 = 0.0;
                #ifdef OPENACC_ENABLED
                #pragma acc loop independent reduction(+:u1, u2)
                #endif
                    for(int k = 0; k < ny; k++) {
                        double x_diff = 2*pi/ L * (x_vals[i] - x_vals[k]);
                        double y_diff = 2*pi/ L * (y_vals[i] - y_vals[k]);
                        double C = cos(x_diff);
                        double S = sin(x_diff);
                        double H = cosh(y_diff);
                        double Sh = sinh(y_diff);
                        double eps_sqr = epsilon * epsilon;
                        double eps_4 = epsilon * epsilon * epsilon * epsilon;
                        double denom_4 = (H - C + eps_sqr)*(H - C + eps_sqr)*(H - C + eps_sqr)*(H - C + eps_sqr);
                        double constant_c = 4 * pi * pi * pi * epsilon * epsilon /L/L/L/L; 
                        u1 += constant_c * (C*C*C +5*C*C*H -5*C*H*H -8*C*H*eps_sqr +2*C*S*S + 10*C*Sh*Sh -C*eps_4
                                -H*H*H + 10*H*S*S +2*H*Sh*Sh + H*eps_4 + 4*S*S*eps_sqr - 4*Sh*Sh*eps_sqr) 
                                / denom_4 * q_ws[k];               
                        u2 += 0;
                    }
                    u1s[i] = u1;
                    u2s[i] = u2;
                }
            break;

        default:
            // cout << "orginal kernel: " <<endl;
            #ifdef OPENACC_ENABLED
            #pragma acc parallel loop independent
            #else
            #pragma omp parallel for
            #endif
                for (int i = 0; i < nx; i++) {
                    double u1 = 0.0;
                    double u2 = 0.0;
                #ifdef OPENACC_ENABLED
                #pragma acc loop independent reduction(+:u1, u2)
                #endif
                    for(int k = 0; k < ny; k++) {
                        // double denom = cosh(2* pi / L * (y_vals[i] - y_vals[k])) - cos(2* pi / L * (x_vals[i] - x_vals[k])) + epsilon * epsilon;
                        // u1 -= 0.5/L * sinh(2 * pi / L * (y_vals[i] - y_vals[k])) / denom * q_ws[k];
                        // u2 += 0.5/L * sin(2 * pi / L * (x_vals[i] - x_vals[k])) / denom * q_ws[k];
                        double x_diff = 2*pi/ L * (x_vals[i] - x_vals[k]);
                        double y_diff = 2*pi/ L * (y_vals[i] - y_vals[k]);
                        double denom = cosh(y_diff) - cos(x_diff) + epsilon * epsilon;
                        u1 -= 0.5/L * sinh(y_diff) / denom * q_ws[k];
                        u2 += 0.5/L * sin(x_diff) / denom * q_ws[k];
                    }
                    u1s[i] = u1;
                    u2s[i] = u2;
                }
            break;
    }

}

void U_DirectSum::print_field_obj() {
    cout << "-------------" << endl;
    cout << "Field object: " << endl;
}


void U_DirectSum::set_mode(KernelMode m){
    mode = m;
}


// =========================
// treecode
// =========================
U_Treecode::U_Treecode() {}
U_Treecode::U_Treecode(double L_, double epsilon_,
                                   double mac_, int degree_, int max_source_,
                                   int max_target_)
    : L(L_), epsilon(epsilon_), mac(mac_), degree(degree_),
      max_source(max_source_), max_target(max_target_), 
      lambda(nullptr), particles_x(nullptr), particles_y(nullptr), velo_tc_reord_x(nullptr), velo_tc_reord_y(nullptr), 
      velo_tc_noreord_x(nullptr), velo_tc_noreord_y(nullptr), 
      tree_members{nullptr, nullptr}, leaf_members{nullptr, nullptr}, iList(nullptr), cList(nullptr)
      {
        #ifdef OPENACC_ENABLED
        #pragma acc enter data copyin(this)
        #endif
        P = degree;
        PP = P + 1;
        Pflat = PP * PP;
        N0 = max_source;
        sq_theta = mac * mac;
}

U_Treecode::~U_Treecode()
    {
    #ifdef OPENACC_ENABLED
        #pragma acc exit data delete(this)
    #endif
    };

void U_Treecode::print_field_obj() {
    std::cout << "[U_Treecode]\n";
    std::cout << "  L=" << L << " epsilon=" << epsilon << endl;
    std::cout << "  P=" << P << " PP=" << PP << " Pflat=" << Pflat << "\n";
    std::cout << "  N0=" << N0 << " theta=" << sqrt(sq_theta) << "\n";
    // std::cout << "  delta=" << delta << "\n";
}

void U_Treecode::set_mode(KernelMode m){
    mode = m;
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

    xyminmax[0] = minval(particles_x, numpars_s);
    xyminmax[1] = maxval(particles_x, numpars_s);
    xyminmax[2] = minval(particles_y, numpars_s);
    xyminmax[3] = maxval(particles_y, numpars_s);

    // reset output buffers
    for (size_t i = 0; i < numpars_s; i++) {
        velo_tc_reord_x[i]   = 0.0;
        velo_tc_reord_y[i]   = 0.0;
        velo_tc_noreord_x[i] = 0.0;
        velo_tc_noreord_y[i] = 0.0;
    }

#if OPENACC_ENABLED
#pragma acc update device(velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s])
#pragma acc update device(velo_tc_noreord_x[0:numpars_s], velo_tc_noreord_y[0:numpars_s])
#pragma acc update device(lambda[0:numpars_s])
#endif


    // temp indices
    int* pt_temp_index = new int[numpars_s];
    int* pt_temp_old_index = new int[numpars_s];
    for (size_t i = 0; i < numpars_s; i++) {
        pt_temp_index[i] = -1;
        pt_temp_old_index[i] = (int)i;
    }


    // tree
    build_tree_init();
    build_tree_2D_Recursive(0, 0, particles_x, particles_y, pt_temp_index, pt_temp_old_index);

#if OPENACC_ENABLED
#pragma acc update device(particles_x[0:numpars_s], particles_y[0:numpars_s])
#pragma acc update device(lambda[0:numpars_s])
#endif

    assert(node_count == tree.size());
    cout << "Tree size = " << node_count << endl;

    assert(leaf_count == leaf.size());
    cout << "Leaf size = " << leaf_count << endl;

    //====================== Compute interaction list =================

    // interaction lists
    std::vector<std::vector<size_t>> Interaction_List_far(leaf_count);
    std::vector<std::vector<size_t>> Interaction_List_near(leaf_count);

    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        build_interaction_list(leaf_index, 0, Interaction_List_far, Interaction_List_near);
    }

    //=============== Copy interaction list and cluster list to the device ===============
    alloc_set_interaction_list(Interaction_List_far, Interaction_List_near);
    alloc_set_cluster_list(particles_x, particles_y);

    // tree_members + leaf_members
    tree_members[0] = new size_t[node_count];
    tree_members[1] = new size_t[node_count];

    for (size_t panel_index = 0; panel_index < node_count; panel_index++) {
        tree_members[0][panel_index] = tree[panel_index].members[0];
        tree_members[1][panel_index] = tree[panel_index].members[1];
    }

    leaf_members[0] = new size_t[leaf_count];
    leaf_members[1] = new size_t[leaf_count];

    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        leaf_members[0][leaf_index] = tree[leaf[leaf_index]].members[0];
        leaf_members[1][leaf_index] = tree[leaf[leaf_index]].members[1];
    }

#if OPENACC_ENABLED
#pragma acc enter data copyin(tree_members[0:2][0:node_count])
#pragma acc enter data copyin(leaf_members[0:2][0:leaf_count])
#endif

    // compute
    Compute_SUM();

#if OPENACC_ENABLED
#pragma acc update host(velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s])
#endif

//========= Change back to original order =====================
    double* lambda_temp = new double[numpars_s];
    for (size_t i = 0; i < numpars_s; i++) {
        velo_tc_noreord_x[pt_temp_old_index[i]] = velo_tc_reord_x[i];
        velo_tc_noreord_y[pt_temp_old_index[i]] = velo_tc_reord_y[i];
        lambda_temp[i] = lambda[i];
    }

    for (size_t i = 0; i < numpars_s; i++) {
        lambda[pt_temp_old_index[i]] = lambda_temp[i];
    }


#if OPENACC_ENABLED
#pragma acc update device(lambda[0:numpars_s])
#pragma acc update device(velo_tc_noreord_x[0:numpars_s], velo_tc_noreord_y[0:numpars_s])

#pragma acc exit data delete(tree_members[0:2][0:node_count])
#pragma acc exit data delete(leaf_members[0:2][0:leaf_count])
#endif

    delete[] lambda_temp;

    // free host members arrays
    delete[] tree_members[0]; tree_members[0] = nullptr;
    delete[] tree_members[1]; tree_members[1] = nullptr;
    delete[] leaf_members[0]; leaf_members[0] = nullptr;
    delete[] leaf_members[1]; leaf_members[1] = nullptr;

    // clear vectors
    tree.clear();
    leaf.clear();

    // free lists
    free_cluster_list();
    free_interaction_list();

    leaf_count = 0;
    node_count = 0;

    delete[] pt_temp_old_index;
    delete[] pt_temp_index;

}


// ===============================
// build tree init
// ===============================
void U_Treecode::build_tree_init() {
    panel temp_panel;

    temp_panel.members[0] = 0;
    temp_panel.members[1] = numpars_s - 1;

    temp_panel.xinterval[0] = xyminmax[0];
    temp_panel.xinterval[1] = xyminmax[1];
    temp_panel.yinterval[0] = xyminmax[2];
    temp_panel.yinterval[1] = xyminmax[3];

    temp_panel.xc = 0.5 * (temp_panel.xinterval[0] + temp_panel.xinterval[1]);
    temp_panel.yc = 0.5 * (temp_panel.yinterval[0] + temp_panel.yinterval[1]);

    double xL = temp_panel.xinterval[1] - temp_panel.xinterval[0];
    double yL = temp_panel.yinterval[1] - temp_panel.yinterval[0];

    double sq_r = 0.25 * (xL * xL + yL * yL);
    temp_panel.MAC = sq_r / sq_theta; // MAC = r^2 / mac^2

    tree.push_back(temp_panel);
    node_count = 1;
}


void U_Treecode::build_tree_2D_Recursive(size_t panel_index, int level,
                                              double* pt_x, double* pt_y,
                                              int* pt_index, int* pt_old_index) 
{
    size_t n = tree[panel_index].members[1] - tree[panel_index].members[0] + 1;
    if (n >= (size_t)N0) {
        split_tree_node(panel_index, pt_x, pt_y, pt_index, pt_old_index);

        for (size_t i = 0; i < tree[panel_index].children.size(); i++) {
            size_t panel_index_new = tree[panel_index].children[i];
            build_tree_2D_Recursive(panel_index_new, level + 1, pt_x, pt_y, pt_index, pt_old_index);
        }
    } else {
        leaf.push_back(panel_index);
        leaf_count++;
    }
}

void U_Treecode::split_tree_node(size_t panel_index,
                                      double* pt_x, double* pt_y,
                                      int* pt_index, int* pt_old_index) {
    double xL = tree[panel_index].xinterval[1] - tree[panel_index].xinterval[0];
    double yL = tree[panel_index].yinterval[1] - tree[panel_index].yinterval[0];

    double L_max = xL;
    if (yL > L_max) L_max = yL;

    int XY_flag = 0;
    const double ratio = std::sqrt(2.0);

    if (xL * ratio > L_max) XY_flag += 2;
    if (yL * ratio > L_max) XY_flag += 1;

    switch (XY_flag) {
        case 1:
        case 2:
            split_2(panel_index, XY_flag, pt_x, pt_y, pt_index, pt_old_index);
            break;
        case 3:
            split_4(panel_index, pt_x, pt_y, pt_index, pt_old_index);
            break;
        default:
            break;
    }
}


void U_Treecode::split_2(size_t panel_index, int split_code,
                              double* pt_x, double* pt_y,
                              int* pt_index, int* pt_old_index) {
    panel child[2];

    /*
     -------------------
     |        |        |
     |        |        |
     |Child 0 |Child 1 |
     |        |        |
     |        |        |
     -----------------------------> axis A
   start     mid      end
     */

    double tp_x0 = tree[panel_index].xinterval[0];
    double tp_x1 = tree[panel_index].xinterval[1];
    double tp_y0 = tree[panel_index].yinterval[0];
    double tp_y1 = tree[panel_index].yinterval[1];

    for (int i = 0; i < 2; i++) {
        child[i].xinterval[0] = tp_x0;
        child[i].xinterval[1] = tp_x1;
        child[i].yinterval[0] = tp_y0;
        child[i].yinterval[1] = tp_y1;
    }

    double xL = tp_x1 - tp_x0;
    double yL = tp_y1 - tp_y0;

    double* intervalA[2] = {NULL, NULL};
    double* coordA = NULL;
    double startpointA = 0.0, midpointA = 0.0, endpointA = 0.0;

    if (split_code == 2) { // XY = 10, A is X
        xL *= 0.5;
        intervalA[0] = child[0].xinterval;
        intervalA[1] = child[1].xinterval;
        coordA = pt_x;
        startpointA = tp_x0;
        endpointA = tp_x1;
    }
    else if (split_code == 1) { // XY = 01, A is Y
        yL *= 0.5;
        intervalA[0] = child[0].yinterval;
        intervalA[1] = child[1].yinterval;
        coordA = pt_y;
        startpointA = tp_y0;
        endpointA = tp_y1;
    }

    midpointA = 0.5 * (startpointA + endpointA);

    // Child 0 ends with mid point on axis A
    intervalA[0][1] = midpointA;

    // Child 1 begins with mid point on axis A
    intervalA[1][0] = midpointA;

    double sq_r = 0.25 * (xL * xL + yL * yL); // r^2
    double MAC = sq_r / sq_theta; // MAC = r^2 / theta^2

    for (int i = 0; i < 2; i++) {
        child[i].xc = 0.5 * (child[i].xinterval[0] + child[i].xinterval[1]);
        child[i].yc = 0.5 * (child[i].yinterval[0] + child[i].yinterval[1]);
        child[i].MAC = MAC;
    }

    vector<size_t> v[2];
    size_t start = tree[panel_index].members[0];
    size_t end = tree[panel_index].members[1];
    size_t* addr_table = new size_t[end - start + 1];

    size_t index;
    for (index = start; index <= end; index++) {
        pt_index[index] = index;
        addr_table[index - start] = index;

        if (coordA[index] <= midpointA)
            v[0].push_back(index);
        else
            v[1].push_back(index);
    }

    size_t seq = start;
    for (size_t j = 0; j < 2; j++) {
        size_t size = v[j].size();

        if (size >= 1) {
            for (size_t k = 0; k < size; k++) {
                if (k == 0)
                    child[j].members[0] = seq;
                if (k == size - 1)
                    child[j].members[1] = seq;

                index = v[j][k];

                // This uses an address table
                size_t pos = addr_table[index - start];
                size_t out = pt_index[seq];
                Swap(pos, seq, pt_x, pt_y, pt_index, pt_old_index);
                addr_table[index - start] = seq;
                addr_table[out - start] = pos;

                seq++;
            }

            node_count++;
            tree[panel_index].children.push_back(node_count - 1);
            tree.push_back(child[j]);
            v[j].clear();
        }
    }

    delete[] addr_table;
}

void U_Treecode::split_4(size_t panel_index,
                              double* pt_x, double* pt_y,
                              int* pt_index, int* pt_old_index) {
    panel child[4];

    /*
      ^ axis y
      |
  end -------------------
      |        |        |
      |Child 2 |Child 3 |
  mid |------------------
      |        |        |
start |Child 0 |Child 1 |
      -----------------------------> axis x
    start     mid      end
     */

    double tp_x0 = tree[panel_index].xinterval[0];
    double tp_x1 = tree[panel_index].xinterval[1];
    double tp_y0 = tree[panel_index].yinterval[0];
    double tp_y1 = tree[panel_index].yinterval[1];

    double xL = 0.5 * (tp_x1 - tp_x0);
    double yL = 0.5 * (tp_y1 - tp_y0);

    double midpointx = 0.5 * (tp_x0 + tp_x1);
    double midpointy = 0.5 * (tp_y0 + tp_y1);

    child[0].xinterval[0] = tp_x0;
    child[0].xinterval[1] = midpointx;
    child[0].yinterval[0] = tp_y0;
    child[0].yinterval[1] = midpointy;

    child[1].xinterval[0] = midpointx;
    child[1].xinterval[1] = tp_x1;
    child[1].yinterval[0] = tp_y0;
    child[1].yinterval[1] = midpointy;

    child[2].xinterval[0] = tp_x0;
    child[2].xinterval[1] = midpointx;
    child[2].yinterval[0] = midpointy;
    child[2].yinterval[1] = tp_y1;

    child[3].xinterval[0] = midpointx;
    child[3].xinterval[1] = tp_x1;
    child[3].yinterval[0] = midpointy;
    child[3].yinterval[1] = tp_y1;

    double sq_r = 0.25 * (xL * xL + yL * yL); // r^2
    double MAC = sq_r / sq_theta; // MAC = r^2 / theta^2

    for (int i = 0; i < 4; i++) {
        child[i].xc = 0.5 * (child[i].xinterval[0] + child[i].xinterval[1]);
        child[i].yc = 0.5 * (child[i].yinterval[0] + child[i].yinterval[1]);
        child[i].MAC = MAC;
    }

    vector<size_t> v[8];
    size_t start = tree[panel_index].members[0];
    size_t end = tree[panel_index].members[1];
    size_t* addr_table = new size_t[end - start + 1];

    size_t index;
    for (index = start; index <= end; index++) {
        pt_index[index] = index;
        addr_table[index - start] = index;

        if (pt_x[index] <= midpointx && pt_y[index] <= midpointy)
            v[0].push_back(index);
        else if (pt_x[index] > midpointx && pt_y[index] <= midpointy)
            v[1].push_back(index);
        else if (pt_x[index] <= midpointx && pt_y[index] > midpointy)
            v[2].push_back(index);
        else if (pt_x[index] > midpointx && pt_y[index] > midpointy)
            v[3].push_back(index);
    }

    size_t seq = start;
    for (size_t j = 0; j < 4; j++) {
        size_t size = v[j].size();

        if (size >= 1) {
            for (size_t k = 0; k < size; k++) {
                if (k == 0)
                    child[j].members[0] = seq;
                if (k == size - 1)
                    child[j].members[1] = seq;

                index = v[j][k];
                // This uses an address table
                size_t pos = addr_table[index - start];
                size_t out = pt_index[seq];
                Swap(pos, seq, pt_x, pt_y, pt_index, pt_old_index);
                addr_table[index - start] = seq;
                addr_table[out - start] = pos;

                seq++;
            }

            node_count++;
            tree[panel_index].children.push_back(node_count - 1);
            tree.push_back(child[j]);
            v[j].clear();
        }
    }
    delete[] addr_table;
}

void U_Treecode::Swap(size_t i, size_t j,
                           double* pt_x, double* pt_y,
                           int* pt_index, int* pt_old_index) {
    if (i == j)
        return;

    double x = pt_x[i];
    double y = pt_y[i];
    size_t index = pt_index[i];
    size_t old_index = pt_old_index[i];

    double temp_lambda;
    temp_lambda = lambda[i];

    pt_x[i] = pt_x[j];
    pt_y[i] = pt_y[j];
    pt_index[i] = pt_index[j];
    pt_old_index[i] = pt_old_index[j];

    lambda[i] = lambda[j];

    pt_x[j] = x;
    pt_y[j] = y;
    pt_index[j] = index;
    pt_old_index[j] = old_index;

    lambda[j] = temp_lambda;
}


void U_Treecode::build_interaction_list(
    size_t leaf_index, size_t panel_index,
    std::vector<std::vector<size_t>>& Interaction_List_far,
    std::vector<std::vector<size_t>>& Interaction_List_near)
{
    double xcb = tree[leaf[leaf_index]].xc;
    double ycb = tree[leaf[leaf_index]].yc;

    double xcp = tree[panel_index].xc;
    double ycp = tree[panel_index].yc;

    double x_vb = tree[leaf[leaf_index]].xinterval[1];
    double y_vb = tree[leaf[leaf_index]].yinterval[1];

    double x_vp = tree[panel_index].xinterval[1];
    double y_vp = tree[panel_index].yinterval[1];

    // double dfbpx = xcb - xcp; 
    // periodic in x 
    double r_dist = fabs(xcb - xcp);
    double dfbpx = fmin(r_dist, L - r_dist);

    double dfbpy = ycb - ycp;

    double dfbx = x_vb - xcb;
    double dfby = y_vb - ycb;

    double dfpx = x_vp - xcp;
    double dfpy = y_vp - ycp;

    double sq_R  = dfbpx * dfbpx + dfbpy * dfbpy;
    double sq_rb = dfbx * dfbx + dfby * dfby;
    double sq_rp = dfpx * dfpx + dfpy * dfpy;

    if (tree[panel_index].children.size() == 0) // A leaf panel
        Interaction_List_near[leaf_index].push_back(panel_index);
    else {
        // std::cout << "leaf=" << leaf_index
        //   << " panel=" << panel_index
        //   << " R=" << sqrt(sq_R)
        //   << " rb=" << sqrt(sq_rb)
        //   << " rp=" << sqrt(sq_rp)
        //   << " ratio=" << (sqrt(sq_rb)+sqrt(sq_rp))/sqrt(sq_R)
        //   << "\n";

        if ((sqrt(sq_rb) + sqrt(sq_rp)) / sqrt(sq_R) < sqrt(sq_theta))
            Interaction_List_far[leaf_index].push_back(panel_index);
        else {
            size_t size = tree[panel_index].children.size();
            if (size == 0) // It is a leaf
                Interaction_List_near[leaf_index].push_back(panel_index);
            else {
                for (size_t i_children = 0; i_children < size; i_children++)
                    build_interaction_list(leaf_index, tree[panel_index].children[i_children], Interaction_List_far, Interaction_List_near);
            }
        }
    }
}


void U_Treecode::alloc_set_interaction_list(
    const std::vector<std::vector<size_t>>& Interaction_List_far,
    const std::vector<std::vector<size_t>>& Interaction_List_near)
{
    iList = new interaction_list[leaf_count];

#if OPENACC_ENABLED
#pragma acc enter data create(iList[0:leaf_count])
#endif

    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t far_size = Interaction_List_far[leaf_index].size();
        iList[leaf_index].far_list_size = far_size;
        iList[leaf_index].far_list = (size_t*) Interaction_List_far[leaf_index].data();

        size_t near_size = Interaction_List_near[leaf_index].size();
        iList[leaf_index].near_list_size = near_size;
        iList[leaf_index].near_list = (size_t*) Interaction_List_near[leaf_index].data();

#if OPENACC_ENABLED
#pragma acc enter data create(iList[leaf_index].far_list[0:far_size], iList[leaf_index].near_list[0:near_size])
#pragma acc update device(iList[leaf_index].far_list_size, iList[leaf_index].far_list[0:far_size], iList[leaf_index].near_list_size, iList[leaf_index].near_list[0:near_size])
#endif
    }
}


void U_Treecode::free_interaction_list()
{
#if OPENACC_ENABLED
    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t far_size = iList[leaf_index].far_list_size;
        size_t near_size = iList[leaf_index].near_list_size;
#pragma acc exit data delete(iList[leaf_index].far_list[0:far_size], iList[leaf_index].near_list[0:near_size])
    }
#endif

#if OPENACC_ENABLED
#pragma acc exit data delete(iList[0:leaf_count])
#endif

    delete[] iList;
    iList = nullptr;
}


void U_Treecode::alloc_set_cluster_list(double* pt_x, double* pt_y)
{
    cList = new cluster_list[node_count];

#if OPENACC_ENABLED
#pragma acc enter data create(cList[0:node_count])
#endif

    for (size_t tree_index = 0; tree_index < node_count; tree_index++) {
        cList[tree_index].t1 = new double[PP];
        cList[tree_index].t2 = new double[PP];

        cList[tree_index].moments = new double[Pflat];
        memset(cList[tree_index].moments, 0, sizeof(double) * Pflat);
    
#if OPENACC_ENABLED
#pragma acc enter data create(cList[tree_index].t1[0:PP], cList[tree_index].t2[0:PP], \
cList[tree_index].moments[0:Pflat])
#endif
    }

    double h = pi / P;
    double t[PP] = {0.0};
    for (int i = 0; i < PP; i++)
        t[i] = cos(i * h); // Chebyshev interpolation points [-1, 1]

    double x1, x2, y1, y2;

    double w1i[PP];
    double dj[PP];
    int i, j, kk;
    dj[0] = 0.5;
    dj[P] = 0.5;
    for (j = 1; j < P; j++)
        dj[j] = 1.0;

    for (j = 0; j < PP; j++)
        w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];

    double a1i[PP];
    double a2j[PP];

    double x, y;
    double dx, dy;
    double SumA1;
    double SumA2;
    double D;
    double s;
    int a1exactind;
    int a2exactind;

    for (size_t tree_index = 1; tree_index < node_count; tree_index++) { // skip roots
        x1 = tree[tree_index].xinterval[0];
        x2 = tree[tree_index].xinterval[1];
        y1 = tree[tree_index].yinterval[0];
        y2 = tree[tree_index].yinterval[1];

        for (j = 0; j < PP; j++) {
            cList[tree_index].t1[j] = x1 + (t[j] + 1.0) * 0.5 * (x2 - x1);
            cList[tree_index].t2[j] = y1 + (t[j] + 1.0) * 0.5 * (y2 - y1);
        }

        size_t tp0 = tree[tree_index].members[0];
        size_t tp1 = tree[tree_index].members[1];
        size_t tp_j;

        for (tp_j = tp0; tp_j <= tp1; tp_j++) {
            x = pt_x[tp_j];
            y = pt_y[tp_j];

            a1exactind = -1;
            a2exactind = -1;

            SumA1 = 0.0;
            SumA2 = 0.0;

            for (j = 0; j < PP; j++) {
                dx = x - cList[tree_index].t1[j];
                dy = y - cList[tree_index].t2[j];
                if (fabs(dx) <= DBL_MIN)
                    a1exactind = j;
                else {
                    a1i[j] = w1i[j] / dx;
                    SumA1 += a1i[j];
                }

                if (fabs(dy) <= DBL_MIN)
                    a2exactind = j;
                else {
                    a2j[j] = w1i[j] / dy;
                    SumA2 += a2j[j];
                }
            } // j

            if (a1exactind > -1) {
                SumA1 = 1.0;
                for (j = 0; j < PP; j++)
                    a1i[j] = 0.0;
                a1i[a1exactind] = 1.0;
            }

            if (a2exactind > -1) {
                SumA2 = 1.0;
                for (j = 0; j < PP; j++)
                    a2j[j] = 0.0;
                a2j[a2exactind] = 1.0;
            }
            D = 1.0 / (SumA1 * SumA2);

            kk = -1;
            for (i = 0; i < PP; i++) {
                for (j = 0; j < PP; j++) {
                        kk++;
                        s = a1i[i] * a2j[j] * D;
                        cList[tree_index].moments[kk] += s * lambda[tp_j];
                } // j
            } // i
        } // tp_j

#if OPENACC_ENABLED
#pragma acc update device(cList[tree_index].t1[0:PP], cList[tree_index].t2[0:PP], cList[tree_index].moments[0:Pflat])
#endif
    } // tree_index
}


void U_Treecode::free_cluster_list()
{
    for (size_t tree_index = 0; tree_index < node_count; tree_index++) {
#if OPENACC_ENABLED
#pragma acc exit data delete(cList[tree_index].t1[0:PP], cList[tree_index].t2[0:PP], \
cList[tree_index].moments[0:Pflat])
#endif

        delete[] cList[tree_index].t1;
        delete[] cList[tree_index].t2;

        delete[] cList[tree_index].moments;
    }

#if OPENACC_ENABLED
#pragma acc exit data delete(cList[0:node_count])
#endif

    delete[] cList;
    cList = nullptr;
}


// =========================
// Original BL, DS
// =========================
void U_Treecode::Call_BL()
{
    // cout << "enter Call_BL" << endl;

#if OPENACC_ENABLED
#pragma acc kernels copyin(leaf_count) \
present(leaf_members[0:2][0:leaf_count], \
particles_x[0:numpars_s], particles_y[0:numpars_s], \
iList[0:leaf_count], \
velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s], \
cList)
{ // Begin acc kernels region
#endif

#if OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t batch_limit_1 = leaf_members[0][leaf_index];
        size_t batch_limit_2 = leaf_members[1][leaf_index];

        size_t far_list_size = iList[leaf_index].far_list_size;

        // cout << "far_list_size for leaf " << leaf_index << ": " << far_list_size <<endl;

#if OPENACC_ENABLED
        #pragma acc loop vector(128) independent
#endif
        for (size_t ii = batch_limit_1; ii <= batch_limit_2; ii++) {
            double tempx = 0.0;
            double tempy = 0.0;

            double p_x = particles_x[ii];
            double p_y = particles_y[ii];

#if OPENACC_ENABLED
            #pragma acc loop seq
#endif
            for (size_t jj = 0; jj < far_list_size; jj++) {
                size_t far_index = iList[leaf_index].far_list[jj];

#if OPENACC_ENABLED
                #pragma acc loop seq
#endif
                for (int kk = 0; kk < Pflat; kk++) {
                    int i = kk / PP;
                    int j = kk % PP;

                    // double xx = p_x - cList[far_index].t1[i];
                    // double yy = p_y - cList[far_index].t2[j];
                    // double R2 = xx * xx + yy * yy;
                    // tempx += cList[far_index].moments[kk] * (-yy / (R2 + delta2));
                    // tempy += cList[far_index].moments[kk] * (xx / (R2 + delta2));

                    double x_diff = p_x - cList[far_index].t1[i];
                    double y_diff = p_y - cList[far_index].t2[j];
                    double denom = cosh(2* pi / L * y_diff) - cos(2* pi / L * x_diff) + epsilon * epsilon;
                    tempx -= 0.5/L * sinh(2 * pi / L *  y_diff) / denom * cList[far_index].moments[kk];
                    tempy += 0.5/L * sin(2 * pi / L * x_diff) / denom * cList[far_index].moments[kk];
                } // kk
            } // jj

            velo_tc_reord_x[ii] += tempx;
            velo_tc_reord_y[ii] += tempy;
        } // ii
    } // leaf_index
#if OPENACC_ENABLED
} // End acc kernels region
#endif
}

void U_Treecode::Call_Ds()
{
    // cout << "enter Call_DS" << endl;
#if OPENACC_ENABLED
#pragma acc kernels copyin(leaf_count) \
present(leaf_members[0:2][0:leaf_count], \
tree_members[0:2][0:node_count], \
particles_x[0:numpars_s], particles_y[0:numpars_s], \
iList[0:leaf_count], \
lambda[0:numpars_s], \
velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s], \
cList)
{ // Begin acc kernels region
#endif

#if OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t limit_1_b = leaf_members[0][leaf_index];
        size_t limit_2_b = leaf_members[1][leaf_index];

        size_t near_list_size = iList[leaf_index].near_list_size;

        // cout << "near_list_size for leaf " << leaf_index << ": " << near_list_size <<endl;

#if OPENACC_ENABLED
        #pragma acc loop vector(128) independent
#endif
        for (size_t ii = limit_1_b; ii <= limit_2_b; ii++) {
            double tempx = 0.0;
            double tempy = 0.0;

#if OPENACC_ENABLED
            #pragma acc loop seq
#endif
            for (size_t kk = 0; kk < near_list_size; kk++) {
                size_t limit_1_c = tree_members[0][iList[leaf_index].near_list[kk]];
                size_t limit_2_c = tree_members[1][iList[leaf_index].near_list[kk]];

#if OPENACC_ENABLED
                #pragma acc loop seq
#endif
                for (size_t jj = limit_1_c; jj <= limit_2_c; jj++) {
                    // double ff0 = lambda[jj];
                    // double ff1 = lambda[jj];

                    // double xx = particles_x[ii] - particles_x[jj];
                    // double yy = particles_y[ii] - particles_y[jj];

                    // double R2 = xx * xx + yy * yy;

                    // tempx += -yy / (R2 + delta2) * ff0;
                    // tempy += xx / (R2 + delta2) * ff1;
                    double x_diff = particles_x[ii] - particles_x[jj];
                    double y_diff = particles_y[ii] - particles_y[jj];
                    double denom = cosh(2* pi / L * y_diff) - cos(2* pi / L * x_diff) + epsilon * epsilon;
                    tempx -= 0.5/L * sinh(2 * pi / L *  y_diff) / denom * lambda[jj];
                    tempy += 0.5/L * sin(2 * pi / L * x_diff) / denom * lambda[jj];
                } // jj
            } // kk

            velo_tc_reord_x[ii] += tempx;
            velo_tc_reord_y[ii] += tempy;
        } // ii
    } // size_t leaf_index
#if OPENACC_ENABLED
} // End acc kernels region
#endif
}


// =========================
// u1_grad BL, DS
// =========================

void U_Treecode::Call_BL_u1_grad()
{
#if OPENACC_ENABLED
#pragma acc kernels copyin(leaf_count) \
present(leaf_members[0:2][0:leaf_count], \
particles_x[0:numpars_s], particles_y[0:numpars_s], \
iList[0:leaf_count], \
velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s], \
cList)
{ // Begin acc kernels region
#endif

#if OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t batch_limit_1 = leaf_members[0][leaf_index];
        size_t batch_limit_2 = leaf_members[1][leaf_index];

        size_t far_list_size = iList[leaf_index].far_list_size;

        // cout << "far_list_size for leaf " << leaf_index << ": " << far_list_size <<endl;

#if OPENACC_ENABLED
        #pragma acc loop vector(128) independent
#endif
        for (size_t ii = batch_limit_1; ii <= batch_limit_2; ii++) {
            double tempx = 0.0;
            double tempy = 0.0;

            double p_x = particles_x[ii];
            double p_y = particles_y[ii];

#if OPENACC_ENABLED
            #pragma acc loop seq
#endif
            for (size_t jj = 0; jj < far_list_size; jj++) {
                size_t far_index = iList[leaf_index].far_list[jj];

#if OPENACC_ENABLED
                #pragma acc loop seq
#endif
                for (int kk = 0; kk < Pflat; kk++) {
                    int i = kk / PP;
                    int j = kk % PP;

                    double x_diff = 2*pi/ L *(p_x - cList[far_index].t1[i]);
                    double y_diff = 2*pi/ L *(p_y - cList[far_index].t2[j]);
                    double denom_sqr = (cosh(y_diff) - cos(x_diff) + epsilon * epsilon) * (cosh(y_diff) - cos(x_diff) + epsilon * epsilon);
                    double constant_c = pi/L/L;
                    tempx += constant_c * sinh(y_diff) * sin(x_diff) / denom_sqr * cList[far_index].moments[kk];
                    tempy += constant_c * (cos(x_diff) * cosh(y_diff) - 1 - epsilon * epsilon * cosh(y_diff)) / denom_sqr * cList[far_index].moments[kk];
                } // kk
            } // jj

            velo_tc_reord_x[ii] += tempx;
            velo_tc_reord_y[ii] += tempy;
        } // ii
    } // leaf_index
#if OPENACC_ENABLED
} // End acc kernels region
#endif
}

void U_Treecode::Call_DS_u1_grad()
{
#if OPENACC_ENABLED
#pragma acc kernels copyin(leaf_count) \
present(leaf_members[0:2][0:leaf_count], \
tree_members[0:2][0:node_count], \
particles_x[0:numpars_s], particles_y[0:numpars_s], \
iList[0:leaf_count], \
lambda[0:numpars_s], \
velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s], \
cList)
{ // Begin acc kernels region
#endif

#if OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t limit_1_b = leaf_members[0][leaf_index];
        size_t limit_2_b = leaf_members[1][leaf_index];

        size_t near_list_size = iList[leaf_index].near_list_size;

        // cout << "near_list_size for leaf " << leaf_index << ": " << near_list_size <<endl;

#if OPENACC_ENABLED
        #pragma acc loop vector(128) independent
#endif
        for (size_t ii = limit_1_b; ii <= limit_2_b; ii++) {
            double tempx = 0.0;
            double tempy = 0.0;

#if OPENACC_ENABLED
            #pragma acc loop seq
#endif
            for (size_t kk = 0; kk < near_list_size; kk++) {
                size_t limit_1_c = tree_members[0][iList[leaf_index].near_list[kk]];
                size_t limit_2_c = tree_members[1][iList[leaf_index].near_list[kk]];

#if OPENACC_ENABLED
                #pragma acc loop seq
#endif
                for (size_t jj = limit_1_c; jj <= limit_2_c; jj++) {
                    double x_diff = 2*pi/ L * (particles_x[ii] - particles_x[jj]);
                    double y_diff = 2*pi/ L * (particles_y[ii] - particles_y[jj]);
                    double denom_sqr = (cosh(y_diff) - cos(x_diff) + epsilon * epsilon) * (cosh(y_diff) - cos(x_diff) + epsilon * epsilon);
                    double constant_c = pi/L/L;
                    tempx += constant_c * sinh(y_diff) * sin(x_diff) / denom_sqr * lambda[jj];
                    tempy += constant_c * (cos(x_diff) * cosh(y_diff) - 1 - epsilon * epsilon * cosh(y_diff)) / denom_sqr * lambda[jj];
                } // jj
            } // kk

            velo_tc_reord_x[ii] += tempx;
            velo_tc_reord_y[ii] += tempy;
        } // ii
    } // size_t leaf_index
#if OPENACC_ENABLED
} // End acc kernels region
#endif
}


// =========================
// u2_grad BL, DS
// =========================

void U_Treecode::Call_BL_u2_grad()
{
#if OPENACC_ENABLED
#pragma acc kernels copyin(leaf_count) \
present(leaf_members[0:2][0:leaf_count], \
particles_x[0:numpars_s], particles_y[0:numpars_s], \
iList[0:leaf_count], \
velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s], \
cList)
{ // Begin acc kernels region
#endif

#if OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t batch_limit_1 = leaf_members[0][leaf_index];
        size_t batch_limit_2 = leaf_members[1][leaf_index];

        size_t far_list_size = iList[leaf_index].far_list_size;

        // cout << "far_list_size for leaf " << leaf_index << ": " << far_list_size <<endl;

#if OPENACC_ENABLED
        #pragma acc loop vector(128) independent
#endif
        for (size_t ii = batch_limit_1; ii <= batch_limit_2; ii++) {
            double tempx = 0.0;
            double tempy = 0.0;

            double p_x = particles_x[ii];
            double p_y = particles_y[ii];

#if OPENACC_ENABLED
            #pragma acc loop seq
#endif
            for (size_t jj = 0; jj < far_list_size; jj++) {
                size_t far_index = iList[leaf_index].far_list[jj];

#if OPENACC_ENABLED
                #pragma acc loop seq
#endif
                for (int kk = 0; kk < Pflat; kk++) {
                    int i = kk / PP;
                    int j = kk % PP;

                    double x_diff = 2*pi/ L *(p_x - cList[far_index].t1[i]);
                    double y_diff = 2*pi/ L *(p_y - cList[far_index].t2[j]);
                    double denom_sqr = (cosh(y_diff) - cos(x_diff) + epsilon * epsilon) * (cosh(y_diff) - cos(x_diff) + epsilon * epsilon);
                    double constant_c = pi/L/L;
                    tempx += constant_c * (cos(x_diff) * cosh(y_diff) - 1 + epsilon * epsilon * cos(x_diff)) / denom_sqr * cList[far_index].moments[kk];
                    tempy += -1 * constant_c * sin(x_diff) * sinh(y_diff) / denom_sqr * cList[far_index].moments[kk];
                } // kk
            } // jj

            velo_tc_reord_x[ii] += tempx;
            velo_tc_reord_y[ii] += tempy;
        } // ii
    } // leaf_index
#if OPENACC_ENABLED
} // End acc kernels region
#endif
}

void U_Treecode::Call_DS_u2_grad()
{
#if OPENACC_ENABLED
#pragma acc kernels copyin(leaf_count) \
present(leaf_members[0:2][0:leaf_count], \
tree_members[0:2][0:node_count], \
particles_x[0:numpars_s], particles_y[0:numpars_s], \
iList[0:leaf_count], \
lambda[0:numpars_s], \
velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s], \
cList)
{ // Begin acc kernels region
#endif

#if OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t limit_1_b = leaf_members[0][leaf_index];
        size_t limit_2_b = leaf_members[1][leaf_index];

        size_t near_list_size = iList[leaf_index].near_list_size;

        // cout << "near_list_size for leaf " << leaf_index << ": " << near_list_size <<endl;

#if OPENACC_ENABLED
        #pragma acc loop vector(128) independent
#endif
        for (size_t ii = limit_1_b; ii <= limit_2_b; ii++) {
            double tempx = 0.0;
            double tempy = 0.0;

#if OPENACC_ENABLED
            #pragma acc loop seq
#endif
            for (size_t kk = 0; kk < near_list_size; kk++) {
                size_t limit_1_c = tree_members[0][iList[leaf_index].near_list[kk]];
                size_t limit_2_c = tree_members[1][iList[leaf_index].near_list[kk]];

#if OPENACC_ENABLED
                #pragma acc loop seq
#endif
                for (size_t jj = limit_1_c; jj <= limit_2_c; jj++) {
                    double x_diff = 2*pi/ L * (particles_x[ii] - particles_x[jj]);
                    double y_diff = 2*pi/ L * (particles_y[ii] - particles_y[jj]);
                    double denom_sqr = (cosh(y_diff) - cos(x_diff) + epsilon * epsilon) * (cosh(y_diff) - cos(x_diff) + epsilon * epsilon);
                    double constant_c = pi/L/L;
                    tempx += constant_c * (cos(x_diff) * cosh(y_diff) - 1 + epsilon * epsilon * cos(x_diff)) / denom_sqr * lambda[jj];
                    tempy += -1 * constant_c * sin(x_diff) * sinh(y_diff) / denom_sqr * lambda[jj];
                } // jj
            } // kk

            velo_tc_reord_x[ii] += tempx;
            velo_tc_reord_y[ii] += tempy;
        } // ii
    } // size_t leaf_index
#if OPENACC_ENABLED
} // End acc kernels region
#endif
}


// =========================
// vorticity_grad BL, DS
// =========================

void U_Treecode::Call_BL_vorticity_grad()
{
#if OPENACC_ENABLED
#pragma acc kernels copyin(leaf_count) \
present(leaf_members[0:2][0:leaf_count], \
particles_x[0:numpars_s], particles_y[0:numpars_s], \
iList[0:leaf_count], \
velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s], \
cList)
{ // Begin acc kernels region
#endif

#if OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t batch_limit_1 = leaf_members[0][leaf_index];
        size_t batch_limit_2 = leaf_members[1][leaf_index];

        size_t far_list_size = iList[leaf_index].far_list_size;

        // cout << "far_list_size for leaf " << leaf_index << ": " << far_list_size <<endl;

#if OPENACC_ENABLED
        #pragma acc loop vector(128) independent
#endif
        for (size_t ii = batch_limit_1; ii <= batch_limit_2; ii++) {
            double tempx = 0.0;
            double tempy = 0.0;

            double p_x = particles_x[ii];
            double p_y = particles_y[ii];

#if OPENACC_ENABLED
            #pragma acc loop seq
#endif
            for (size_t jj = 0; jj < far_list_size; jj++) {
                size_t far_index = iList[leaf_index].far_list[jj];

#if OPENACC_ENABLED
                #pragma acc loop seq
#endif
                for (int kk = 0; kk < Pflat; kk++) {
                    int i = kk / PP;
                    int j = kk % PP;

                    double x_diff = 2*pi/ L *(p_x - cList[far_index].t1[i]);
                    double y_diff = 2*pi/ L *(p_y - cList[far_index].t2[j]);
                    double denom_cube = (cosh(y_diff) - cos(x_diff) + epsilon * epsilon) * (cosh(y_diff) - cos(x_diff) + epsilon * epsilon) * (cosh(y_diff) - cos(x_diff) + epsilon * epsilon);
                    double constant_c = 2 * pi * pi * epsilon * epsilon /L/L/L; 
                    tempx += constant_c * sin(x_diff) *(-3*cosh(y_diff) - cos(x_diff) - epsilon * epsilon) / denom_cube * cList[far_index].moments[kk];
                    tempy += constant_c * sinh(y_diff) *(-cosh(y_diff) - 3 * cos(x_diff) + epsilon * epsilon) / denom_cube * cList[far_index].moments[kk];
                } // kk
            } // jj

            velo_tc_reord_x[ii] += tempx;
            velo_tc_reord_y[ii] += tempy;
        } // ii
    } // leaf_index
#if OPENACC_ENABLED
} // End acc kernels region
#endif
}

void U_Treecode::Call_DS_vorticity_grad()
{
#if OPENACC_ENABLED
#pragma acc kernels copyin(leaf_count) \
present(leaf_members[0:2][0:leaf_count], \
tree_members[0:2][0:node_count], \
particles_x[0:numpars_s], particles_y[0:numpars_s], \
iList[0:leaf_count], \
lambda[0:numpars_s], \
velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s], \
cList)
{ // Begin acc kernels region
#endif

#if OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t limit_1_b = leaf_members[0][leaf_index];
        size_t limit_2_b = leaf_members[1][leaf_index];

        size_t near_list_size = iList[leaf_index].near_list_size;

        // cout << "near_list_size for leaf " << leaf_index << ": " << near_list_size <<endl;

#if OPENACC_ENABLED
        #pragma acc loop vector(128) independent
#endif
        for (size_t ii = limit_1_b; ii <= limit_2_b; ii++) {
            double tempx = 0.0;
            double tempy = 0.0;

#if OPENACC_ENABLED
            #pragma acc loop seq
#endif
            for (size_t kk = 0; kk < near_list_size; kk++) {
                size_t limit_1_c = tree_members[0][iList[leaf_index].near_list[kk]];
                size_t limit_2_c = tree_members[1][iList[leaf_index].near_list[kk]];

#if OPENACC_ENABLED
                #pragma acc loop seq
#endif
                for (size_t jj = limit_1_c; jj <= limit_2_c; jj++) {
                    double x_diff = 2*pi/ L * (particles_x[ii] - particles_x[jj]);
                    double y_diff = 2*pi/ L * (particles_y[ii] - particles_y[jj]);
                    double denom_cube = (cosh(y_diff) - cos(x_diff) + epsilon * epsilon) * (cosh(y_diff) - cos(x_diff) + epsilon * epsilon) * (cosh(y_diff) - cos(x_diff) + epsilon * epsilon);
                    double constant_c = 2 * pi * pi * epsilon * epsilon /L/L/L; 
                    tempx += constant_c * sin(x_diff) *(-3*cosh(y_diff) - cos(x_diff) - epsilon * epsilon) / denom_cube * lambda[jj];
                    tempy += constant_c * sinh(y_diff) *(-cosh(y_diff) - 3 * cos(x_diff) + epsilon * epsilon) / denom_cube * lambda[jj];
                } // jj
            } // kk

            velo_tc_reord_x[ii] += tempx;
            velo_tc_reord_y[ii] += tempy;
        } // ii
    } // size_t leaf_index
#if OPENACC_ENABLED
} // End acc kernels region
#endif
}


// =========================
// laplacian BL, DS
// =========================

void U_Treecode::Call_BL_laplacian()
{
#if OPENACC_ENABLED
#pragma acc kernels copyin(leaf_count) \
present(leaf_members[0:2][0:leaf_count], \
particles_x[0:numpars_s], particles_y[0:numpars_s], \
iList[0:leaf_count], \
velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s], \
cList)
{ // Begin acc kernels region
#endif

#if OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t batch_limit_1 = leaf_members[0][leaf_index];
        size_t batch_limit_2 = leaf_members[1][leaf_index];

        size_t far_list_size = iList[leaf_index].far_list_size;

        // cout << "far_list_size for leaf " << leaf_index << ": " << far_list_size <<endl;

#if OPENACC_ENABLED
        #pragma acc loop vector(128) independent
#endif
        for (size_t ii = batch_limit_1; ii <= batch_limit_2; ii++) {
            double tempx = 0.0;
            double tempy = 0.0;

            double p_x = particles_x[ii];
            double p_y = particles_y[ii];

#if OPENACC_ENABLED
            #pragma acc loop seq
#endif
            for (size_t jj = 0; jj < far_list_size; jj++) {
                size_t far_index = iList[leaf_index].far_list[jj];

#if OPENACC_ENABLED
                #pragma acc loop seq
#endif
                for (int kk = 0; kk < Pflat; kk++) {
                    int i = kk / PP;
                    int j = kk % PP;

                    double x_diff = 2*pi/ L *(p_x - cList[far_index].t1[i]);
                    double y_diff = 2*pi/ L *(p_y - cList[far_index].t2[j]);
                    double C = cos(x_diff);
                    double S = sin(x_diff);
                    double H = cosh(y_diff);
                    double Sh = sinh(y_diff);
                    double eps_sqr = epsilon * epsilon;
                    double eps_4 = epsilon * epsilon * epsilon * epsilon;
                    double denom_4 = (H - C + eps_sqr)*(H - C + eps_sqr)*(H - C + eps_sqr)*(H - C + eps_sqr);
                    double constant_c = 4 * pi * pi * pi * epsilon * epsilon /L/L/L/L; 
                    tempx += constant_c * (C*C*C +5*C*C*H -5*C*H*H -8*C*H*eps_sqr +2*C*S*S + 10*C*Sh*Sh -C*eps_4
                            -H*H*H + 10*H*S*S +2*H*Sh*Sh + H*eps_4 + 4*S*S*eps_sqr - 4*Sh*Sh*eps_sqr) 
                            / denom_4 * cList[far_index].moments[kk];               
                    tempy += 0;
                } // kk
            } // jj

            velo_tc_reord_x[ii] += tempx;
            velo_tc_reord_y[ii] += tempy;
        } // ii
    } // leaf_index
#if OPENACC_ENABLED
} // End acc kernels region
#endif
}

void U_Treecode::Call_DS_laplacian()
{
#if OPENACC_ENABLED
#pragma acc kernels copyin(leaf_count) \
present(leaf_members[0:2][0:leaf_count], \
tree_members[0:2][0:node_count], \
particles_x[0:numpars_s], particles_y[0:numpars_s], \
iList[0:leaf_count], \
lambda[0:numpars_s], \
velo_tc_reord_x[0:numpars_s], velo_tc_reord_y[0:numpars_s], \
cList)
{ // Begin acc kernels region
#endif

#if OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t limit_1_b = leaf_members[0][leaf_index];
        size_t limit_2_b = leaf_members[1][leaf_index];

        size_t near_list_size = iList[leaf_index].near_list_size;

        // cout << "near_list_size for leaf " << leaf_index << ": " << near_list_size <<endl;

#if OPENACC_ENABLED
        #pragma acc loop vector(128) independent
#endif
        for (size_t ii = limit_1_b; ii <= limit_2_b; ii++) {
            double tempx = 0.0;
            double tempy = 0.0;

#if OPENACC_ENABLED
            #pragma acc loop seq
#endif
            for (size_t kk = 0; kk < near_list_size; kk++) {
                size_t limit_1_c = tree_members[0][iList[leaf_index].near_list[kk]];
                size_t limit_2_c = tree_members[1][iList[leaf_index].near_list[kk]];

#if OPENACC_ENABLED
                #pragma acc loop seq
#endif
                for (size_t jj = limit_1_c; jj <= limit_2_c; jj++) {
                    double x_diff = 2*pi/ L * (particles_x[ii] - particles_x[jj]);
                    double y_diff = 2*pi/ L * (particles_y[ii] - particles_y[jj]);
                    double C = cos(x_diff);
                    double S = sin(x_diff);
                    double H = cosh(y_diff);
                    double Sh = sinh(y_diff);
                    double eps_sqr = epsilon * epsilon;
                    double eps_4 = epsilon * epsilon * epsilon * epsilon;
                    double denom_4 = (H - C + eps_sqr)*(H - C + eps_sqr)*(H - C + eps_sqr)*(H - C + eps_sqr);
                    double constant_c = 4 * pi * pi * pi * epsilon * epsilon /L/L/L/L; 

                    tempx += constant_c * (C*C*C +5*C*C*H -5*C*H*H -8*C*H*eps_sqr +2*C*S*S + 10*C*Sh*Sh -C*eps_4
                            -H*H*H + 10*H*S*S +2*H*Sh*Sh + H*eps_4 + 4*S*S*eps_sqr - 4*Sh*Sh*eps_sqr) 
                            / denom_4 * lambda[jj];
                    tempy += 0;
                } // jj
            } // kk

            velo_tc_reord_x[ii] += tempx;
            velo_tc_reord_y[ii] += tempy;
        } // ii
    } // size_t leaf_index
#if OPENACC_ENABLED
} // End acc kernels region
#endif
}





void U_Treecode::Compute_SUM()
{
    switch (mode)
    {
        case original:
            // cout << "orginal kernel: " <<endl;
            Call_BL();
            Call_Ds();
            break;

        case u1_grad: // for u1s_grad_x, u1s_grad_y, b1s_grad_x, b1s_grad_y
            // cout << "u1_grad kernel: " <<endl;
            Call_BL_u1_grad();
            Call_DS_u1_grad();
            break;

        case u2_grad: // for u2s_grad_x, u2s_grad_y, b2s_grad_x, b2s_grad_y
            // cout << "u2_grad kernel: " <<endl;
            Call_BL_u2_grad();
            Call_DS_u2_grad();
            break;

        case vorticity_grad: // for vorticity_grad_x, vorticity_grad_x, j_grad_x, j_grad_y
            // cout << "vorticity_grad kernel: " <<endl;
            Call_BL_vorticity_grad();
            Call_DS_vorticity_grad();
            break;

        case laplacian: // for vorticity_laplacian, j_laplacian;
            // cout << "laplacian kernel: " <<endl;
            Call_BL_laplacian();
            Call_DS_laplacian();
            break;

        default:
            // cout << "orginal kernel: " <<endl;
            Call_BL();
            Call_Ds();
            break;
    }
}



// ===============================
// helpers
// ===============================

long U_Treecode::getTickCount() {
    tms tm;
    return times(&tm);
}

double U_Treecode::minval(const double* x, size_t len) {
    double MinVal = x[0];
    for (size_t i = 1; i < len; i++) {
        if (MinVal > x[i]) MinVal = x[i];
    }
    return MinVal - 0.1234; 
}

double U_Treecode::maxval(const double* x, size_t len) {
    double MaxVal = x[0];
    for (size_t i = 1; i < len; i++) {
        if (MaxVal < x[i]) MaxVal = x[i];
    }
    return MaxVal + 0.001; 
}

