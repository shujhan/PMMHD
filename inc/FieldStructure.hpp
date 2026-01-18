#ifndef FIELD_STRUCTURE_HPP
#define FIELD_STRUCTURE_HPP
#include <math.h>
#include <vector>
#include <cstring>
using namespace std;

class Field {
    public: 
        virtual void operator()     (double* e1s, double* e2s, double* x_vals, int nx, 
                                    double* y_vals, double* q_ws, int ny) = 0;
        virtual void print_field_obj() = 0;
        virtual ~Field();
};

class U_DirectSum : public Field {
    public:
        U_DirectSum();
        U_DirectSum(double L, double epsilon);
        double epsilon;
        double L;
        void operator() (double* e1s, double* e2s, double* x_vals, int nx, 
                        double* y_vals, double* q_ws, int ny);
        void print_field_obj();
        ~U_DirectSum();
};

class U_Treecode : public Field {
    public:
        U_Treecode();
        U_Treecode(double L, double epsilon,
                double mac, int degree, int max_source, int max_target);
        void operator() (double* e1s, double* e2s, double* x_vals, int nx, 
                        double* y_vals, double* q_ws, int ny);
        void print_field_obj();
        ~U_Treecode();


    private:
        struct panel
        {
            size_t members[2];
            double xinterval[2];
            double yinterval[2];
            double xc; // Panel center x coordinate
            double yc; // Panel center y coordinate
            double MAC; // r^2 / theta^2
            vector<size_t> children;

            // Initialization
            panel() : xc(0.0), yc(0.0), MAC(0.0)
            {
                memset(members, 0, sizeof(members));
                memset(xinterval, 0, sizeof(xinterval));
                memset(yinterval, 0, sizeof(yinterval));
            }
        };

        struct cluster_list
        {
            double* t1; // Interpolation points in x direction
            double* t2;
            double* moments;

            // Initialization
            cluster_list(): t1(NULL), t2(NULL)
            {
                moments = NULL;
            }
        };

        struct interaction_list
        {
            size_t* far_list;
            size_t far_list_size;
            size_t* near_list;
            size_t near_list_size;

            // Initialization
            interaction_list(): far_list(NULL), far_list_size(0), near_list(NULL), near_list_size(0){}
        };

        // =========================
        // functions
        // =========================
        void compute_RHS_BLTC();   // build tree + compute far/near sum (treecode)
        // void compute_RHS_DS();     // direct sum (debug/verification)
        void cleanup();            // free everything

        // =========================
        // helpers
        // =========================
        long getTickCount();

        double minval(const double* x, size_t len);
        double maxval(const double* x, size_t len);

        // =========================
        // tree building
        // =========================
        void build_tree_init();
        void build_tree_2D_Recursive(size_t panel_index, int level,
                                    double* pt_x, double* pt_y,
                                    int* pt_index, int* pt_old_index);

        void split_tree_node(size_t panel_index,
                            double* pt_x, double* pt_y,
                            int* pt_index, int* pt_old_index);

        void split_2(size_t panel_index, int split_code,
                    double* pt_x, double* pt_y,
                    int* pt_index, int* pt_old_index);

        void split_4(size_t panel_index,
                    double* pt_x, double* pt_y,
                    int* pt_index, int* pt_old_index);

        void Swap(size_t i, size_t j,
                double* pt_x, double* pt_y,
                int* pt_index, int* pt_old_index);

        // =========================
        // interaction list + cluster list
        // =========================
        void build_interaction_list(size_t leaf_index, size_t panel_index,
                                    std::vector<std::vector<size_t>>& Interaction_List_far,
                                    std::vector<std::vector<size_t>>& Interaction_List_near);

        void alloc_set_interaction_list(const std::vector<std::vector<size_t>>& Interaction_List_far,
                                        const std::vector<std::vector<size_t>>& Interaction_List_near);
        void free_interaction_list();

        void alloc_set_cluster_list(double* pt_x, double* pt_y);
        void free_cluster_list();

        // =========================
        // summation kernels (same math as your code)
        // =========================
        void Compute_SUM();
        void Call_BL(); // far field
        void Call_Ds(); // near field direct


        // =========================
        // Parameters
        // =========================
        double L;
        double epsilon;

        double mac;
        int degree;
        int max_source;
        int max_target;

        // Tree parameters
        size_t N0; //max_

        // Chebyshev / interpolation
        int P;
        int PP;
        int Pflat;

        // MAC theta^2
        double sq_theta;
        const double pi = 3.14159265358979323846;

        size_t numpars_s = 0;

        // =========================
        // Data arrays
        // =========================
        double* lambda;       // weights
        double* particles_x;  // sources x (reordered)
        double* particles_y;  // sources y (reordered)

        // // direct sum result
        // double* velo_ds_x;
        // double* velo_ds_y;

        // treecode velocity reordered
        double* velo_tc_reord_x;
        double* velo_tc_reord_y;

        // treecode velocity in original order
        double* velo_tc_noreord_x;
        double* velo_tc_noreord_y;

        // tree data
        size_t node_count;
        size_t leaf_count;

        std::vector<panel> tree;
        std::vector<size_t> leaf;

        // members arrays for device kernels
        size_t* tree_members[2];
        size_t* leaf_members[2];

        // interaction list and cluster list
        interaction_list* iList;
        cluster_list* cList;

        // bounding box
        double xyminmax[4];
};


#endif /* FIELD_STRUCTURE_HPP */
