#ifndef PANEL
#define PANEL


#include <iostream> // cout, endl
using namespace std;
    /**
     * @brief ordering of points and child indices
     * 
     * Stored in the format
     *  2 ----- 5 ----- 8
     *  |  [1]  |  [3]  |
     *  1 ----- 4 ----- 7
     *  |  [0]  |  [2]  |
     *  0 ----- 3 ----- 6
     * 
     *  y refined panel
     * Stored in the format
     *  2 ----- 5 ----- 8
     *  |      [1]      |
     *  1 ----- 4 ----- 7
     *  |      [0]      |
     *  0 ----- 3 ----- 6
     * 
     *
     * -----
     * Root panel `parent_ind` and `which_child` are -1
     * if panel is not refined, `child_ind_start` is -1
     * Neighbor indices are -1 if not determined yet
     * Neighbor indices are -2 if panel is a boundary and can't have neighbors
    */

struct Panel {
    int panel_ind;
    int point_inds[9]; // point index in vector xs and ys 
    int level;
    int parent_ind;
    int which_child;
    int left_nbr_ind, top_nbr_ind, right_nbr_ind, bottom_nbr_ind;
    bool is_left_bdry, is_right_bdry;
    bool is_top_bdry, is_bottom_bdry;
    bool needs_refinement;
    bool is_refined_xy;
    bool is_refined_y;
    int child_inds_start;



    // Constructor definitions
    Panel();
    Panel(int panel_ind);

    Panel(int panel_ind, int level, int parent_ind, int which_child,
        int ln_ind, int tn_ind, 
        int rn_ind, int bn_ind);

    Panel(int panel_ind, int level, int parent_ind, int which_child,
        int p0, int p1, int p2, int p3, int p4, int p5, 
        int p6, int p7, int p8, 
        int ln_ind, int tn_ind, 
        int rn_ind, int bn_ind,
        bool is_left_bdry, bool is_right_bdry,
        bool is_top_bdry, bool is_bottom_bdry);


    // End Constructors
    void set_point_inds(int p0, int p1, int p2, int p3, int p4, 
                        int p5, int p6, int p7, int p8);

    void set_child_inds_start(int c0);
    void set_child_inds_start(int c0, bool refined_y);


    void print_panel() const;

};

std::ostream& operator<<(std::ostream& os, const Panel& panel);

#endif