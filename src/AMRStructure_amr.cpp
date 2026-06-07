#include "AMRStructure.hpp"


void AMRStructure::generate_mesh(std::function<double (double,double)> f0, 
                                std::function<double (double,double)> f1, 
                                 bool do_adaptively_refine_vorticity, bool do_adaptively_refine_j, bool is_initial_step) 
{
    bool verbose=false;

    // auto start = high_resolution_clock::now();
    create_prerefined_mesh(is_initial_step);
    // auto stop = high_resolution_clock::now();
    // add_time(tree_build_time,  duration_cast<duration<double>>(stop - start) );


    if (is_initial_step) {
        for (int ii = 0; ii < xs.size(); ii++) {
            w0s[ii] = (*w0)(xs[ii],ys[ii]);
            j0s[ii] = (*j0)(xs[ii],ys[ii]);
            q_plus[ii] = w0s[ii] + j0s[ii];
            q_minus[ii] = w0s[ii] - j0s[ii];
        }
    } else {
        int nx_points = 2*npanels_x + 1;
        int ny_points = 2*npanels_y + 1;

        #ifdef DEBUG
        cout << "interpolating to grid " << endl;
        #endif
        
        // Remesh each deformed Lagrangian copy back onto the fresh uniform grid.
        // Source = deformed mesh (old_panels + old_xs/old_ys + old_q0s),
        // target = the freshly built uniform grid points (xs, ys).
        // interpolate q_plus and q_minus back to uniform panel mesh so that they are at the same location, then calculate w0 and j0
        old_xs = old_xs_plus;  old_ys = old_ys_plus;  old_q0s = old_q_plus;
        interpolate_to_initial_xys(q_plus, xs, ys, nx_points, ny_points);
        old_xs = old_xs_minus; old_ys = old_ys_minus; old_q0s = old_q_minus;
        interpolate_to_initial_xys(q_minus, xs, ys, nx_points, ny_points);
        for (int ii = 0; ii < xs.size(); ii++) {
            w0s[ii] = 0.5 * (q_plus[ii] + q_minus[ii]);
            j0s[ii] = 0.5 * (q_plus[ii] - q_minus[ii]);
        }
        #ifdef DEBUG
        cout << "done interpolating to grid" << endl;
        #endif
    }

    int num_panels_pre_refine = panels.size();



    if (do_adaptively_refine_vorticity || do_adaptively_refine_j) {
        // start = high_resolution_clock::now();

        for (int ii = minimum_unrefined_index; ii < panels.size(); ++ii) {
            test_panel(ii, verbose);
        }
        // stop = high_resolution_clock::now();
        // add_time(panel_test_time,  duration_cast<duration<double>>(stop - start) );


        while (need_further_refinement) {
            need_further_refinement = false;
            // auto amr_start = high_resolution_clock::now();
            refine_panels(f0, true, is_initial_step);
            if (is_initial_step) {
                for (int ii = 0; ii < xs.size(); ii++) {
                    w0s[ii] = (*w0)(xs[ii],ys[ii]);
                    j0s[ii] = (*j0)(xs[ii],ys[ii]);
                    q_plus[ii] = w0s[ii] + j0s[ii];
                    q_minus[ii] = w0s[ii] - j0s[ii];
                }
            } else {
                // int nx_points = 2*npanels_x + 1;
                // int ny_points = 2*npanels_y + 1;
                // old_xs = old_xs_plus;  old_ys = old_ys_plus;  old_q0s = old_q_plus;
                // interpolate_to_initial_xys(q_plus, xs, ys, nx_points, ny_points);
                // old_xs = old_xs_minus; old_ys = old_ys_minus; old_q0s = old_q_minus;
                // interpolate_to_initial_xys(q_minus, xs, ys, nx_points, ny_points);
                // for (int ii = 0; ii < xs.size(); ii++) {
                //     w0s[ii] = 0.5 * (q_plus[ii] + q_minus[ii]);
                //     j0s[ii] = 0.5 * (q_plus[ii] - q_minus[ii]);
                // }

                // // Robustly value the points created this pass (and re-value the
                // // rest, harmlessly) through the corner-safe neighbor-walk, so
                // // test_panel sees correct w0/j0 and the final mesh is already
                // // filled -- no separate post-loop pass needed.
                // old_xs = old_xs_plus;  old_ys = old_ys_plus;  old_q0s = old_q_plus;
                // interpolate_q_scattered(q_plus);
                // old_xs = old_xs_minus; old_ys = old_ys_minus; old_q0s = old_q_minus;
                // interpolate_q_scattered(q_minus);
                // for (int ii = 0; ii < xs.size(); ii++) {
                //     w0s[ii] = 0.5 * (q_plus[ii] + q_minus[ii]);
                //     j0s[ii] = 0.5 * (q_plus[ii] - q_minus[ii]);
                // }
                for (int ii = 0; ii < xs.size(); ii++) {
                    old_xs = old_xs_plus;  old_ys = old_ys_plus;  old_q0s = old_q_plus;
                    double qp = interpolate_from_mesh(new_xs.at(ii), new_ys.at(ii), false);
                    old_xs = old_xs_minus; old_ys = old_ys_minus; old_q0s = old_q_minus;
                    double qm = interpolate_from_mesh(new_xs.at(ii), new_ys.at(ii), false);
                    q_plus[ii] = qp;
                    q_minus[ii] = qm;
                    w0s[ii] = 0.5 * (q_plus[ii] + q_minus[ii]);
                    j0s[ii] = 0.5 * (q_plus[ii] - q_minus[ii]);
                }
            }
            // auto amr_stop = high_resolution_clock::now();
            // add_time(amr_refine_time, duration_cast<duration<double>>(amr_stop-amr_start) );

            // amr_start = high_resolution_clock::now();
            for (int ii = minimum_unrefined_index; ii < panels.size(); ++ii) {
                if (!panels[ii].is_refined_xy) {
                    test_panel(ii, verbose);
                }
            }
            // amr_stop = high_resolution_clock::now();
            // add_time(amr_test_time,  duration_cast<duration<double>>(amr_stop - amr_start) );
        }
        // stop = high_resolution_clock::now();
        // add_time(tree_build_time,  duration_cast<duration<double>>(stop - start) );
  
    }

    set_leaves_weights();
}


int AMRStructure::create_prerefined_mesh(bool is_initial_step) {
    // create a level 1 panel, 2^2 + 1 = 5 points in x and y
    // in total 5 * 5 = 25 points 

    if (initial_height < 1) {
        throw std::invalid_argument("height must be greater than 1");
    }

    double dx = Lx / 4;
    double dy = Ly / 4;
    std::vector<double> xs_init, ys_init; 
    // get the x, y value for level 0 panel
    for(int ii = 0; ii < 5; ++ii) {
        xs_init.push_back(x_min + ii * dx);
        ys_init.push_back(y_min + ii * dy);
    }

    panels.clear();
    xs.clear();
    ys.clear();
    w0s.clear();
    j0s.clear();
    q_plus.clear();
    q_minus.clear();
    xs.reserve(25);
    ys.reserve(25);
    // put 9 points to xs and ys: level 0 panel 
    for (int ii = 0; ii < 5; ii += 2) {
        for (int jj = 0; jj < 5; jj+=2) { 
            xs.push_back(xs_init[ii]); 
            ys.push_back(ys_init[jj]);
        }
    }
    // put points for level 1 panel 
    for (int ii = 0; ii < 2; ++ii) {
        xs.push_back(xs_init[2*ii]); xs.push_back(xs_init[2*ii]);
        ys.push_back(ys_init[1]); ys.push_back(ys_init[3]);
        for (int jj = 0; jj < 5; ++jj) {
            xs.push_back(xs_init[1 + 2*ii]);
            ys.push_back(ys_init[jj]);
        }
    }
    // finish level 1 panel 
    for (int jj = 1; jj < 5; jj+=2) {
        xs.push_back(xs_init[4]);
        ys.push_back(ys_init[jj]);
    }

    w0s = std::vector<double>(25, 1.0);
    j0s = std::vector<double>(25, 1.0);
    q_plus = std::vector<double>(25, 1.0);
    q_minus = std::vector<double>(25, 1.0);


    if (bcs == 1) { // open in y 
        panels.push_back(Panel{});
        panels[0].is_left_bdry = true;
        panels[0].is_right_bdry = true;
        panels.push_back(Panel{1, 1, 0, 0, 3, 2, 3, -2});
        panels[1].is_left_bdry = true;
        panels.push_back(Panel{2, 1, 0, 1, 4, -2, 4, 1});
        panels[2].is_left_bdry = true;
        panels.push_back(Panel{3, 1, 0, 2, 1, 4, 1, -2});
        panels[3].is_right_bdry = true;
        panels.push_back(Panel{4, 1, 0, 3, 2,-2,2,3});
        panels[4].is_right_bdry = true;
    } else if (bcs == 0) { // periodic in y 
        panels.push_back(Panel{});
        panels[0].is_left_bdry = true;
        panels[0].is_right_bdry = true;
        panels[0].is_top_bdry = true;
        panels[0].is_bottom_bdry = true;
        panels.push_back(Panel{1, 1, 0, 0, 3, 2, 3, 2});
        panels[1].is_left_bdry = true;
        panels[1].is_bottom_bdry = true;
        panels.push_back(Panel{2, 1, 0, 1, 4, 1, 4, 1});
        panels[2].is_left_bdry = true;
        panels[2].is_top_bdry = true;
        panels.push_back(Panel{3, 1, 0, 2, 1, 4, 1, 4});
        panels[3].is_right_bdry = true;
        panels[3].is_bottom_bdry = true;
        panels.push_back(Panel{4, 1, 0, 3, 2, 3, 2,3});
        panels[4].is_right_bdry = true;
        panels[4].is_top_bdry = true;
    }

    panels[1].set_point_inds(0,9,1,11,12,13,3,16,4);
    panels[1].needs_refinement = true;
    panels[2].set_point_inds(1,10,2,13,14,15,4,17,5);
    panels[2].needs_refinement = true;
    panels[3].set_point_inds(3,16,4,18,19,20,6,23,7);
    panels[3].needs_refinement = true;
    panels[4].set_point_inds(4,17,5,20,21,22,7,24,8);
    panels[4].needs_refinement = true;

    panels[0].set_child_inds_start(1);
    minimum_unrefined_index = 1;

    // finish level 1, 25 points and 5 panels in total 


    // call refine
    // first do uniform mesh 
    for (int level = 1; level < initial_height; ++level) {
        int num_panels_pre_refine = panels.size();

        for (auto panel_it = panels.begin() + minimum_unrefined_index; panel_it != panels.end(); ++panel_it) {
            panel_it->needs_refinement = true;
        }
        refine_panels( [] (double x, double y) {return 1.0;} , false, is_initial_step);
        minimum_unrefined_index = num_panels_pre_refine;
    }

    for (int level = 0; level < y_height; ++level) {
        int num_panels_pre_refine = panels.size();

        for (auto panel_it = panels.begin() + minimum_unrefined_index; panel_it != panels.end(); ++panel_it) {
            panel_it->needs_refinement = true;
        }
        refine_panels_refine_v( [] (double x, double v) {return 1.0;} , false, is_initial_step);
        minimum_unrefined_index = num_panels_pre_refine;
    }

    is_initial_mesh_set = true;
    return 0;

}


void AMRStructure::refine_panels_refine_v(std::function<double (double,double)> f, bool do_adaptive_refine, bool is_initial_step) {

    // Note: this assumes that we are refining in v uniformly before any xp refinement;
    // No compatibility with xp refined panels is guaranteed
    std::vector <double> new_xs;
    std::vector <double> new_ys;
    std::vector <double> new_w0s;
    std::vector <double> new_j0s;
    std::vector <double> new_q_plus;
    std::vector <double> new_q_minus;
    std::vector <int> prospective_leaf_inds;
    // int new_vert_ind = particles.size();
    int new_vert_ind = xs.size();
    int num_panels_before_this_iter = panels.size();
    


    for (int jj = minimum_unrefined_index; jj < num_panels_before_this_iter; ++jj) {
        Panel* panel= &(panels[jj]);
        
        if (panel->needs_refinement ) {
            std::vector<double> panel_xs;
            std::vector<double> panel_ys;
            double dx, dy;

            const int* panel_points = panel->point_inds;
            for (int ii = 0; ii < 9; ++ii) {
                int point_ind = panel_points[ii];
                panel_xs.push_back(xs[point_ind]);
                panel_ys.push_back(ys[point_ind]);
            }
            dx = panel_xs[3] - panel_xs[0];
            dy = panel_ys[1] - panel_ys[0];
            double sub_dy = 0.5 * dy;
            double sub_dx = 0.5 * dx;

            int num_new_panels = panels.size();
            double subpanel_xs[5], subpanel_ys[5];

            for (int ii = 0; ii < 5; ii ++) {
                subpanel_ys[ii] = panel_ys[0] + sub_dy * ii;
                subpanel_xs[ii] = panel_xs[0] + sub_dx * ii;
            }
            //   2        5         8 (2 by periodic bcs)
            
            //   10      12 [2]    14 (10 by pbcs)
            
            //   1        4         7 (1 by periodic bcs)
            
            //   9       11 [1]    13 (9 by pbcs)
            
            //   0        3         6 (0 by periodic bcs)

            int point_9_ind, point_10_ind, point_13_ind, point_14_ind;
            int child_0_bottom_nbr_ind = -1;
            int child_0_left_nbr_ind = num_new_panels; 
            int child_0_right_nbr_ind = child_0_left_nbr_ind;
            int child_1_left_nbr_ind = num_new_panels + 1;
            int child_1_top_nbr_ind = -1;
            int child_1_right_nbr_ind = child_1_left_nbr_ind;

            // generate new vertices

            Panel* panel_parent;
            // check left neighbor
            if (panel->left_nbr_ind == -2) {
                child_0_left_nbr_ind = -2;
                child_1_left_nbr_ind = -2;
                point_9_ind = new_vert_ind++;
                point_10_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
            } else if (panel->left_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                Panel* parent_left = &(panels[panel_parent->left_nbr_ind]);
                if (! (parent_left->is_refined_xy || parent_left->is_refined_y) ) {
                    parent_left->needs_refinement = true;
                    need_further_refinement = true;
                    // cout << "refine: setting refinement flag in panel " << jj << endl;
                }
                point_9_ind = new_vert_ind;
                // point_10_ind = new_vert_ind++;
                point_10_ind = point_9_ind + 1;
                new_vert_ind += 2;
                new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
            } else {
                Panel* panel_left = &(panels[panel->left_nbr_ind]);
                if (! (panel_left->is_refined_xy || panel_left->is_refined_y) ) {
                    point_9_ind = new_vert_ind++;
                    point_10_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                    new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
                }
                else {
                    if (panel_left->is_refined_xy) {
                        child_0_left_nbr_ind = panel_left->child_inds_start +2;
                        child_1_left_nbr_ind = panel_left->child_inds_start + 3;
                    } else { // panel_left is refined in v
                        child_0_left_nbr_ind = panel_left->child_inds_start +0;
                        child_1_left_nbr_ind = panel_left->child_inds_start + 1;
                    }
                    Panel* child_0_left_nbr = &(panels[child_0_left_nbr_ind]);
                    child_0_left_nbr->right_nbr_ind = num_new_panels;
                    Panel* child_1_left_nbr = &(panels[child_1_left_nbr_ind]);
                    child_1_left_nbr->right_nbr_ind = num_new_panels + 1;
                    // if (panel->is_left_bdry && bcs==periodic_bcs) {
                    if (panel->is_left_bdry) {
                        point_9_ind = new_vert_ind++;
                        point_10_ind = new_vert_ind++;
                        new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                        new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
                    } else {
                        point_9_ind = child_0_left_nbr->point_inds[7];
                        point_10_ind = child_1_left_nbr->point_inds[7];
                    }
                }
            }
            
            // check bottom neighbor
            int bottom_nbr_ind = panel->bottom_nbr_ind;
            if (bottom_nbr_ind == -2) {
                // Hard boundary: no periodicity in y; child 0 has no bottom neighbor.
                // Note: refine_v on its own does not introduce new bottom-edge vertices
                // because the bottom edge of the parent IS the bottom edge of child 0
                // (no new mid-x vertex needed below the parent).
                child_0_bottom_nbr_ind = -2;
            } else if (bottom_nbr_ind == -1) {
                // Bottom neighbor at this level is not yet refined: walk up to parent
                // and flag the parent's bottom neighbor for refinement so the cascade
                // can continue. This mirrors the behavior in refine_panels.
                panel_parent = &(panels[panel->parent_ind]);
                Panel* parent_bottom = &(panels[panel_parent->bottom_nbr_ind]);
                if (!(parent_bottom->is_refined_xy || parent_bottom->is_refined_y)) {
                    parent_bottom->needs_refinement = true;
                    need_further_refinement = true;
                    #ifdef DEBUG
                    cout << "refine_v: setting refinement flag on parent's bottom nbr, from panel " << jj << endl;
                    #endif
                }
                child_0_bottom_nbr_ind = -1;
            } else {
                Panel* panel_bottom = &(panels[bottom_nbr_ind]);
                if (panel_bottom->is_refined_xy) {
                    // Design intent: don't mix xy-refined neighbors with v-only refinement.
                    cout << "Shouldn't be allowed to call refine in v if bottom neighbor is refined in x and v!" << endl;
                    child_0_bottom_nbr_ind = bottom_nbr_ind;  // best-effort fallback
                }
                else {
                    if (!panel_bottom->is_refined_y) {
                        // Neighbor exists at our level but isn't subdivided.
                        // The new child 0's bottom neighbor is the whole panel_bottom.
                        // Update panel_bottom's top neighbor to point at our new child 0.
                        child_0_bottom_nbr_ind = bottom_nbr_ind;
                        // Only repoint if we're not crossing a periodic seam --
                        // i.e., if our panel is NOT on the bottom boundary, OR if BCs are not periodic.
                        // (If we ARE on the periodic bottom seam, panel_bottom is the wrapped panel
                        //  at the top of the domain; its top_nbr_ind should keep pointing into us
                        //  the same way it always did. The same line works either way.)
                        panel_bottom->top_nbr_ind = num_new_panels;
                    } else {
                        // Neighbor is itself y-refined: reuse its top child (index +1) as our child 0's bottom neighbor.
                        child_0_bottom_nbr_ind = panel_bottom->child_inds_start + 1;
                        Panel* child_0_bottom_nbr = &(panels[child_0_bottom_nbr_ind]);
                        child_0_bottom_nbr->top_nbr_ind = num_new_panels;
                    }
                }
            }

            // check top neighbor
            int top_nbr_ind = panel->top_nbr_ind;
            if (top_nbr_ind == -2) {
                child_1_top_nbr_ind = -2;
            } else if (top_nbr_ind == -1) {
                // Top neighbor at this level isn't refined: cascade up to the parent.
                panel_parent = &(panels[panel->parent_ind]);
                Panel* parent_top = &(panels[panel_parent->top_nbr_ind]);
                if (!(parent_top->is_refined_xy || parent_top->is_refined_y)) {
                    parent_top->needs_refinement = true;
                    need_further_refinement = true;
                    #ifdef DEBUG
                    cout << "refine_v: setting refinement flag on parent's top nbr, from panel " << jj << endl;
                    #endif
                }
                child_1_top_nbr_ind = -1;
            } else {
                Panel* panel_top = &(panels[top_nbr_ind]);
                if (panel_top->is_refined_xy) {
                    cout << "Shouldn't be allowed to call refine in v if top neighbor is refined in x and v!" << endl;
                    child_1_top_nbr_ind = top_nbr_ind;  // best-effort fallback
                } else {
                    if (!panel_top->is_refined_y) {
                        // Top neighbor exists but isn't subdivided.
                        // Make child 1's top neighbor the whole panel_top, and update panel_top's bottom_nbr_ind.
                        child_1_top_nbr_ind = top_nbr_ind;
                        panel_top->bottom_nbr_ind = num_new_panels + 1;
                    }
                    else {
                        // Top neighbor is y-refined: reuse its bottom child (index +0) as our child 1's top neighbor.
                        child_1_top_nbr_ind = panel_top->child_inds_start;
                        Panel* child_1_top_nbr = &(panels[child_1_top_nbr_ind]);
                        child_1_top_nbr->bottom_nbr_ind = num_new_panels + 1;
                    }
                }
            }

            // check right neighbor
            if (panel->right_nbr_ind == -2) {
                child_0_right_nbr_ind = -2;
                child_1_right_nbr_ind = -2;
                point_13_ind = new_vert_ind++;
                point_14_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
            } else if (panel->right_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                Panel* parent_right = &(panels[panel_parent->right_nbr_ind]);
                if (!(parent_right->is_refined_xy || parent_right->is_refined_y) ) {
                    parent_right->needs_refinement = true;
                    need_further_refinement = true;
                    // cout << "refine: setting refinement flag in panel " << jj << endl;
                }
                point_13_ind = new_vert_ind++;
                point_14_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
            } else {
                Panel* panel_right = &(panels[panel->right_nbr_ind]);
                if (! (panel_right->is_refined_xy || panel_right->is_refined_y) ) {
                    point_13_ind = new_vert_ind++;
                    point_14_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                    new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
                }
                else {
                    child_0_right_nbr_ind = panel_right->child_inds_start;
                    child_1_right_nbr_ind = panel_right->child_inds_start + 1;
                    Panel* child_0_right_nbr = &(panels[child_0_right_nbr_ind]);
                    child_0_right_nbr->left_nbr_ind = num_new_panels;
                    Panel* child_1_right_nbr = &(panels[child_1_right_nbr_ind]);
                    child_1_right_nbr->left_nbr_ind = num_new_panels + 1;
                    if (panel->is_right_bdry && bcs==periodic_bcs) {
                        point_13_ind = new_vert_ind++;
                        point_14_ind = new_vert_ind++;
                        new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                        new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
                    } else {
                        point_13_ind = child_0_right_nbr->point_inds[1];
                        point_14_ind = child_1_right_nbr->point_inds[1];
                    }
                }
            } // end check right neighbor

            // add interior points
            int point_11_ind = new_vert_ind;
            for (int ii = 0; ii < 2; ii++) {
                new_xs.push_back(subpanel_xs[2]);
                new_ys.push_back(subpanel_ys[1+2*ii]);
            }
            new_vert_ind += 2;

            // generate new panels
            // add these to list of prospective_panel_indices
            if (do_adaptive_refine) {
                for (int ii = num_new_panels; ii < num_new_panels + 4; ++ii) {
                    prospective_leaf_inds.push_back(ii);
                }
            }
            // panel->child_inds_start = num_new_panels;
            bool refining_in_y = true;
            panel->set_child_inds_start(num_new_panels, refining_in_y);
            // printf("post refinement, panel looks like:\n");
            // panel->print_panel();
            int child_level = panel->level + 1;
            int panel_ind = panel->panel_ind;
            int* point_inds = panel->point_inds;
            // for (int ii = 0; ii < ; ii++) {
            //     panel_vertex_inds[ii] = panel->point_inds[ii];
            // }
            panels.push_back(Panel {num_new_panels, child_level, panel_ind, 0, 
                    point_inds[0], point_9_ind, point_inds[1],
                    point_inds[3], point_11_ind, point_inds[4],
                    point_inds[6], point_13_ind, point_inds[7],
                    child_0_left_nbr_ind, num_new_panels + 1, 
                    child_0_right_nbr_ind, child_0_bottom_nbr_ind,
                    panel->is_left_bdry, panel->is_right_bdry,
                    false, panel->is_bottom_bdry});
            panels.push_back(Panel {num_new_panels+1, child_level, panel_ind, 1, 
                    point_inds[1], point_10_ind, point_inds[2],
                    point_inds[4], point_11_ind+1, point_inds[5],
                    point_inds[7], point_14_ind, point_inds[8],
                    child_1_left_nbr_ind, child_1_top_nbr_ind,
                    child_1_right_nbr_ind, num_new_panels,
                    panel->is_left_bdry, panel->is_right_bdry,
                    panel->is_top_bdry, false});

        } // end if panel is flagged
        
    } //end for loop through panels

    // set fs
    new_w0s.reserve(new_xs.size() );
    new_j0s.reserve(new_xs.size() );
    new_q_plus.reserve(new_xs.size() );
    new_q_minus.reserve(new_xs.size() );
    for (int ii = 0; ii < new_xs.size(); ++ii) {
        new_w0s.push_back( f(new_xs.at(ii), new_ys.at(ii)) );
        new_j0s.push_back( f(new_xs.at(ii), new_ys.at(ii)) );
        new_q_plus.push_back( f(new_xs.at(ii), new_ys.at(ii)) );
        new_q_minus.push_back( f(new_xs.at(ii), new_ys.at(ii)) );
    }

    for (int ii = 0; ii < new_xs.size(); ++ii) {
        xs.push_back(new_xs[ii]); 
        ys.push_back(new_ys[ii]); 
        w0s.push_back(new_w0s[ii]);
        j0s.push_back(new_j0s[ii]);
        q_plus.push_back(new_q_plus[ii]);
        q_minus.push_back(new_q_minus[ii]);
    }

}









void AMRStructure::refine_panels(std::function<double (double,double)> f, bool do_adaptive_refine, bool is_initial_step) {
    std::vector <double> new_xs;
    std::vector <double> new_ys;
    std::vector <double> new_w0s;
    std::vector <double> new_j0s;
    std::vector <double> new_q_plus;
    std::vector <double> new_q_minus;

    std::vector <int> prospective_leaf_inds;
    // int new_vert_ind = particles.size();
    int new_vert_ind = xs.size();
    int num_panels_before_this_iter = panels.size();
    

    // for each panel in this height 
    for (int jj = minimum_unrefined_index; jj < num_panels_before_this_iter; ++jj) {
        Panel* panel= &(panels[jj]); 
        if (panel->needs_refinement ) {
            // printf("refining panel %i\n", panel->get_panel_ind() );
            std::vector<double> panel_xs;
            std::vector<double> panel_ys;
            double dx, dy;

            const int* panel_points = panel->point_inds;
            for (int ii = 0; ii < 9; ++ii) {
                int point_ind = panel_points[ii];
                panel_xs.push_back(xs[point_ind]);
                panel_ys.push_back(ys[point_ind]);
            }
            dx = panel_xs[3] - panel_xs[0];
            dy = panel_ys[1] - panel_ys[0];
            double sub_dx = 0.5 * dx;
            double sub_dy = 0.5 * dy;

            int num_new_panels = panels.size();
            // double x_left = panel_xs[0];
            // double x_mid = x_left + .5 * dx;
            // double x_right = panel_xs[2];
            // double p_bottom = panel_ps[0];
            // double p_mid = p_bottom + .5 * dp;
            // double p_top = panel_ps[1];
            double subpanel_xs[5], subpanel_ys[5];
            for (int ii = 0; ii < 5; ii ++) {
                subpanel_xs[ii] = panel_xs[0] + sub_dx * ii;
                subpanel_ys[ii] = panel_ys[0] + sub_dy * ii;
            }

            // int left_vert_ind, bottom_vert_ind, mid_vert_ind, top_vert_ind, right_vert_ind;
            int point_9_ind, point_10_ind, point_11_ind, point_15_ind, 
                point_18_ind, point_22_ind, point_23_ind, point_24_ind;
            int child_0_bottom_nbr_ind = -1;
            int child_0_left_nbr_ind = -1;
            int child_1_left_nbr_ind = -1;
            int child_1_top_nbr_ind = -1;
            int child_2_bottom_nbr_ind = -1;
            int child_2_right_nbr_ind = -1;
            int child_3_right_nbr_ind = -1;
            int child_3_top_nbr_ind = -1;

            // generate new vertices

            Panel* panel_parent;
            // check left neighbor
            if (panel->left_nbr_ind == -2) { // no left neighbor 
                child_0_left_nbr_ind = -2;
                child_1_left_nbr_ind = -2;
                point_9_ind = new_vert_ind++;
                point_10_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
            } 
            else if (panel->left_nbr_ind == -1) { // Not refined in parent panel, need go back to parent 
                panel_parent = &(panels[panel->parent_ind]);
                Panel* parent_left = &(panels[panel_parent->left_nbr_ind]);
                if (!parent_left->is_refined_xy) {
                    parent_left->needs_refinement = true;
                    need_further_refinement = true;
                    // cout << "refine: setting refinement flag in panel " << jj << endl;
                }
                point_9_ind = new_vert_ind;
                // point_10_ind = new_vert_ind++;
                point_10_ind = point_9_ind + 1;
                new_vert_ind += 2;
                new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
            } 
            else {
                Panel* panel_left = &(panels[panel->left_nbr_ind]);
                if (! panel_left->is_refined_xy) {
                    point_9_ind = new_vert_ind++;
                    point_10_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                    new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
                }
                else {
                    child_0_left_nbr_ind = panel_left->child_inds_start +2;
                    Panel* child_0_left_nbr = &(panels[child_0_left_nbr_ind]);
                    child_0_left_nbr->right_nbr_ind = num_new_panels;
                    child_1_left_nbr_ind = panel_left->child_inds_start + 3;
                    Panel* child_1_left_nbr = &(panels[child_1_left_nbr_ind]);
                    child_1_left_nbr->right_nbr_ind = num_new_panels + 1;
                    // if (panel->is_left_bdry && bcs==0) {
                    if (panel->is_left_bdry) {
                        point_9_ind = new_vert_ind++;
                        point_10_ind = new_vert_ind++;
                        new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                        new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
                    } else {
                        point_9_ind = child_0_left_nbr->point_inds[7];
                        point_10_ind = child_1_left_nbr->point_inds[7];
                    }
                }
            }
            // check bottom neighbor
            if (panel->bottom_nbr_ind == -2) {
                child_0_bottom_nbr_ind = -2;
                child_2_bottom_nbr_ind = -2;
                point_11_ind = new_vert_ind++;
                point_18_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                new_ys.push_back(subpanel_ys[0]); new_ys.push_back(subpanel_ys[0]);
            } else if (panel->bottom_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                
                Panel* parent_bottom = &(panels[panel_parent->bottom_nbr_ind]);
                if (!parent_bottom->is_refined_xy ) {
                    parent_bottom->needs_refinement = true;
                    need_further_refinement = true;
                    #ifdef DEBUG
                    cout << "refine: setting refinement flag in panel " << jj << endl;
                    #endif
                }
                point_11_ind = new_vert_ind++;
                point_18_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                new_ys.push_back(subpanel_ys[0]); new_ys.push_back(subpanel_ys[0]);
            } else {
                Panel* panel_bottom = &(panels[panel->bottom_nbr_ind]);
                if (! panel_bottom->is_refined_xy ) {
                    point_11_ind = new_vert_ind++;
                    point_18_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                    new_ys.push_back(subpanel_ys[0]); new_ys.push_back(subpanel_ys[0]);
                }
                else {
                    child_0_bottom_nbr_ind = panel_bottom->child_inds_start + 1;
                    Panel* child_0_bottom_nbr = &(panels[child_0_bottom_nbr_ind]);
                    child_0_bottom_nbr->top_nbr_ind = num_new_panels;
                    child_2_bottom_nbr_ind = panel_bottom->child_inds_start +3;
                    Panel* child_2_bottom_nbr = &(panels[child_2_bottom_nbr_ind]);
                    child_2_bottom_nbr->top_nbr_ind = num_new_panels + 2;

                    if (panel->is_bottom_bdry && bcs == periodic_bcs) {
                        point_11_ind = new_vert_ind++;
                        point_18_ind = new_vert_ind++;
                        new_xs.push_back(subpanel_xs[1]);
                        new_xs.push_back(subpanel_xs[3]);
                        new_ys.push_back(subpanel_ys[0]);
                        new_ys.push_back(subpanel_ys[0]);
                    }
                    else {
                        point_11_ind = child_0_bottom_nbr->point_inds[5];
                        point_18_ind = child_2_bottom_nbr->point_inds[5];
                    }
                }
            }

            // check top neighbor
            if (panel->top_nbr_ind == -2) {
                child_1_top_nbr_ind = -2;
                child_3_top_nbr_ind = -2;
                point_15_ind = new_vert_ind++;
                point_22_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                new_ys.push_back(subpanel_ys[4]); new_ys.push_back(subpanel_ys[4]);
            } else if (panel->top_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                
                Panel* parent_top = &(panels[panel_parent->top_nbr_ind]);
                if (!parent_top->is_refined_xy ) {
                    parent_top->needs_refinement = true;
                    need_further_refinement = true;
                    // cout << "refine: setting refinement flag in panel " << jj << endl;
                }
                point_15_ind = new_vert_ind++;
                point_22_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                new_ys.push_back(subpanel_ys[4]); new_ys.push_back(subpanel_ys[4]);
            } else {
                Panel* panel_top = &(panels[panel->top_nbr_ind]);
                if (! panel_top->is_refined_xy ) {
                    point_15_ind = new_vert_ind++;
                    point_22_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                    new_ys.push_back(subpanel_ys[4]); new_ys.push_back(subpanel_ys[4]);
                }
                else {
                    child_1_top_nbr_ind = panel_top->child_inds_start;
                    Panel* child_1_top_nbr = &(panels[child_1_top_nbr_ind]);
                    child_1_top_nbr->bottom_nbr_ind = num_new_panels + 1;
                    child_3_top_nbr_ind = panel_top->child_inds_start + 2;
                    Panel* child_3_top_nbr = &(panels[child_3_top_nbr_ind]);
                    child_3_top_nbr->bottom_nbr_ind = num_new_panels + 3;
                    if (panel->is_top_bdry && bcs == periodic_bcs) {
                        point_15_ind = new_vert_ind++;
                        point_22_ind = new_vert_ind++;

                        new_xs.push_back(subpanel_xs[1]);
                        new_xs.push_back(subpanel_xs[3]);
                        new_ys.push_back(subpanel_ys[4]);
                        new_ys.push_back(subpanel_ys[4]);
                    }
                    else {
                        point_15_ind = child_1_top_nbr->point_inds[3];
                        point_22_ind = child_3_top_nbr->point_inds[3];
                    }
                }
            }

            // check right neighbor
            if (panel->right_nbr_ind == -2) {
                child_2_right_nbr_ind = -2;
                child_3_right_nbr_ind = -2;
                point_23_ind = new_vert_ind++;
                point_24_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
            } else if (panel->right_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                Panel* parent_right = &(panels[panel_parent->right_nbr_ind]);
                if (!parent_right->is_refined_xy ) {
                    parent_right->needs_refinement = true;
                    need_further_refinement = true;
                    // cout << "refine: setting refinement flag in panel " << jj << endl;
                }
                point_23_ind = new_vert_ind++;
                point_24_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
            } else {
                Panel* panel_right = &(panels[panel->right_nbr_ind]);
                if (! panel_right->is_refined_xy) {
                    point_23_ind = new_vert_ind++;
                    point_24_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                    new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
                    if (panel->right_nbr_ind == jj /* panel_ind = jj */) {
                        child_2_right_nbr_ind = num_new_panels;
                        child_0_left_nbr_ind = num_new_panels + 2;
                        child_3_right_nbr_ind = num_new_panels + 1;
                        child_1_left_nbr_ind = num_new_panels + 3;
                    }
                }
                else {
                    child_2_right_nbr_ind = panel_right->child_inds_start;
                    Panel* child_2_right_nbr = &(panels[child_2_right_nbr_ind]);
                    child_2_right_nbr->left_nbr_ind = num_new_panels+2;
                    child_3_right_nbr_ind = panel_right->child_inds_start + 1;
                    Panel* child_3_right_nbr = &(panels[child_3_right_nbr_ind]);
                    child_3_right_nbr->left_nbr_ind = num_new_panels + 3;

                    // if (panel->is_right_bdry and bcs==periodic_bcs) {
                    if (panel->is_right_bdry) {
                        point_23_ind = new_vert_ind++;
                        point_24_ind = new_vert_ind++;
                        new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                        new_ys.push_back(subpanel_ys[1]); new_ys.push_back(subpanel_ys[3]);
                    } else {
                        point_23_ind = child_2_right_nbr->point_inds[1];
                        point_24_ind = child_3_right_nbr->point_inds[1];
                    }
                }
            } // end check right neighbor

            // add interior points
            int point_12_ind = new_vert_ind;
            int point_16_ind = point_12_ind + 3;
            int point_19_ind = point_16_ind + 2;
            for (int ii = 0; ii < 3; ii++) { // point 12, 13, 14
                new_xs.push_back(subpanel_xs[1]);
                new_ys.push_back(subpanel_ys[1+ii]);
            }
            for (int ii = 0; ii < 2; ii++) { // point 16, 17
                new_xs.push_back(subpanel_xs[2]);
                new_ys.push_back(subpanel_ys[1+2*ii]);
            }
            for (int ii = 0; ii < 3; ii++) {// point 19, 20, 21
                new_xs.push_back(subpanel_xs[3]);
                new_ys.push_back(subpanel_ys[1+ii]);
            }
            new_vert_ind += 8;

            // generate new panels
            // add these to list of prospective_panel_indices
            if (do_adaptive_refine) {
                for (int ii = num_new_panels; ii < num_new_panels + 4; ++ii) {
                    prospective_leaf_inds.push_back(ii);
                }
            }
            // panel->child_inds_start = num_new_panels;
            panel->set_child_inds_start(num_new_panels);
            // printf("post refinement, panel looks like:\n");
            // panel->print_panel();
            int child_level = panel->level + 1;
            int panel_ind = panel->panel_ind;
            int* point_inds = panel->point_inds;
            // for (int ii = 0; ii < ; ii++) {
            //     panel_vertex_inds[ii] = panel->point_inds[ii];
            // }
            panels.push_back(Panel {num_new_panels, child_level, panel_ind, 0, 
                    point_inds[0], point_9_ind, point_inds[1],
                    point_11_ind, point_12_ind, point_12_ind + 1,
                    point_inds[3], point_16_ind, point_inds[4],
                    child_0_left_nbr_ind, num_new_panels + 1, 
                    num_new_panels + 2, child_0_bottom_nbr_ind,
                    panel->is_left_bdry, false, false, panel->is_bottom_bdry});
            panels.push_back(Panel {num_new_panels+1, child_level, panel_ind, 1, 
                    point_inds[1], point_10_ind, point_inds[2],
                    point_12_ind+1, point_12_ind+2, point_15_ind,
                    point_inds[4], point_16_ind+1, point_inds[5],
                    child_1_left_nbr_ind, child_1_top_nbr_ind,
                    num_new_panels + 3, num_new_panels,
                    panel->is_left_bdry, false, panel->is_top_bdry, false});
            panels.push_back(Panel {num_new_panels+2, child_level, panel_ind, 2, 
                    point_inds[3], point_16_ind, point_inds[4],
                    point_18_ind, point_19_ind, point_19_ind+1,
                    point_inds[6], point_23_ind, point_inds[7],
                    num_new_panels, num_new_panels+3, 
                    child_2_right_nbr_ind, child_2_bottom_nbr_ind,
                    false, panel->is_right_bdry, false, panel->is_bottom_bdry});
            panels.push_back(Panel {num_new_panels+3, child_level, panel_ind, 3,
                    point_inds[4], point_16_ind+1, point_inds[5],
                    point_19_ind+1, point_19_ind+2, point_22_ind,
                    point_inds[7], point_24_ind, point_inds[8],
                    num_new_panels+1, child_3_top_nbr_ind,
                    child_3_right_nbr_ind, num_new_panels+2,
                    false, panel->is_right_bdry, panel->is_top_bdry, false});

        } // end if panel is flagged
    } //end for loop through panels

    // set fs
    new_w0s.reserve(new_xs.size() );
    new_j0s.reserve(new_xs.size() );
    new_q_plus.reserve(new_xs.size() );
    new_q_minus.reserve(new_xs.size() );

    // refine_panels is now purely topological for values: new points get
    // placeholders here and are filled by the value pass in generate_mesh
    // (analytic ICs on the initial step; otherwise interpolate_to_initial_xys for
    // the base grid and interpolate_q_scattered after each adaptive pass). This
    // drops the per-point recursive interpolation whose result was overwritten
    // anyway, and routes every value through the corner-robust neighbor-walk.
    for (int ii = 0; ii < new_xs.size(); ++ii) {
        new_w0s.push_back(0.0);
        new_j0s.push_back(0.0);
        new_q_plus.push_back(0.0);
        new_q_minus.push_back(0.0);
    }


    for (int ii = 0; ii < new_xs.size(); ++ii) {
        xs.push_back(new_xs[ii]); 
        ys.push_back(new_ys[ii]); 
        w0s.push_back(new_w0s[ii]);
        j0s.push_back(new_j0s[ii]);
        q_plus.push_back(new_q_plus[ii]);
        q_minus.push_back(new_q_minus[ii]);
    }

}



void AMRStructure::test_panel(int panel_ind, bool verbose) {
    auto panel_it = panels.begin() + panel_ind;

    double panel_w0s[9];
    double panel_j0s[9];
    for (int ii = 0; ii < 9; ++ii) {
        panel_w0s[ii] = w0s[panel_it->point_inds[ii]];
        panel_j0s[ii] = j0s[panel_it->point_inds[ii]];
    }

    double max_w0 = panel_w0s[0], min_w0 = panel_w0s[0];
    double max_j0 = panel_j0s[0], min_j0 = panel_j0s[0];
    for (int ii = 1; ii < 9; ++ii) {
        if (panel_w0s[ii] > max_w0) max_w0 = panel_w0s[ii];
        if (panel_w0s[ii] < min_w0) min_w0 = panel_w0s[ii];
        if (panel_j0s[ii] > max_j0) max_j0 = panel_j0s[ii];
        if (panel_j0s[ii] < min_j0) min_j0 = panel_j0s[ii];
    }

    // amr_epsilons[0] is for vorticity, amr_epsilons[1] is for j
    bool refine_criteria_met =
        (max_w0 - min_w0 > amr_epsilons_vorticity) ||
        (max_j0 - min_j0 > amr_epsilons_j); 

    if (panel_it->level < max_height && refine_criteria_met) {
        panel_it->needs_refinement = true;
        need_further_refinement = true;
        if (verbose) {
            cout << "panel " << panel_ind << " is level " << panel_it->level
                 << ", max height " << max_height
                 << ", and is flagged for refinement" << endl;
        }
    } else if (verbose) {
        cout << "panel " << panel_ind << " is level " << panel_it->level
             << ", max height " << max_height
             << ", and is not flagged for refinement" << endl;
    }
}



void AMRStructure::set_leaves_weights() {
    leaf_inds = std::vector<int> ();
    weights = std::vector<double> (xs.size());
    recursively_set_leaves_weights(0);

    u_weights = std::vector<double> (weights.size(), 0.0);
    b_weights = std::vector<double> (weights.size(), 0.0);

    for (int ii = 0; ii < xs.size(); ++ii) {
        u_weights[ii] = weights[ii] * w0s[ii];
        b_weights[ii] = weights[ii] * j0s[ii];
    }
}

void AMRStructure::recursively_set_leaves_weights(int panel_ind) {
    auto panel_it = panels.begin() + panel_ind;
    if (panel_it->is_refined_y) {
        int child_start = panel_it->child_inds_start;
        for (int ii = 0; ii < 2; ii++) {
            recursively_set_leaves_weights(child_start + ii);
        }
    } else if (panel_it->is_refined_xy) {
        int child_start = panel_it->child_inds_start;
        for (int ii = 0; ii < 4; ii++) {
            recursively_set_leaves_weights(child_start + ii);
        }
    }
    else {
        leaf_inds.push_back(panel_ind);
        double dx = xs[panel_it->point_inds[3]] - xs[panel_it->point_inds[0]];
        double y0 = ys[panel_it->point_inds[0]];
        double y1 = ys[panel_it->point_inds[1]];
        double dy = y1 - y0;
        switch (quad) {
            case simpsons : {
                double dxdp9 = dx * dy / 9.0;
                double quad_weights[9] = {1.0,4.0,1.0, 4.0, 16.0,4.0,1.0,4.0,1.0};
                for (int ii = 0; ii < 9; ii++) {
                    weights[panel_it->point_inds[ii]] += dxdp9 * quad_weights[ii];
                }
                break;
            }
            default : {// trap 
                double dxdp4 = dx * dy / 4.0;
                double quad_weights[9] = {1.0,2.0,1.0,2.0,4.0,2.0,1.0, 2.0, 1.0};
                for (int ii = 0; ii < 9; ii++) {
                    weights[panel_it->point_inds[ii]] += dxdp4 * quad_weights[ii];
                }
                break;
            }
        }  
    }
}



void AMRStructure::remesh()
{
    old_panels   = panels;           
    old_xs_plus  = xs_plus;   old_ys_plus  = ys_plus;   old_q_plus  = q_plus;
    old_xs_minus = xs_minus;  old_ys_minus = ys_minus;  old_q_minus = q_minus;
 
    bool is_initial_step = false;
    generate_mesh([&](double x, double y) { return (*w0)(x, y); },
                  [&](double x, double y) { return (*j0)(x, y); },
                  do_adaptively_refine_vorticity, do_adaptively_refine_j, is_initial_step);
}