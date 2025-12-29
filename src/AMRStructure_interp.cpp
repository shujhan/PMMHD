#include "AMRStructure.hpp"

void AMRStructure::interpolate_to_initial_xys(
    std::vector<double>& fs, std::vector<double>& xs, std::vector<double>& ys, 
    int nx, int ny, bool verbose) 
{

    std::vector<double> shifted_xs(xs.size() );
    shift_xs(shifted_xs, xs, ys);
    // cout << "shift xs" << endl;
    // std::copy(shifted_xs.begin(), shifted_xs.end(), std::ostream_iterator<double>(cout, ", "));
    // cout << endl;

#ifdef DEBUG
cout << "sorting " << endl;
#endif
    std::vector<int> sort_indices(xs.size());
    for (int ii = 0; ii < xs.size(); ii++ ) { sort_indices[ii] = ii; }
    // std::iota(sort_indices.begin(), sort_indices.end(), 0);
    bool do_sort = true;
    double sort_threshold = initial_dy / 10.0;
    if (do_sort) {
        std::sort(sort_indices.begin(), sort_indices.end(),
            [&] (int a, int b) 
            { 
                if (fabs(ys[a] - ys[b]) >= sort_threshold) { return ys[a] < ys[b]; }
                else {
                    return shifted_xs[a] < shifted_xs[b];
                }
            });
    }
#ifdef DEBUG
cout << "Done sorting" << endl;
#endif
    std::vector<double> sortxs(shifted_xs.size()), sortys(ys.size());
    for (int ii = 0; ii < xs.size(); ii++) {
        sortxs[ii] = shifted_xs[sort_indices[ii]];
        sortys[ii] = ys[sort_indices[ii]];
    }


    std::vector<double> sortfs(xs.size() );

    std::vector<int> leaf_panel_of_points(xs.size() );
    std::vector<std::vector<int> > point_in_leaf_panels_by_inds(old_panels.size() );

    bool beyond_boundary = false;
    int leaf_ind = find_leaf_containing_xy_recursively(sortxs[0], sortys[0], beyond_boundary, 0, verbose);
    std::vector<int> first_column_leaf_inds(ny);




}


void AMRStructure::shift_xs(std::vector<double>& shifted_xs, const std::vector<double>& xs, const std::vector<double>& ys) {
    bool verbose = false;

    double x_bl, x_tl, x_br, x_tr;
    double y_bl, y_tl, y_br, y_tr;
    x_bl = this->old_xs[0]; y_bl = this->old_ys[0];
    x_tl = this->old_xs[2]; y_tl = this->old_ys[2];
    x_br = this->old_xs[6]; y_br = this->old_ys[6];
    x_tr = this->old_xs[8]; y_tr = this->old_ys[8];
    if (verbose) {
        cout << "(x,y)_bl (" << x_bl << ", " << y_bl << ")" << endl;
        cout << "(x,y)_tl (" << x_tl << ", " << y_tl << ")" << endl;
        cout << "(x,y)_br (" << x_br << ", " << y_br << ")" << endl;
        cout << "(x,y)_tr (" << x_tr << ", " << y_tr << ")" << endl;
    }

    for (int ii = 0; ii < xs.size(); ++ii) {
        double y = ys[ii];
        double x_temp = xs[ii];
        // trouble shooting in amr
        // if ( (fabs(x_temp) < 0.1 || fabs(x_temp -12.5664) < 0.1) && fabs(v -1.125) < 0.1) { verbose = true; }
        // else {verbose = false; }
        //end troubleshoot
        bool ineq_00_left = (x_tl - x_bl) * (y - y_bl) <= (y_tl - y_bl) * (x_temp - x_bl);


        if (verbose) {
            cout << "point " << ii << ": (x,y)= (" << x_temp << ", " << y << ")" << endl;
            cout << "ineq_00_left, " << ineq_00_left << endl; 
        }
        int counter = 0;
        while (not ineq_00_left) {
            x_temp += Lx;
            ineq_00_left = (x_tl - x_bl) * (y - y_bl) <= (y_tl - y_bl) * (x_temp - x_bl);
            if(verbose) {
                cout << "post shift x= (" << x_temp << ", ineq_00_left, " << ineq_00_left << endl; 
            }
            counter++;
            if (counter > 10) {
                throw std::runtime_error("too many shifts!");
            }
        }
        bool ineq_00_right = (x_tr - x_br) * (y - y_br) > (y_tr - y_br) * (x_temp - (x_br));
        if (verbose) {
            cout << "ineq_00_right, " << ineq_00_left << endl; 
        }
        counter = 0;
        while (not ineq_00_right) {
            x_temp -= Lx;
            ineq_00_right = (x_tr - x_br) * (y - y_br) > (y_tr - y_br) * (x_temp - (x_br));
            if(verbose) {
                cout << "post shift x= (" << x_temp << ", ineq_00_right, " << ineq_00_right << endl; 
            }
            counter++;
            if (counter > 10) {
                throw std::runtime_error("too many shifts!");
            }
        }
        shifted_xs[ii] = x_temp;
    }
}




int AMRStructure::find_leaf_containing_xy_recursively(double  &x, const double &y, bool& beyond_boundary, int panel_ind, bool verbose) {
    int leaf_ind;
    int subpanel_ind;
    int child_inds_start;

    #ifdef DEBUG
    verbose = true;
    #endif
    // double x_temp = x;
    // if ( fabs(x +1.57) < 0.5 && fabs(v +4.125) < 0.5) { verbose = true; }
    // else {verbose = false; }

    //trouble shooting
    
    Panel* panel = &(old_panels[panel_ind]);
    child_inds_start = panel->child_inds_start;

    if (verbose) {
        cout << "In panel " << panel_ind << endl;
        cout << *panel << endl;
        cout << "testing (x,y)=(" << x << ", " << y << ")" << endl;
    }
    if (! (panel->is_refined_xy || panel->is_refined_y ) ) {
        leaf_ind = panel_ind;
        if (verbose) {
            cout << "leaf panel!" << endl;
        }

        if (!allow_boundary_extrapolation) {
            double x_bl = old_xs[panel->point_inds[0]]; double y_bl = old_ys[panel->point_inds[0]];
            double x_tl = old_xs[panel->point_inds[2]]; double y_tl = old_ys[panel->point_inds[2]];
            double x_br = old_xs[panel->point_inds[6]]; double y_br = old_ys[panel->point_inds[6]];
            double x_tr = old_xs[panel->point_inds[8]]; double y_tr = old_ys[panel->point_inds[8]];
            bool ineq_right = (x_tr - x_br) * (y - y_br) > (y_tr - y_br) * (x - (x_br));
            bool ineq_left = (x_tl - x_bl) * (y - y_bl) <= (y_tl - y_bl) * (x - x_bl);
            bool ineq_top = (x_tr - x_tl) * (y - y_tl) < (y_tr - y_tl) * (x - x_tl);
            bool ineq_bottom = (x_br - x_bl) * (y - y_bl) >= (y_br - y_bl) * (x - x_bl);
            bool boundary_extrapolating_right = !ineq_right && panel->right_nbr_ind==-2;
            bool boundary_extrapolating_left = !ineq_left && panel->left_nbr_ind==-2;
            bool boundary_extrapolating_top = !ineq_top && panel->top_nbr_ind==-2;
            bool boundary_extrapolating_bottom = !ineq_bottom && panel->bottom_nbr_ind==-2;
            if (boundary_extrapolating_left || boundary_extrapolating_right || 
                boundary_extrapolating_top || boundary_extrapolating_bottom) {
                if (verbose) {
                    cout << "Not allowing boundary extrapolation and this point is flagged for beyond the boundary" << endl;
                }
                beyond_boundary = true;
            }
        }
    } else {
        double x_bl = old_xs[panel->point_inds[0]]; double y_bl = old_ys[panel->point_inds[0]];
        double x_ml = old_xs[panel->point_inds[1]]; double y_ml = old_ys[panel->point_inds[1]];
        double x_tl = old_xs[panel->point_inds[2]]; double y_tl = old_ys[panel->point_inds[2]];
        double x_bm = old_xs[panel->point_inds[3]]; double y_bm = old_ys[panel->point_inds[3]];
        double x_mm = old_xs[panel->point_inds[4]]; double y_mm = old_ys[panel->point_inds[4]];
        double x_tm = old_xs[panel->point_inds[5]]; double y_tm = old_ys[panel->point_inds[5]];
        double x_br = old_xs[panel->point_inds[6]]; double y_br = old_ys[panel->point_inds[6]];
        double x_mr = old_xs[panel->point_inds[7]]; double y_mr = old_ys[panel->point_inds[7]];
        double x_tr = old_xs[panel->point_inds[8]]; double y_tr = old_ys[panel->point_inds[8]];
        // left_shear_ineq = (x_tl - x_bl) * (v - v_bl) > (v_tl - v_bl) * (x_temp - x_bl)

        bool ineq_1_bottom = (x_mm - x_ml) * (y - y_ml) >= (y_mm - y_ml) * (x - x_ml);
        bool ineq_1_right = (x_tm - x_mm) * (y - y_mm) >= (y_tm - y_mm) * (x - x_mm);
        bool ineq_3_bottom = (x_mr - x_mm) * (y - y_mm) >= (y_mr - y_mm) * (x - x_mm);

        if (verbose) {
            std::cout << "ineq_1_bottom = " << ineq_1_bottom << endl;
            std::cout << "ineq_1_right = " << ineq_1_right << endl;
            std::cout << "ineq_3_bottom = " << ineq_3_bottom << endl;
        }

        if (ineq_1_bottom && ineq_1_right) {

            bool ineq_1_top = (x_tm - x_tl) * (y - y_tl) <= (y_tm - y_tl) * (x - x_tl);
            Panel* child_1 = &old_panels[child_inds_start+1];
            int child_1_top_nbr_ind = child_1->top_nbr_ind;
            if (ineq_1_top ||  child_1_top_nbr_ind < 0) {
                bool ineq_1_left = (x_tl - x_ml) * (y - y_ml) <= (y_tl - y_ml) * (x - x_ml);
                int child_1_left_nbr_ind = child_1->left_nbr_ind;
                if (ineq_1_left || child_1_left_nbr_ind < 0) {
                    subpanel_ind = child_inds_start + 1;
                    if (verbose) {
                        cout << "in child 1, panel " << subpanel_ind << endl;
                    }
                }
                else {
                    subpanel_ind = child_1_left_nbr_ind;
                    if (panel->is_left_bdry and bcs==0) {
                        x += Lx;
                    }
                }
            } else {
                subpanel_ind = child_1_top_nbr_ind;
            }
            
        } else if (ineq_3_bottom && !ineq_1_right)
        {
            bool ineq_3_top = (x_tr - x_tm) * (y - y_tm) <= (y_tr - y_tm) * (x - x_tm);
            Panel* child_3;
            if (panel->is_refined_y) {
                child_3 = &old_panels[child_inds_start +1];
            } else { // panel is refined in xv
                child_3 = &old_panels[child_inds_start+3];
            }
            int child_3_top_nbr_ind = child_3->top_nbr_ind;
            if (ineq_3_top || child_3_top_nbr_ind < 0) {
                bool ineq_3_right = (x_tr - x_mr) * (y - y_mr) <= (y_tr - y_mr) * (x - x_mr);
                int child_3_right_nbr_ind = child_3->right_nbr_ind;
                if (!ineq_3_right || child_3_right_nbr_ind < 0) {
                    if (panel->is_refined_y) {
                        subpanel_ind = child_inds_start + 1;
                    }
                    else {
                        subpanel_ind = child_inds_start + 3;
                    }
                    if (verbose) {
                        cout << "in child 3, panel " << subpanel_ind << endl;
                    }
                }
                else {
                    subpanel_ind = child_3_right_nbr_ind;
                    if (panel->is_right_bdry && bcs==0) {
                        x -= Lx;
                    }
                }
            } else {
                subpanel_ind = child_3_top_nbr_ind;
            }
        } else {
            bool ineq_0_right = (x_mm - x_bm) * (y - y_bm) >= (y_mm - y_bm) * (x - x_bm);
            if (verbose) {
                std::cout << "ineq_0_right = " << ineq_0_right << endl;
            }
            if (ineq_0_right) {
                bool ineq_0_bottom = (x_bm - x_bl) * (y - y_bl) >= (y_bm - y_bl) * (x - x_bl);
                Panel* child_0 = &old_panels[child_inds_start];
                int child_0_bottom_nbr_ind = child_0->bottom_nbr_ind;
                if (ineq_0_bottom || child_0_bottom_nbr_ind < 0) {
                    bool ineq_0_left = (x_ml - x_bl) * (y - y_bl) <= (y_ml - y_bl) * (x - x_bl);
                    int child_0_left_nbr_ind = child_0->left_nbr_ind;
                    if (ineq_0_left || child_0_left_nbr_ind < 0) {
                        
                        subpanel_ind = child_inds_start;
                        if (verbose) {
                            cout << "in child 0, panel " << subpanel_ind << endl;
                        }
                    }
                    else {
                        subpanel_ind = child_0_left_nbr_ind;
                        if (panel->is_left_bdry && bcs==0) 
                        { x += Lx; }
                    }
                } else {
                    subpanel_ind = child_0_bottom_nbr_ind;
                }
            } else
            {
                bool ineq_2_bottom = (x_br - x_bm) * (y - y_bm) >= (y_br - y_bm) * (x - x_bm);
                Panel* child_2;
                if (panel->is_refined_y) {
                     child_2 = &old_panels[child_inds_start];
                } else { // panel is refined in x and y
                     child_2 = &old_panels[child_inds_start+2];
                }
                int child_2_bottom_nbr_ind = child_2->bottom_nbr_ind;
                if (ineq_2_bottom || child_2_bottom_nbr_ind < 0) {
                    bool ineq_2_right = (x_mr - x_br) * (y - y_br) <= (y_mr - y_br) * (x - x_br);
                    int child_2_right_nbr_ind = child_2->right_nbr_ind;
                    if (! ineq_2_right || child_2_right_nbr_ind < 0) {
                        if (panel->is_refined_y) {
                            subpanel_ind = child_inds_start;
                        } else {
                            subpanel_ind = child_inds_start + 2;
                        }
                        if (verbose) {
                            cout << "in child 2, panel " << subpanel_ind << endl;
                        }
                    }
                    else {
                        subpanel_ind = child_2_right_nbr_ind;
                        if (panel->is_right_bdry && bcs == 0) {
                            x -= Lx;
                        }
                    }
                } else {
                    subpanel_ind = child_2_bottom_nbr_ind;
                }
            }
            
        }

        leaf_ind = find_leaf_containing_xy_recursively(x, y, beyond_boundary, subpanel_ind, verbose);

    }
    return leaf_ind;
}

