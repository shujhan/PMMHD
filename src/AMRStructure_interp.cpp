#include "AMRStructure.hpp"

void AMRStructure::interpolate_to_initial_xys(
    std::vector<double>& q0s, std::vector<double>& xs, std::vector<double>& ys, 
    int nx, int ny) 
{
    std::vector<double> shifted_xs(xs.size());
    shift_xs(shifted_xs, xs, ys); // shift xs to the bounded domain 
    std::vector<double> shifted_ys(ys.size());
    for (int ii = 0; ii < ys.size(); ++ii) {
        shifted_ys[ii] = ys[ii];
    }
    // also shift y back into the principal periodic strip when y is periodic
    if (bcs == periodic_bcs) {
        double x_bl = this->old_xs[0], y_bl = this->old_ys[0];
        double x_tl = this->old_xs[2], y_tl = this->old_ys[2];
        double x_br = this->old_xs[6], y_br = this->old_ys[6];
        double x_tr = this->old_xs[8], y_tr = this->old_ys[8];

        for (int ii = 0; ii < shifted_ys.size(); ++ii) {
            double x = shifted_xs[ii];
            double y_temp = shifted_ys[ii];

            bool ineq_bottom = (x_br - x_bl) * (y_temp - y_bl) >= (y_br - y_bl) * (x - x_bl);
            int counter = 0;
            while (!ineq_bottom) {
                y_temp += Ly;
                ineq_bottom = (x_br - x_bl) * (y_temp - y_bl) >= (y_br - y_bl) * (x - x_bl);
                counter++;
                if (counter > 10) {
                    throw std::runtime_error("too many y shifts at bottom!");
                }
            }

            bool ineq_top = (x_tr - x_tl) * (y_temp - y_tl) <= (y_tr - y_tl) * (x - x_tl);
            counter = 0;
            while (!ineq_top) {
                y_temp -= Ly;
                ineq_top = (x_tr - x_tl) * (y_temp - y_tl) <= (y_tr - y_tl) * (x - x_tl);
                counter++;
                if (counter > 10) {
                    throw std::runtime_error("too many y shifts at top!");
                }
            }

            shifted_ys[ii] = y_temp;
        }
    }

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
                if (fabs(shifted_ys[a] - shifted_ys[b]) >= sort_threshold) {
                    return shifted_ys[a] < shifted_ys[b];
                }
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
        sortys[ii] = shifted_ys[sort_indices[ii]];
    }


    // std::vector<double> sortw0s(xs.size());
    // std::vector<double> sortj0s(xs.size());
    std::vector<double> sortq0s(xs.size());

    std::vector<int> leaf_panel_of_points(xs.size() );
    std::vector<std::vector<int> > point_in_leaf_panels_by_inds(old_panels.size() );

    bool beyond_boundary = false;
    // find the leaf index for the first point first 
    int leaf_ind = find_leaf_containing_xy_recursively(sortxs[0], sortys[0], beyond_boundary, 0);
    std::vector<int> first_column_leaf_inds(ny);

    // utilize periodic bc in x, find leaf index for the first point in each row from the previous first point 
    for (int ii = 0; ii < ny; ++ii) {
        beyond_boundary = false;
        int point_ind = ii * nx;
        std::set<int> history;
        history.emplace(leaf_ind);

        leaf_ind = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortys[point_ind], beyond_boundary, 
                                                                leaf_ind, history);
        first_column_leaf_inds[ii] = leaf_ind;
        if (beyond_boundary) {
            leaf_panel_of_points[point_ind] = 0;
        } 
        else {
            leaf_panel_of_points[point_ind] = leaf_ind;
        }
    }

    // find leaf index for each row, and since we know all leaf index for the first point in each row, we can use neighbor to speed up searching 
    for (int ii = 0; ii < ny; ++ii) {
        int jj0 = 0;
        int point_ind = ii * nx + jj0;
        int leaf_ind_c = first_column_leaf_inds[ii];
        for (int jj = jj0+1; jj < nx; ++jj) {
            beyond_boundary = false;
            point_ind++;
            std::set<int> history;
            history.emplace(leaf_ind_c);
            leaf_ind_c = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortys[point_ind], beyond_boundary, leaf_ind_c, history);
            if (beyond_boundary) {
                leaf_panel_of_points[point_ind] = 0;
            } else {
                leaf_panel_of_points[point_ind] = leaf_ind_c;
            }
        }
    }    


    // Repair neighbor-walk failures before depositing points into panels.
    // find_leaf_containing_point_from_neighbor can terminate on a cycle (history
    // hit, line ~655) and return a leaf that does NOT contain the target point.
    // This happens where the source mesh is sheared/folded (current sheets) or
    // mixed-level (coarse-fine interface under AMR). That point would then be
    // interpolated from the wrong panel and evaluated far outside its node
    // range -> an isolated outlier; RK4 feeds those back into the stage slopes,
    // so the specks compound step over step. For any point not inside its
    // assigned leaf, redo the search with a recursive descent from the root
    // (the same method that seeds the first point), which always returns a
    // containing leaf (or flags beyond_boundary).
    auto point_in_old_leaf = [&](double tx, double ty, int pl) -> bool {
        const Panel& P = old_panels[pl];
        double x_bl=old_xs[P.point_inds[0]], y_bl=old_ys[P.point_inds[0]];
        double x_tl=old_xs[P.point_inds[2]], y_tl=old_ys[P.point_inds[2]];
        double x_mid=old_xs[P.point_inds[4]], y_mid=old_ys[P.point_inds[4]];
        double x_br=old_xs[P.point_inds[6]], y_br=old_ys[P.point_inds[6]];
        double x_tr=old_xs[P.point_inds[8]], y_tr=old_ys[P.point_inds[8]];
        if (tx - x_mid >= Lx/2) tx -= Lx;
        if (tx - x_mid < -Lx/2) tx += Lx;
        if (bcs == periodic_bcs) {
            if (ty - y_mid >= Ly/2) ty -= Ly;
            if (ty - y_mid < -Ly/2) ty += Ly;
        }
        bool r =(x_tr-x_br)*(ty-y_br) >= (y_tr-y_br)*(tx-x_br);
        bool l =(x_tl-x_bl)*(ty-y_bl) <= (y_tl-y_bl)*(tx-x_bl);
        bool tp=(x_tr-x_tl)*(ty-y_tl) <= (y_tr-y_tl)*(tx-x_tl);
        bool bt=(x_br-x_bl)*(ty-y_bl) >= (y_br-y_bl)*(tx-x_bl);
        return r && l && tp && bt;
    };
    for (int ii = 0; ii < (int)leaf_panel_of_points.size(); ++ii) {
        // pl<=0 is a point intentionally flagged beyond the boundary; leave it.
        if (leaf_panel_of_points[ii] <= 0) continue;
        if (!point_in_old_leaf(sortxs[ii], sortys[ii], leaf_panel_of_points[ii])) {
            bool bb = false;
            double rx = sortxs[ii], ry = sortys[ii];
            int rl = find_leaf_containing_xy_recursively(rx, ry, bb, 0);
            if (!bb && rl > 0 && rl < (int)old_panels.size()) {
                leaf_panel_of_points[ii] = rl;
            }
        }
    }




    for (int ii = 0; ii < leaf_panel_of_points.size(); ++ii) {
        point_in_leaf_panels_by_inds[leaf_panel_of_points[ii]].push_back(ii);
    }

    // eval interpolant 

    for (int panel_ind = 0; panel_ind < old_panels.size(); panel_ind++) {
        if (point_in_leaf_panels_by_inds[panel_ind].size() > 0) {
            interpolate_from_panel_to_points(sortq0s,sortxs,sortys,point_in_leaf_panels_by_inds[panel_ind], panel_ind, use_limiter, limit_val);
        }
    }
    for (int ii = 0; ii < q0s.size(); ii++) {
        q0s[sort_indices[ii]] = sortq0s[ii];
    }




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





int AMRStructure::find_leaf_containing_xy_recursively(double  &x, double &y, bool& beyond_boundary, int panel_ind) {
    int leaf_ind;
    int subpanel_ind;
    int child_inds_start;
    bool verbose = false;

    #ifdef DEBUG
    verbose = true;
    #endif
    // double x_temp = x;
    // if ( fabs(x +1.57) < 0.5 && fabs(v +4.125) < 0.5) { verbose = true; }
    // else {verbose = false; }

    //trouble shooting
    
    Panel* panel = &(old_panels[panel_ind]);
    child_inds_start = panel->child_inds_start;


    if (! (panel->is_refined_xy || panel->is_refined_y ) ) {
        leaf_ind = panel_ind;
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

        // if (verbose) {
        //     std::cout << "ineq_1_bottom = " << ineq_1_bottom << endl;
        //     std::cout << "ineq_1_right = " << ineq_1_right << endl;
        //     std::cout << "ineq_3_bottom = " << ineq_3_bottom << endl;
        // }

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
                    if (panel->is_left_bdry) {
                        x += Lx;
                    }
                }
            } else {
                subpanel_ind = child_1_top_nbr_ind;
                if (panel->is_top_bdry && bcs == periodic_bcs) {
                    y -= Ly;
                }
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
                    // if (verbose) {
                    //     cout << "in child 3, panel " << subpanel_ind << endl;
                    // }
                }
                else {
                    subpanel_ind = child_3_right_nbr_ind;
                    if (panel->is_right_bdry) {
                        x -= Lx;
                    }
                }
            } else {
                subpanel_ind = child_3_top_nbr_ind;
                if (panel->is_top_bdry && bcs == periodic_bcs) {
                    y -= Ly;
                }
            }
        } else {
            bool ineq_0_right = (x_mm - x_bm) * (y - y_bm) >= (y_mm - y_bm) * (x - x_bm);
            // if (verbose) {
            //     std::cout << "ineq_0_right = " << ineq_0_right << endl;
            // }
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
                        if (panel->is_left_bdry) 
                        { x += Lx; }
                    }
                } else {
                    subpanel_ind = child_0_bottom_nbr_ind;
                    if (panel->is_bottom_bdry && bcs == periodic_bcs) {
                        y += Ly;
                    }
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
                        // if (verbose) {
                        //     cout << "in child 2, panel " << subpanel_ind << endl;
                        // }
                    }
                    else {
                        subpanel_ind = child_2_right_nbr_ind;
                        if (panel->is_right_bdry) {
                            x -= Lx;
                        }
                    }
                } else {
                    subpanel_ind = child_2_bottom_nbr_ind;
                    if (panel->is_bottom_bdry && bcs == periodic_bcs) {
                        y += Ly;
                    }
                }
            }
            
        }

        leaf_ind = find_leaf_containing_xy_recursively(x, y, beyond_boundary, subpanel_ind);

    }
    return leaf_ind;
}



int AMRStructure::find_leaf_containing_point_from_neighbor(double& tx, double& ty, bool& beyond_boundary, int leaf_ind, std::set<int>& history) {

    bool verbose = false;
    if (leaf_ind == 0) {
        history.emplace(leaf_ind);
        cout << "If you see this then you are sending in a leaf ind 0 somewhere you hoped not to." << endl;
        return leaf_ind;
    } else {
        Panel* panel = &(old_panels[leaf_ind]);
        double x_bl = old_xs[panel->point_inds[0]]; double y_bl = old_ys[panel->point_inds[0]];
        double x_tl = old_xs[panel->point_inds[2]]; double y_tl = old_ys[panel->point_inds[2]];
        double x_mid = old_xs[panel->point_inds[4]];
        double y_mid = old_ys[panel->point_inds[4]];
        double x_br = old_xs[panel->point_inds[6]]; double y_br = old_ys[panel->point_inds[6]];
        double x_tr = old_xs[panel->point_inds[8]]; double y_tr = old_ys[panel->point_inds[8]];
        // if (verbose) {
        //     cout << "(x,y)_bl = (" << x_bl << ", " << y_bl << ")" << endl;
        //     cout << "(x,y)_tl = (" << x_tl << ", " << y_tl << ")" << endl;
        //     cout << "(x,y)_br = (" << x_br << ", " << y_br << ")" << endl;
        //     cout << "(x,y)_tr = (" << x_tr << ", " << y_tr << ")" << endl;
        // }
        // need to correct periodic distance
        if (tx - x_mid >= Lx/2) { 
            tx -= Lx; 
            // if (verbose) {
            //     cout << "shifting across boundary, tx= " << tx << endl;
            // }
        }
        if (tx - x_mid < -Lx/2) { 
            tx += Lx; 
            // if (verbose) {
            //     cout << "shifting across boundary, tx= " << tx << endl;
            // }
        }

        // periodic correction in y
        if (ty - y_mid >= Ly/2 && bcs == periodic_bcs) {
            ty -= Ly;
        }
        if (ty - y_mid < -Ly/2 && bcs == periodic_bcs) {
            ty += Ly;
        }


        bool ineq_right = (x_tr - x_br) * (ty - y_br) >= (y_tr - y_br) * (tx - (x_br));
        bool ineq_left = (x_tl - x_bl) * (ty - y_bl) <= (y_tl - y_bl) * (tx - x_bl);
        bool ineq_top = (x_tr - x_tl) * (ty - y_tl) <= (y_tr - y_tl) * (tx - x_tl);
        bool ineq_bottom = (x_br - x_bl) * (ty - y_bl) >= (y_br - y_bl) * (tx - x_bl);
        int new_leaf_ind = leaf_ind;
        // if (verbose) {
        //     cout << "(tx,ty) = (" << tx << ", " << ty << ")" << endl;
        //     cout << "testing leaf panel " << leaf_ind << " for containment" << endl;
        //     cout << "ineq_right " << ineq_right << endl;
        //     cout << "(x_tr - x_tl) * (tv - v_tl)" << (x_tr - x_tl) * (ty - y_tl) << endl;
        //     cout << "(v_tr - v_tl) * (tx - x_tl)" << (y_tr - y_tl) * (tx - x_tl) << endl;
        //     cout << "ineq_top " << ineq_top << endl;
        //     cout << "ineq_bottom " << ineq_bottom << endl;
        //     cout << "ineq_left " << ineq_left << endl;
        // }
        if (! ineq_right && panel->right_nbr_ind != -2) {
            if (panel->is_right_bdry) { 
                tx -= Lx; 
                if (verbose) {
                    cout << "shifting across boundary, tx= " << tx << endl;
                }
            }

            if (panel->right_nbr_ind == -1) {
                new_leaf_ind = old_panels[panel->parent_ind].right_nbr_ind;
                if (verbose) {
                    cout << "in parent right" << endl;
                    cout << "next leaf test " << new_leaf_ind << endl;
                }
            } else {
                Panel* panel_right = &old_panels[panel->right_nbr_ind];
                if (panel_right->is_refined_xy) { 
                    new_leaf_ind = panel_right->child_inds_start; 
                    if (verbose) {
                        cout << "in panel right children" << endl;
                        cout << "next leaf test " << new_leaf_ind << endl;
                    }
                }
                else { 
                    new_leaf_ind = panel->right_nbr_ind; 
                    // if (verbose) {
                    //     cout << "in panel right" << endl;
                    //     cout << "next leaf test " << new_leaf_ind << endl;
                    // }
                }
            }
        } else {
            if (!ineq_top && panel->top_nbr_ind != -2) {
                if (panel->is_top_bdry && bcs == periodic_bcs) {
                    ty -= Ly;
                }
                if (panel->top_nbr_ind == -1) {
                    new_leaf_ind = old_panels[panel->parent_ind].top_nbr_ind;
                    // if (verbose) {
                    //     cout << "in parent top" << endl;
                    //     cout << "next leaf test " << new_leaf_ind << endl;
                    // }
                } else {
                    Panel* panel_top = &old_panels[panel->top_nbr_ind];
                    if (panel_top->is_refined_xy) {
                        new_leaf_ind = panel_top->child_inds_start;
                        // if (verbose) {
                        //     cout << "in top children" << endl;
                        //     cout << "next leaf test " << new_leaf_ind << endl;
                        // }
                    } else {
                        new_leaf_ind = panel->top_nbr_ind;
                        // if (verbose) {
                        //     cout << "in panel top" << endl;
                        //     cout << "next leaf test " << new_leaf_ind << endl;
                        // }
                    }
                }
            } // end if checking top boundary 
            else {
                if (!ineq_bottom && panel->bottom_nbr_ind != -2) {
                    if (panel->is_bottom_bdry && bcs == periodic_bcs) {
                        ty += Ly;
                    }
                    if (panel->bottom_nbr_ind == -1) {
                        new_leaf_ind = old_panels[panel->parent_ind].bottom_nbr_ind;
                        // if (verbose) {
                        //     cout << "in parent bottom" << endl;
                        //     cout << "next leaf test " << new_leaf_ind << endl;
                        // }
                    } else {
                        Panel* panel_bottom = &old_panels[panel->bottom_nbr_ind];
                        if (panel_bottom -> is_refined_xy) {
                            new_leaf_ind = panel_bottom->child_inds_start + 1;
                            // if (verbose) {
                            //     cout << "in parent bottom children" << endl;
                            //     cout << "next leaf test " << new_leaf_ind << endl;
                            // }
                        } else {
                            new_leaf_ind = panel->bottom_nbr_ind;
                            // if (verbose) {
                            //     cout << "in parent bottom" << endl;
                            //     cout << "next leaf test " << new_leaf_ind << endl;
                            // }
                        }
                    }
                } // end if checking bottom boundary 
                else {
                    if (! ineq_left && panel->left_nbr_ind != -2) {
                        if (panel->is_left_bdry) { 
                            tx += Lx; 
                            // if (verbose) {
                            //     cout << "shifting across boundary, now tx= " << tx << endl;
                            // }
                        }
                        if (panel->left_nbr_ind == -1) {
#ifdef DEBUG
cout << "panel left: " << panel->left_nbr_ind << endl;
cout <<"length of panels_list " << old_panels.size() << endl;
#endif
                            new_leaf_ind = old_panels[panel->parent_ind].left_nbr_ind;
                            // if (verbose) {
                            //     cout << "in parent left" << endl;
                            //     cout << "next leaf test " << new_leaf_ind << endl;
                            // }
                        } else {
                            Panel* panel_left = &old_panels[panel->left_nbr_ind];
                            if (panel_left->is_refined_xy) {
                                new_leaf_ind = panel_left->child_inds_start+2;
                                // if (verbose) {
                                //     cout << "in left children" << endl;
                                //     cout << "next leaf test " << new_leaf_ind << endl;
                                // }
                            }
                            else {
                                new_leaf_ind = panel->left_nbr_ind;
                                // if (verbose) {
                                //     cout << "in panel left" << endl;
                                //     cout << "next leaf test " << new_leaf_ind << endl;
                                // }
                            }
                        }
                    } // end if checking left boundary
                } // end else for left
            } //end else for bottom/left
        } // end else for top/bottom/left
        // if (verbose) {
        //     cout << "History" << endl;
        //     for (std::set<int>::iterator it = history.begin(); it != history.end(); ++it) {
        //         cout << *it << " ";
        //     }
        //     cout << endl;
        // }
        if (history.find(new_leaf_ind) == history.end()) {
            history.emplace(new_leaf_ind);
            // if (verbose) {
            //     cout << "running find_leaf again for (x,y)=(" << tx << ", " << ty << ")" << endl;
            // }
            new_leaf_ind = find_leaf_containing_point_from_neighbor(tx,ty,beyond_boundary, new_leaf_ind, history);
        } 
        // else if (verbose)
        // {
        //     cout << "done searching." << endl;
        // }
        if (!allow_boundary_extrapolation) {
            bool boundary_extrapolating_right = !ineq_right && panel->right_nbr_ind==-2;
            bool boundary_extrapolating_left = !ineq_left && panel->left_nbr_ind==-2;
            bool boundary_extrapolating_top = !ineq_top && panel->top_nbr_ind==-2;
            bool boundary_extrapolating_bottom = !ineq_bottom && panel->bottom_nbr_ind==-2;
            if (boundary_extrapolating_left || boundary_extrapolating_right || 
                boundary_extrapolating_top || boundary_extrapolating_bottom) {
                // if (verbose) {
                //     cout << "Not allowing boundary extrapolation and this point is flagged as beyond the boundary" << endl;
                // }
                beyond_boundary = true;
            }
        }
        // if (verbose) {
        //     cout << "Point (" << tx << ", " << ty << ") is in panel " << new_leaf_ind << endl;
        // }
        return new_leaf_ind;
    }
}


void AMRStructure::interpolate_from_panel_to_points(
    std::vector<double>& values_q0,std::vector<double>& xs, std::vector<double>& ys,
    std::vector<int>& point_inds, int panel_ind, bool use_limiter, double limit_val) 
{
    if (panel_ind == 0) { // if we are extrapolating beyond boundaries; assume 0
        for (int ii = 0; ii < point_inds.size(); ++ii) {
            // values_w0[point_inds[ii]] = f_beyond_boundary;
            // values_j0[point_inds[ii]] = f_beyond_boundary;
            values_q0[point_inds[ii]] = q0_beyond_boundary;
        }
    }
    else {
        Panel* panel = &(old_panels[panel_ind]);
        const int* panel_point_inds = panel->point_inds;
        double panel_xs[9], panel_ys[9];
        double panel_q0s[9];

        for (int ii = 0; ii < 9; ++ii) {
            int pind = panel_point_inds[ii];
            panel_xs[ii] = old_xs[pind];
            panel_ys[ii] = old_ys[pind];
            panel_q0s[ii] = old_q0s[pind];
        }

        // if (do_unshear) {
        //     for (int ii = 0; ii < 9; ++ii) {
        //         panel_xs[ii] -= dt * (panel_ys[ii] - panel_ys[4]); // assumes remesh frequency = 1
        //     }
        //     for (int ii = 0; ii < xs.size(); ++ii) {
        //         xs[ii] -= dt * (ys[ii] - panel_ys[4]);

        //     }
        // }

        double panel_dx[9], panel_dy[9];
        std::vector<double> dxs(point_inds.size()), dys(point_inds.size());
        for (int ii = 0; ii < point_inds.size(); ++ii) {
            int pind = point_inds[ii];
            dxs[ii] = xs[pind] - panel_xs[4];
            dys[ii] = ys[pind] - panel_ys[4];
        }
        for (int ii = 0; ii < 9; ii ++) {
            panel_dx[ii] = panel_xs[ii] - panel_xs[4];
            panel_dy[ii] = panel_ys[ii] - panel_ys[4];
        }

        // if (verbose) {
        //     std::cout << "test point distance from midpoint:" << std::endl;
        //     for (int ii = 0; ii < point_inds.size(); ++ii) {
        //         std::cout << ii <<": dx=" << dxs[ii] <<", dy=" << dys[ii] << std::endl;
        //     } 
        //     std::cout << "panel vertex distances from midpoint:" << std::endl;
        //     for (int ii = 0; ii < 9; ii++ ) {
        //         std::cout << ii << ": " << panel_dx[ii] << ", " << panel_dy[ii] << std::endl;
        //     }
        // }

        Eigen::Matrix<double,9,9> A;
        for (int ii = 0; ii < 9; ++ii) {
            A(ii,0) = 1; A(ii,1) = panel_dx[ii];
            A(ii,2) = panel_dx[ii] * panel_dy[ii];
            A(ii,3) = panel_dy[ii];
            A(ii,4) = panel_dx[ii] * panel_dx[ii];
            A(ii,5) =  panel_dx[ii] * panel_dx[ii] * panel_dy[ii];
            A(ii,6) = panel_dx[ii] * panel_dx[ii] * panel_dy[ii] * panel_dy[ii];
            A(ii,7) = panel_dx[ii] * panel_dy[ii] * panel_dy[ii];
            A(ii,8) = panel_dy[ii] * panel_dy[ii];
        }
        // create RHS vector for q0
        Eigen::Map<Eigen::Matrix<double,9,1>> f_q0(panel_q0s);
        // solve for q0
        Eigen::Matrix<double,9,1> c_q0 = A.lu().solve(f_q0);

        // if (verbose) {
        //     std::cout << "Here is the matrix A:\n" << A << std::endl;
        //     std::cout << "Here is the f vector b:\n" << b << std::endl;
        //     std::cout << "Here is the coefficient vector c:\n" << c << std::endl;
        // }

        Eigen::Matrix<double, Dynamic, Dynamic> Dx(point_inds.size(),9);
        for (int ii = 0; ii < point_inds.size(); ++ii) {
            double dxsq = dxs[ii] * dxs[ii];
            double dysq = dys[ii] * dys[ii];
            Dx(ii,0) = 1; Dx(ii,1) = dxs[ii];
            Dx(ii,2) = dxs[ii] * dys[ii];
            Dx(ii,3) = dys[ii];
            Dx(ii,4) = dxsq;
            Dx(ii,5) =  dxsq * dys[ii];
            Dx(ii,6) = dxsq * dysq;
            Dx(ii,7) = dxs[ii] * dysq;
            Dx(ii,8) = dysq;
        }

        Eigen::Matrix<double, Dynamic,1> interp_vals_q0 = Dx * c_q0;

        // Extrapolation backstop. A target that lands outside this panel's node
        // extent is being extrapolated, and the biquadratic's high-order terms
        // (dx^2 dy^2, ...) can diverge there -> an isolated outlier (the late-
        // time salt-and-pepper). Bound any extrapolated value to the range
        // spanned by the panel's 9 nodes (widened by extrap_margin*range so a
        // legitimate gradient continuing just past the panel edge is not
        // flattened). Interior targets are left at full accuracy. This cannot
        // create a new extremum from outside the panel, so it cannot inject the
        // spurious energy seen at t~1.0.
        const double extrap_margin = 0.5;
        double hx = 0.0, hy = 0.0, qlo = panel_q0s[0], qhi = panel_q0s[0];
        for (int ii = 0; ii < 9; ++ii) {
            if (fabs(panel_dx[ii]) > hx) hx = fabs(panel_dx[ii]);
            if (fabs(panel_dy[ii]) > hy) hy = fabs(panel_dy[ii]);
            if (panel_q0s[ii] < qlo) qlo = panel_q0s[ii];
            if (panel_q0s[ii] > qhi) qhi = panel_q0s[ii];
        }
        double pad = extrap_margin * (qhi - qlo);
        double qlo_b = qlo - pad, qhi_b = qhi + pad;

        for (int ii = 0; ii < point_inds.size(); ++ii) {
            double v = interp_vals_q0(ii);
            if (fabs(dxs[ii]) > hx || fabs(dys[ii]) > hy) {   // extrapolating
                if (v < qlo_b) v = qlo_b;
                if (v > qhi_b) v = qhi_b;
            }
            values_q0[point_inds[ii]] = v;
        }

        if (use_limiter) {
            for (int ii = 0; ii < values_q0.size(); ++ii) {
                // if (values[ii] < 0) { values[ii] = 0; } //min_f; }
                if (values_q0[ii] < 0.0) values_q0[ii] = 0.0; //min_w0;
            }
        }
    }
}


// Robust remesh of the entire current mesh (base grid + AMR-refined points) from
// the deformed source copy. Mirrors interpolate_to_initial_xys: shift into the
// principal periodic strip, sort by (y,x), bootstrap the first point with the
// recursive search, then resolve every subsequent point by the neighbor-walk
// search seeded from the previous (spatially adjacent) point. The neighbor-walk
// re-centers the periodic image relative to the local leaf in BOTH x and y, so it
// resolves the doubly-periodic, sheared corners that the per-point recursive
// descent misroutes -- which is what produced the corner speckle on q+/q-.
void AMRStructure::interpolate_q_scattered(std::vector<double>& q0s) {
    int n = xs.size();

    std::vector<double> shifted_xs(n);
    shift_xs(shifted_xs, xs, ys);
    std::vector<double> shifted_ys(ys);

    if (bcs == periodic_bcs) {
        double x_bl = old_xs[0], y_bl = old_ys[0];
        double x_tl = old_xs[2], y_tl = old_ys[2];
        double x_br = old_xs[6], y_br = old_ys[6];
        double x_tr = old_xs[8], y_tr = old_ys[8];
        for (int ii = 0; ii < n; ++ii) {
            double x = shifted_xs[ii];
            double yt = shifted_ys[ii];
            bool ineq_bottom = (x_br - x_bl) * (yt - y_bl) >= (y_br - y_bl) * (x - x_bl);
            int counter = 0;
            while (!ineq_bottom) {
                yt += Ly;
                ineq_bottom = (x_br - x_bl) * (yt - y_bl) >= (y_br - y_bl) * (x - x_bl);
                if (++counter > 10) { throw std::runtime_error("too many y shifts at bottom (scattered)!"); }
            }
            bool ineq_top = (x_tr - x_tl) * (yt - y_tl) <= (y_tr - y_tl) * (x - x_tl);
            counter = 0;
            while (!ineq_top) {
                yt -= Ly;
                ineq_top = (x_tr - x_tl) * (yt - y_tl) <= (y_tr - y_tl) * (x - x_tl);
                if (++counter > 10) { throw std::runtime_error("too many y shifts at top (scattered)!"); }
            }
            shifted_ys[ii] = yt;
        }
    }

    std::vector<int> sort_indices(n);
    for (int ii = 0; ii < n; ++ii) { sort_indices[ii] = ii; }
    double sort_threshold = initial_dy / 10.0;
    std::sort(sort_indices.begin(), sort_indices.end(),
        [&](int a, int b) {
            if (fabs(shifted_ys[a] - shifted_ys[b]) >= sort_threshold) { return shifted_ys[a] < shifted_ys[b]; }
            return shifted_xs[a] < shifted_xs[b];
        });

    std::vector<double> sortxs(n), sortys(n), sortq0s(n);
    for (int ii = 0; ii < n; ++ii) { sortxs[ii] = shifted_xs[sort_indices[ii]]; sortys[ii] = shifted_ys[sort_indices[ii]]; }

    std::vector<int> leaf_of(n);
    bool beyond_boundary = false;
    int leaf = find_leaf_containing_xy_recursively(sortxs[0], sortys[0], beyond_boundary, 0);
    leaf_of[0] = beyond_boundary ? 0 : leaf;

    for (int ii = 1; ii < n; ++ii) {
        beyond_boundary = false;
        int seed = leaf_of[ii-1];
        if (seed <= 0) {
            // previous point fell back to root: re-bootstrap with the recursive search
            int lf = find_leaf_containing_xy_recursively(sortxs[ii], sortys[ii], beyond_boundary, 0);
            leaf_of[ii] = beyond_boundary ? 0 : lf;
        } else {
            std::set<int> history;
            history.emplace(seed);
            int lf = find_leaf_containing_point_from_neighbor(sortxs[ii], sortys[ii], beyond_boundary, seed, history);
            leaf_of[ii] = beyond_boundary ? 0 : lf;
        }
    }

    std::vector<std::vector<int> > point_in_leaf_panels_by_inds(old_panels.size());
    for (int ii = 0; ii < n; ++ii) { point_in_leaf_panels_by_inds[leaf_of[ii]].push_back(ii); }

    for (int panel_ind = 0; panel_ind < (int)old_panels.size(); ++panel_ind) {
        if (point_in_leaf_panels_by_inds[panel_ind].size() > 0) {
            interpolate_from_panel_to_points(sortq0s, sortxs, sortys,
                point_in_leaf_panels_by_inds[panel_ind], panel_ind, use_limiter, limit_val);
        }
    }
    for (int ii = 0; ii < n; ++ii) { q0s[sort_indices[ii]] = sortq0s[ii]; }
}

double AMRStructure::interpolate_from_mesh(double x, double y, bool verbose) {

    // shift x into the principal periodic strip
    std::vector<double> xs(1,x);
    std::vector<double> shifted_xs(1,x);
    std::vector<double> ys(1,y);
    shift_xs(shifted_xs, xs, ys);
    double shifted_x = shifted_xs[0];
    double shifted_y = y;

    // also pull y back into the principal periodic strip when y is periodic,
    // using the corners of the deformed copy currently in old_xs/old_ys --
    // mirrors the per-point y-shift in interpolate_to_initial_xys so that the
    // leaf search starts from the correct image before descending.
    if (bcs == periodic_bcs) {
        double x_bl = this->old_xs[0], y_bl = this->old_ys[0];
        double x_tl = this->old_xs[2], y_tl = this->old_ys[2];
        double x_br = this->old_xs[6], y_br = this->old_ys[6];
        double x_tr = this->old_xs[8], y_tr = this->old_ys[8];

        bool ineq_bottom = (x_br - x_bl) * (shifted_y - y_bl) >= (y_br - y_bl) * (shifted_x - x_bl);
        int counter = 0;
        while (!ineq_bottom) {
            shifted_y += Ly;
            ineq_bottom = (x_br - x_bl) * (shifted_y - y_bl) >= (y_br - y_bl) * (shifted_x - x_bl);
            if (++counter > 10) {
                throw std::runtime_error("too many y shifts at bottom (interpolate_from_mesh)!");
            }
        }

        bool ineq_top = (x_tr - x_tl) * (shifted_y - y_tl) <= (y_tr - y_tl) * (shifted_x - x_tl);
        counter = 0;
        while (!ineq_top) {
            shifted_y -= Ly;
            ineq_top = (x_tr - x_tl) * (shifted_y - y_tl) <= (y_tr - y_tl) * (shifted_x - x_tl);
            if (++counter > 10) {
                throw std::runtime_error("too many y shifts at top (interpolate_from_mesh)!");
            }
        }
    }

    bool beyond_boundary = false;
    // shifted_x / shifted_y are passed by reference and may be wrapped further
    // by the recursive descent across periodic boundaries; the same (wrapped)
    // values are then handed to interpolate_from_panel so the point and the
    // leaf's vertices live in the same frame.
    int leaf_containing = find_leaf_containing_xy_recursively(shifted_x, shifted_y, beyond_boundary, 0);
    if (beyond_boundary) {
        leaf_containing = 0;
    }

    double val = interpolate_from_panel(shifted_x, shifted_y, leaf_containing, verbose);
    if (verbose) {
        cout << "(" << shifted_x << ", " << shifted_y << ") is in panel " << leaf_containing << ", f_interpolated(x,y) = " << val << endl;
    }
    return val;
}


double AMRStructure::interpolate_from_panel(double x, double y, int panel_ind, bool verbose) {
    if (panel_ind == 0) { return w0_beyond_boundary; }
    else {
        Panel* panel = &(old_panels[panel_ind]);
        const int* point_inds = panel->point_inds;
        double panel_xs[9], panel_ys[9], panel_q0s[9];

        for (int ii = 0; ii < 9; ++ii) {
            int pind = point_inds[ii];
            panel_xs[ii] = old_xs[pind];
            panel_ys[ii] = old_ys[pind];
            panel_q0s[ii] = old_q0s[pind];
        }
//         if (do_unshear) {

//             double gamma = sqrt(1 + p*p);
//             double v = p * q / qm / gamma; 
//             double panel_gammas[9], panel_vs[9];
//             for (int ii = 0; ii < 9; ++ii) {
//                 panel_gammas[ii] = sqrt(1 + panel_ps[ii]*panel_ps[ii]);
//                 panel_vs[ii] = panel_ps[ii] * q / qm / panel_gammas[ii];
//                 panel_xs[ii] -= dt * (panel_vs[ii] - panel_vs[4]); // assumes remesh frequency = 1
// #ifdef DEBUG

// cout << "unshear panel x " << panel_xs[ii] << ", dt " << dt << ", p " << p << ", pmid " << panel_ps[4] << endl;
// #endif
//             }
//             x -= dt * (v - panel_vs[4]);
// #ifdef DEBUG

// cout << "unshear x " << x << ", dt " << dt << ", p " << p << ", pmid " << panel_ps[4] << endl;
// #endif

//         }


        if (verbose) {
            std::cout << "Interpolating from " << std::endl;
            for (int ii = 0; ii < 9; ++ii) {
                std::cout << "(x,y,f)=(" << panel_xs[ii] << ", " << panel_ys[ii] << ", " << panel_q0s[ii] << ")" << std::endl;
            }
        
            std::cout << "onto (" << x << ", " << y << ")\n";
        }

        double dx, dy, panel_dx[9], panel_dy[9];
        dx = x - panel_xs[4];
        dy = y - panel_ys[4];
        for (int ii = 0; ii < 9; ii ++) {
            panel_dx[ii] = panel_xs[ii] - panel_xs[4];
            panel_dy[ii] = panel_ys[ii] - panel_ys[4];
        }

        // if (verbose) {
        //     std::cout << "test point distance from midpoint: dx=" << dx <<", dp=" << dp << std::endl;
        //     std::cout << "panel vertex distances from midpoint:" << std::endl;
        //     for (int ii = 0; ii < 9; ii++ ) {
        //         std::cout << ii << ": " << panel_dx[ii] << ", " << panel_dp[ii] << std::endl;
        //     }
        // }

        Eigen::Matrix<double,9,9> A;
        for (int ii = 0; ii < 9; ++ii) {
            A(ii,0) = 1; A(ii,1) = panel_dx[ii];
            A(ii,2) = panel_dx[ii] * panel_dy[ii];
            A(ii,3) = panel_dy[ii];
            A(ii,4) = panel_dx[ii] * panel_dx[ii];
            A(ii,5) =  panel_dx[ii] * panel_dx[ii] * panel_dy[ii];
            A(ii,6) = panel_dx[ii] * panel_dx[ii] * panel_dy[ii] * panel_dy[ii];
            A(ii,7) = panel_dx[ii] * panel_dy[ii] * panel_dy[ii];
            A(ii,8) = panel_dy[ii] * panel_dy[ii];
        }
        Eigen::Map<Eigen::Matrix<double,9,1>> b(panel_q0s);
        Eigen::Matrix<double,9,1> c = A.lu().solve(b);


        if (verbose) {
            std::cout << "Here is the matrix A:\n" << A << std::endl;
            std::cout << "Here is the f vector b:\n" << b << std::endl;
            std::cout << "Here is the coefficient vector c:\n" << c << std::endl;
        }

        double val = c(0) + c(1)*dx + c(2) * dx*dy + c(3) * dy +
                c(4) * dx*dx + c(5) * dx*dx*dy + c(6) * dx*dx*dy*dy +
                c(7) * dx*dy*dy + c(8) * dy*dy;

        // Extrapolation backstop (mirrors interpolate_from_panel_to_points):
        // bound an extrapolated target to the panel's 9-node range so a
        // misrouted / far-outside point cannot diverge into an outlier.
        const double extrap_margin = 0.5;
        double hx = 0.0, hy = 0.0, qlo = panel_q0s[0], qhi = panel_q0s[0];
        for (int ii = 0; ii < 9; ++ii) {
            if (fabs(panel_dx[ii]) > hx) hx = fabs(panel_dx[ii]);
            if (fabs(panel_dy[ii]) > hy) hy = fabs(panel_dy[ii]);
            if (panel_q0s[ii] < qlo) qlo = panel_q0s[ii];
            if (panel_q0s[ii] > qhi) qhi = panel_q0s[ii];
        }
        if (fabs(dx) > hx || fabs(dy) > hy) {   // extrapolating
            double pad = extrap_margin * (qhi - qlo);
            if (val < qlo - pad) val = qlo - pad;
            if (val > qhi + pad) val = qhi + pad;
        }

        if (verbose) { cout << "Result = " << val << endl; }
        return val;
    }
}