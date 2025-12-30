#include "AMRStructure.hpp"

void AMRStructure::interpolate_to_initial_xys(
    std::vector<double>& fs, std::vector<double>& xs, std::vector<double>& ys, 
    int nx, int ny, bool verbose) 
{
    std::vector<double> shifted_xs(xs.size() );
    shift_xs(shifted_xs, xs, ys); // shift xs to the bounded domain 

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
    // find the leaf index for the first point first 
    int leaf_ind = find_leaf_containing_xy_recursively(sortxs[0], sortys[0], beyond_boundary, 0, verbose);
    std::vector<int> first_column_leaf_inds(ny);

    // utilize periodic bc in x, find leaf index for the first point in each row from the previous first point 
    for (int ii = 0; ii < ny; ++ii) {
        beyond_boundary = false;
        int point_ind = ii * nx;
        std::set<int> history;
        history.emplace(leaf_ind);

        leaf_ind = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortys[point_ind], beyond_boundary, 
                                                                leaf_ind, history, verbose);
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
            leaf_ind_c = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortys[point_ind], beyond_boundary, leaf_ind_c, history, verbose);
            if (beyond_boundary) {
                leaf_panel_of_points[point_ind] = 0;
            } else {
                leaf_panel_of_points[point_ind] = leaf_ind_c;
            }
        }
    }    

    for (int ii = 0; ii < leaf_panel_of_points.size(); ++ii) {
        point_in_leaf_panels_by_inds[leaf_panel_of_points[ii]].push_back(ii);
    }

    // eval interpolant 

    for (int panel_ind = 0; panel_ind < old_panels.size(); panel_ind++) {
        if (point_in_leaf_panels_by_inds[panel_ind].size() > 0) {
            interpolate_from_panel_to_points(sortfs,sortxs,sortys,point_in_leaf_panels_by_inds[panel_ind], panel_ind, use_limiter, limit_val);
        }
    }
    for (int ii = 0; ii < fs.size(); ii++) {
        fs[sort_indices[ii]] = sortfs[ii];
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



int AMRStructure::find_leaf_containing_point_from_neighbor(double& tx, double& ty, bool& beyond_boundary, int leaf_ind, std::set<int>& history, bool verbose) {

    if (leaf_ind == 0) {
        history.emplace(leaf_ind);
        cout << "If you see this then you are sending in a leaf ind 0 somewhere you hoped not to." << endl;
        return leaf_ind;
    } else {
        Panel* panel = &(old_panels[leaf_ind]);
        double x_bl = old_xs[panel->point_inds[0]]; double y_bl = old_ys[panel->point_inds[0]];
        double x_tl = old_xs[panel->point_inds[2]]; double y_tl = old_ys[panel->point_inds[2]];
        double x_mid = old_xs[panel->point_inds[4]];
        double x_br = old_xs[panel->point_inds[6]]; double y_br = old_ys[panel->point_inds[6]];
        double x_tr = old_xs[panel->point_inds[8]]; double y_tr = old_ys[panel->point_inds[8]];
        if (verbose) {
            cout << "(x,y)_bl = (" << x_bl << ", " << y_bl << ")" << endl;
            cout << "(x,y)_tl = (" << x_tl << ", " << y_tl << ")" << endl;
            cout << "(x,y)_br = (" << x_br << ", " << y_br << ")" << endl;
            cout << "(x,y)_tr = (" << x_tr << ", " << y_tr << ")" << endl;
        }
        // need to correct periodic distance
        if (tx - x_mid >= Lx/2 && bcs == 0) { 
            tx -= Lx; 
            if (verbose) {
                cout << "shifting across boundary, tx= " << tx << endl;
            }
        }
        if (tx - x_mid < -Lx/2 && bcs ==0) { 
            tx += Lx; 
            if (verbose) {
                cout << "shifting across boundary, tx= " << tx << endl;
            }
        }


        bool ineq_right = (x_tr - x_br) * (ty - y_br) > (y_tr - y_br) * (tx - (x_br));
        bool ineq_left = (x_tl - x_bl) * (ty - y_bl) <= (y_tl - y_bl) * (tx - x_bl);
        bool ineq_top = (x_tr - x_tl) * (ty - y_tl) < (y_tr - y_tl) * (tx - x_tl);
        bool ineq_bottom = (x_br - x_bl) * (ty - y_bl) >= (y_br - y_bl) * (tx - x_bl);
        int new_leaf_ind = leaf_ind;
        if (verbose) {
            cout << "(tx,ty) = (" << tx << ", " << ty << ")" << endl;
            cout << "testing leaf panel " << leaf_ind << " for containment" << endl;
            cout << "ineq_right " << ineq_right << endl;
            cout << "(x_tr - x_tl) * (tv - v_tl)" << (x_tr - x_tl) * (ty - y_tl) << endl;
            cout << "(v_tr - v_tl) * (tx - x_tl)" << (y_tr - y_tl) * (tx - x_tl) << endl;
            cout << "ineq_top " << ineq_top << endl;
            cout << "ineq_bottom " << ineq_bottom << endl;
            cout << "ineq_left " << ineq_left << endl;
        }
        if (! ineq_right && panel->right_nbr_ind != -2) {
            if (panel->is_right_bdry && bcs == 0) { 
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
                    if (verbose) {
                        cout << "in panel right" << endl;
                        cout << "next leaf test " << new_leaf_ind << endl;
                    }
                }
            }
        } else {
            if (!ineq_top && panel->top_nbr_ind != -2) {
                if (panel->top_nbr_ind == -1) {
                    new_leaf_ind = old_panels[panel->parent_ind].top_nbr_ind;
                    if (verbose) {
                        cout << "in parent top" << endl;
                        cout << "next leaf test " << new_leaf_ind << endl;
                    }
                } else {
                    Panel* panel_top = &old_panels[panel->top_nbr_ind];
                    if (panel_top->is_refined_xy) {
                        new_leaf_ind = panel_top->child_inds_start;
                        if (verbose) {
                            cout << "in top children" << endl;
                            cout << "next leaf test " << new_leaf_ind << endl;
                        }
                    } else {
                        new_leaf_ind = panel->top_nbr_ind;
                        if (verbose) {
                            cout << "in panel top" << endl;
                            cout << "next leaf test " << new_leaf_ind << endl;
                        }
                    }
                }
            } // end if checking top boundary 
            else {
                if (!ineq_bottom && panel->bottom_nbr_ind != -2) {
                    if (panel->bottom_nbr_ind == -1) {
                        new_leaf_ind = old_panels[panel->parent_ind].bottom_nbr_ind;
                        if (verbose) {
                            cout << "in parent bottom" << endl;
                            cout << "next leaf test " << new_leaf_ind << endl;
                        }
                    } else {
                        Panel* panel_bottom = &old_panels[panel->bottom_nbr_ind];
                        if (panel_bottom -> is_refined_xy) {
                            new_leaf_ind = panel_bottom->child_inds_start + 1;
                            if (verbose) {
                                cout << "in parent bottom children" << endl;
                                cout << "next leaf test " << new_leaf_ind << endl;
                            }
                        } else {
                            new_leaf_ind = panel->bottom_nbr_ind;
                            if (verbose) {
                                cout << "in parent bottom" << endl;
                                cout << "next leaf test " << new_leaf_ind << endl;
                            }
                        }
                    }
                } // end if checking bottom boundary 
                else {
                    if (! ineq_left && panel->left_nbr_ind != -2) {
                        if (panel->is_left_bdry && bcs == 0) { 
                            tx += Lx; 
                            if (verbose) {
                                cout << "shifting across boundary, now tx= " << tx << endl;
                            }
                        }
                        if (panel->left_nbr_ind == -1) {
#ifdef DEBUG
cout << "panel left: " << panel->left_nbr_ind << endl;
cout <<"length of panels_list " << old_panels.size() << endl;
#endif
                            new_leaf_ind = old_panels[panel->parent_ind].left_nbr_ind;
                            if (verbose) {
                                cout << "in parent left" << endl;
                                cout << "next leaf test " << new_leaf_ind << endl;
                            }
                        } else {
                            Panel* panel_left = &old_panels[panel->left_nbr_ind];
                            if (panel_left->is_refined_xy) {
                                new_leaf_ind = panel_left->child_inds_start+2;
                                if (verbose) {
                                    cout << "in left children" << endl;
                                    cout << "next leaf test " << new_leaf_ind << endl;
                                }
                            }
                            else {
                                new_leaf_ind = panel->left_nbr_ind;
                                if (verbose) {
                                    cout << "in panel left" << endl;
                                    cout << "next leaf test " << new_leaf_ind << endl;
                                }
                            }
                        }
                    } // end if checking left boundary
                } // end else for left
            } //end else for bottom/left
        } // end else for top/bottom/left
        if (verbose) {
            cout << "History" << endl;
            for (std::set<int>::iterator it = history.begin(); it != history.end(); ++it) {
                cout << *it << " ";
            }
            cout << endl;
        }
        if (history.find(new_leaf_ind) == history.end()) {
            history.emplace(new_leaf_ind);
            if (verbose) {
                cout << "running find_leaf again for (x,y)=(" << tx << ", " << ty << ")" << endl;
            }
            new_leaf_ind = find_leaf_containing_point_from_neighbor(tx,tv,beyond_boundary, new_leaf_ind, history, verbose);
        } else if (verbose)
        {
            cout << "done searching." << endl;
        }
        if (!allow_boundary_extrapolation) {
            bool boundary_extrapolating_right = !ineq_right && panel->right_nbr_ind==-2;
            bool boundary_extrapolating_left = !ineq_left && panel->left_nbr_ind==-2;
            bool boundary_extrapolating_top = !ineq_top && panel->top_nbr_ind==-2;
            bool boundary_extrapolating_bottom = !ineq_bottom && panel->bottom_nbr_ind==-2;
            if (boundary_extrapolating_left || boundary_extrapolating_right || 
                boundary_extrapolating_top || boundary_extrapolating_bottom) {
                if (verbose) {
                    cout << "Not allowing boundary extrapolation and this point is flagged as beyond the boundary" << endl;
                }
                beyond_boundary = true;
            }
        }
        if (verbose) {
            cout << "Point (" << tx << ", " << tv << ") is in panel " << new_leaf_ind << endl;
        }
        return new_leaf_ind;
    }
}


void AMRStructure::interpolate_from_panel_to_points(
    std::vector<double>& values, std::vector<double>& xs, std::vector<double>& ys,
    std::vector<int>& point_inds, int panel_ind, bool use_limiter, double limit_val) 
{
    if (panel_ind == 0) { // if we are extrapolating beyond boundaries; assume 0
        for (int ii = 0; ii < point_inds.size(); ++ii) {
            values[point_inds[ii]] = f_beyond_boundary;
        }
    }
    else {
        bool verbose = false;
        Panel* panel = &(old_panels[panel_ind]);
        const int* panel_point_inds = panel->point_inds;
        double panel_xs[9], panel_ys[9];
        // , panel_fs[9];
        double panel_w0s[9], panel_j0s[9];

        for (int ii = 0; ii < 9; ++ii) {
            int pind = panel_point_inds[ii];
            panel_xs[ii] = old_xs[pind];
            panel_ys[ii] = old_ys[pind];
            // panel_fs[ii] = old_fs[pind];
            panel_w0s[ii] = old_w0s[pind];
            panel_j0s[ii] = old_j0s[pind];

        }

        if (do_unshear) {
            for (int ii = 0; ii < 9; ++ii) {
                panel_xs[ii] -= dt * (panel_ys[ii] - panel_ys[4]); // assumes remesh frequency = 1
            }
            for (int ii = 0; ii < xs.size(); ++ii) {
                xs[ii] -= dt * (ys[ii] - panel_ys[4]);

            }
        }

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

        if (verbose) {
            std::cout << "test point distance from midpoint:" << std::endl;
            for (int ii = 0; ii < point_inds.size(); ++ii) {
                std::cout << ii <<": dx=" << dxs[ii] <<", dy=" << dys[ii] << std::endl;
            } 
            std::cout << "panel vertex distances from midpoint:" << std::endl;
            for (int ii = 0; ii < 9; ii++ ) {
                std::cout << ii << ": " << panel_dx[ii] << ", " << panel_dy[ii] << std::endl;
            }
        }

        Eigen::Matrix<double,9,9> A;
        for (int ii = 0; ii < 9; ++ii) {
            A(ii,0) = 1; A(ii,1) = panel_dx[ii];
            A(ii,2) = panel_dx[ii] * panel_dy[ii];
            A(ii,3) = panel_dy[ii];
            A(ii,4) = panel_dy[ii] * panel_dy[ii];
            A(ii,5) =  panel_dy[ii] * panel_dy[ii] * panel_dy[ii];
            A(ii,6) = panel_dy[ii] * panel_dy[ii] * panel_dy[ii] * panel_dy[ii];
            A(ii,7) = panel_dy[ii] * panel_dy[ii] * panel_dy[ii];
            A(ii,8) = panel_dy[ii] * panel_dy[ii];
        }
        Eigen::Map<Eigen::Matrix<double,9,1>> b(panel_fs);
        Eigen::Matrix<double,9,1> c = A.lu().solve(b);


        if (verbose) {
            std::cout << "Here is the matrix A:\n" << A << std::endl;
            std::cout << "Here is the f vector b:\n" << b << std::endl;
            std::cout << "Here is the coefficient vector c:\n" << c << std::endl;
        }
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
        Eigen::Matrix<double, Dynamic,1> interp_vals = Dx * c;
        if (verbose) {
            std::cout << "Here is the result:" << std::endl << interp_vals << std::endl;
        }
        for (int ii = 0; ii < point_inds.size(); ++ii) {
            values[point_inds[ii]] = interp_vals(ii);
        }

        if (use_limiter) {
            for (int ii = 0; ii < values.size(); ++ii) {
                if (values[ii] < 0) { values[ii] = 0; } //min_f; }
            }
        }
    }
}