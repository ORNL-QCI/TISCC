#include <TISCC/gridmanager.hpp>
#include <TISCC/plaquette.hpp>

#include <iomanip>
#include <limits>
#include <cassert>
#include <map>
#include <algorithm>
#include <iostream>

namespace TISCC 
{
    // Constructor for GridManager
    GridManager::GridManager(unsigned int nrows, unsigned int ncols) {
        nrows_ = nrows;
        ncols_ = ncols;
        grid_ = new SiteType[nrows*ncols*7];
        for (unsigned int i=0; i<nrows; i++) {
            for (unsigned int j=0; j<ncols; j++) {
                grid_[(i*ncols+j)*7+0] = SiteType::QSite_Memory;
                grid_[(i*ncols+j)*7+1] = SiteType::QSite_Memory_and_Ops;
                grid_[(i*ncols+j)*7+2] = SiteType::QSite_Memory;
                grid_[(i*ncols+j)*7+3] = SiteType::Junction;
                grid_[(i*ncols+j)*7+4] = SiteType::QSite_Memory;
                grid_[(i*ncols+j)*7+5] = SiteType::QSite_Memory_and_Ops;
                grid_[(i*ncols+j)*7+6] = SiteType::QSite_Memory;
            }
        }
    }

    // Providing read-only access to array elements from coordinates
    const GridManager::SiteType& GridManager::from_coords(unsigned int row, unsigned int col, unsigned int idx) const {
        if ((row<nrows_) && (col<ncols_) && (idx<7)) {
            return grid_[(row*ncols_+col)*7+idx];
        }
        else {
            std::cerr << "GridManager::from_coords: Invalid coordinates given." << std::endl;
            abort();
        }
    } 
    
    // Providing array index from coords
    unsigned int GridManager::index_from_coords(unsigned int row, unsigned int col, unsigned int idx) const {
        if ((row<nrows_) && (col<ncols_) && (idx<7)) {
            return (row*ncols_+col)*7+idx;
        }
        else {
            std::cerr << "GridManager::index_from_coords: Invalid coordinates given." << std::endl;
            abort();
        }
    }  

    // Get a qsite shifted on the grid (index_from_coords handles validity checking)
    unsigned int GridManager::shift_qsite(unsigned int i, int nrows, int ncols) const {
        unsigned int uint_max = std::numeric_limits<unsigned int>::max(); 
        if (i == uint_max) {
            return i;
        }
        else {
            return index_from_coords(get_row(i) + nrows, get_col(i) + ncols, get_idx(i));
        }
    }

    // Provide a path from one 'O' site to the closest site next to a surrounding 'O' site (includes junctions)
    std::vector<unsigned int> GridManager::get_path(unsigned int site1, unsigned int site2) const {
        
        // Constrain ourselves to the case that these are both 'O' sites
        assert((grid_[site1] == 'O') && (grid_[site2] == 'O'));

        std::vector<unsigned int> seq; 
        seq.push_back(site1);

        unsigned int idx1 = get_idx(site1);
        unsigned int col1 = get_col(site1); 
        unsigned int row1 = get_row(site1);

        unsigned int idx2 = get_idx(site2);
        unsigned int col2 = get_col(site2); 
        unsigned int row2 = get_row(site2);

        if ((idx1 == 1) && (idx2 == 5)) {
            // NE 
            if ((col2 == col1) && (row1 == row2)) {
                seq.push_back(index_from_coords(row1, col1, idx1+1));
                seq.push_back(index_from_coords(row1, col1, idx1+2));
                seq.push_back(index_from_coords(row1, col1, idx1+3));
            }
            // NW
            else if ((col2 == col1 - 1) && (row1 == row2)) {
                seq.push_back(index_from_coords(row1, col1, idx1+1));
                seq.push_back(index_from_coords(row1, col1, idx1+2));
                seq.push_back(index_from_coords(row1, col1-1, 6));
            }
            // SE
            else if ((col2 == col1) && (row2 == row1 + 1)) {
                seq.push_back(index_from_coords(row1, col1, idx1-1));
                seq.push_back(index_from_coords(row1+1, col1, 3));
                seq.push_back(index_from_coords(row1+1, col1, 4));
            }
            // SW
            else if ((col2 == col1 - 1) && (row2 == row1 + 1)) {
                seq.push_back(index_from_coords(row1, col1, idx1-1));
                seq.push_back(index_from_coords(row1+1, col1, 3));
                seq.push_back(index_from_coords(row1+1, col1-1, 6));

            }
            else {
                std::cerr << "GridManager::get_path: Only paths to neighbor Op sites starting from the [1] position are implemented." << std::endl;
                abort();
            }
        }

        return seq;
    }

    // Return all sites adjacent to a given site
    std::set<unsigned int> GridManager::get_adjacent(unsigned int site) const {
        std::set<unsigned int> adjacent_sites;

        unsigned int idx = get_idx(site);
        unsigned int col = get_col(site);
        unsigned int row = get_row(site);

        if ((idx != 0) && (idx != 6)) {
            adjacent_sites.insert(index_from_coords(row, col, idx-1));
            adjacent_sites.insert(index_from_coords(row, col, idx+1));
        }
        else if (idx==0) {
            adjacent_sites.insert(index_from_coords(row, col, idx+1));
            if (row < nrows_ - 1) {
                adjacent_sites.insert(index_from_coords(row+1, col, 3));
            }
        }
        else if (idx==6) {
            adjacent_sites.insert(index_from_coords(row, col, idx-1));
            if (col < ncols_ - 1) {
                adjacent_sites.insert(index_from_coords(row, col+1, 3));
            }
        }
        if (idx == 3) {
            if (row > 0) {
                adjacent_sites.insert(index_from_coords(row-1, col, 0));
            }
            if (col > 0) {
                adjacent_sites.insert(index_from_coords(row, col-1, 6));   
            }
        }

        return adjacent_sites;
    }

    // Flip occupation state of two sites i.e. ``move a qubit'' (correctness relies on there only ever being one qubit per site)
    unsigned int GridManager::move_qubit(unsigned int site1, unsigned int site2) {

        // Make sure it does not move to a junction
        if (grid_[site2] == 'J') {
            std::cerr << "GridManager::move_qubit: A qubit cannot be moved to a junction." << std::endl;
            abort();
        }

        // Check whether the target site is either (a) adjacent or (b) adjacent to an adjacent junction (and if so, track that junction)
        unsigned int uint_max = std::numeric_limits<unsigned int>::max(); 
        unsigned int junction = uint_max;
        std::set<unsigned int> adjacent = get_adjacent(site1);
        if (adjacent.find(site2) == adjacent.end()) {
            bool found_J = false;
            for (unsigned int site : adjacent) {
                if (grid_[site] == 'J') {
                    found_J = true;
                    junction = site;
                    std::set<unsigned int> adj_J = get_adjacent(site);
                    if (adj_J.find(site2) == adj_J.end()) {
                        std::cerr << "GridManager::move_qubit: Attempted to move to a site both not adjacent to starting site and not adjacent to the starting site's adjacent junction." << std::endl;
                        abort();
                    }
                    else {continue;}
                }
            }
            if (!found_J) {
                std::cerr << "GridManager::move_qubit: Attempted move to non-adjacent site." << std::endl;
                abort();
            }
        }

        // Check validity of operation based on current state of occupied_sites. If valid, erase site 1 and add site2.
        if (is_occupied(site1) && !is_occupied(site2)) {
            occupied_sites.erase(site1); 
            occupied_sites.insert(site2);
        }
        else {
            std::cerr << "GridManager::move_qubit: Operation inconsistent with state of occupied_sites set." << std::endl;
            abort();
        }

        // Return qsite corresponding to the junction passed through (if applicable)
        return junction;

    }

    // Provide a plaquette object ``pinned" at a particular grid point 
    Plaquette GridManager::get_plaquette(unsigned int row, unsigned int col, char shape, char type) {
        /* Notes: 
            - shape defines either a four-qubit plaquette ('f') or a two-qubit plaquette with directionality {'n','s','e','w'}.
            - type defines whether it is an 'X' or a 'Z' plaquette. */

        // Validity check on inputs
        if (!((row<nrows_) && (col<ncols_)) && !((type=='X') || (type=='Z')) &&
            !((shape=='f') || (shape=='n') || (shape=='s') || (shape=='e') || (shape=='w'))) {
            std::cerr << "GridManager::get_plaquette: Invalid input given." << std::endl;
            abort();
        }

        // Validity checks on consistency of row, col, and shape
        if ((((col==0) || (row==nrows_-1)) && ((shape=='f') || (shape=='n') || (shape=='e'))) ||
            ((row==nrows_-1) && (shape=='w')) ||
            ((col==0) && (shape=='s')))
        {
            std::cerr << "GridManager::get_plaquette: Cannot place plaquette." << std::endl;
            abort();
        }

        // Assign index of measure qubit 
        unsigned int m = index_from_coords(row, col, 1);

        // Use the maximum possible unsigned int to designate non-existent qubits
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();   

        // Construct Plaquettes and record the relevant sites as occupied
        if (shape== 'f') {
            unsigned int a = index_from_coords(row, col-1, 5);
            unsigned int b = index_from_coords(row, col, 5);
            unsigned int c = index_from_coords(row+1, col-1, 5);
            unsigned int d = index_from_coords(row+1, col, 5);
            std::set<unsigned int> sites{a, b, c, d, m};
            occupied_sites.insert(sites.begin(), sites.end());
            return Plaquette(a,b,c,d,m,row,col,shape,type,this);
        }
        else if (shape== 'n') {
            unsigned int c = index_from_coords(row+1, col-1, 5);
            unsigned int d = index_from_coords(row+1, col, 5);
            std::set<unsigned int> sites{c, d, m};
            occupied_sites.insert(sites.begin(), sites.end());
            return Plaquette(uint_max,uint_max,c,d,m,row,col,shape,type,this);
        }
        else if (shape== 's') {
            unsigned int a = index_from_coords(row, col-1, 5);
            unsigned int b = index_from_coords(row, col, 5);
            std::set<unsigned int> sites{a, b, m};
            occupied_sites.insert(sites.begin(), sites.end());
            return Plaquette(a,b,uint_max,uint_max,m,row,col,shape,type,this);
        }
        else if (shape== 'e') {
            unsigned int a = index_from_coords(row, col-1, 5);
            unsigned int c = index_from_coords(row+1, col-1, 5);
            std::set<unsigned int> sites{a, c, m};
            occupied_sites.insert(sites.begin(), sites.end());
            return Plaquette(a,uint_max,c,uint_max,m,row,col,shape,type,this);
        }
        else if (shape== 'w') {
            unsigned int b = index_from_coords(row, col, 5);
            unsigned int d = index_from_coords(row+1, col, 5);
            std::set<unsigned int> sites{b, d, m};
            occupied_sites.insert(sites.begin(), sites.end());
            return Plaquette(uint_max,b,uint_max,d,m,row,col,shape,type,this);
        }
        else {std::cerr << "GridManager::get_plaquette: Invalid input given." << std::endl; abort();}
    }

    // Routine that enforces the validity of hardware instructions
    void GridManager::enforce_hw_master_validity(std::vector<HW_Instruction>& hw_master) {
        
        // Instantiate HardwareModel
        HardwareModel TI_model;

        // Initialize time counter and define a time_offset that will be used to resolve junction conflicts
        double time = 0; 
        double time_offset = TI_model.get_ops().at("Move") + TI_model.get_ops().at("Junction");

        // Initialize junction tracker 
        unsigned int uint_max = std::numeric_limits<unsigned int>::max(); 
        unsigned int junction = uint_max;
        std::set<unsigned int> junctions;

        // We shift time according to how many junction conflicts have been seen, with maximum one shift per time slice
        unsigned int num_shifts = 0;
        bool shift_yet = false;

        // Loop over instructions
        for (HW_Instruction& instruction : hw_master) {

            // Update time (and update/reset other variables) when we reach a new time-slice
            if (instruction.get_time() != time) {
                time = instruction.get_time();
                if (shift_yet == true) {
                    num_shifts += 1;
                    shift_yet = false;
                }
                junctions.clear();
            }

            // Perform the move_qubit operation (which contains validity checks) and obtain junction qsite if applicable
            if (instruction.get_name() == "Move") {
                junction = move_qubit(instruction.get_site1(), instruction.get_site2());
            }
            else {
                junction = uint_max;
            }

            // Check for junction conflicts in this time step and note a time shift if one is found
            if (junction != uint_max) {
                if (junctions.find(junction) == junctions.end()) {
                    junctions.insert(junction);
                }
                else {
                    if (!shift_yet) {
                        shift_yet = true;
                    }
                    // If a collision is found, shift the time for this instruction
                    instruction = HW_Instruction(instruction, (num_shifts+1)*time_offset);
                    // std::cerr << "GridManager::check_hw_master_validity: Junction " << junction << " is being used more than once at t = " << time + num_shifts*time_offset << " us." << std::endl;
                    continue;
                }
            }

            // Update time by copying with time_offset
            if (num_shifts != 0) {
                instruction = HW_Instruction(instruction, num_shifts*time_offset); 
            }
        }

        // Sort updated hw_master
        std::stable_sort(hw_master.begin(), hw_master.end());
    }

    // Routine that counts resources
    // TODO: Re-think whether this belongs as a GridManager member function
    void GridManager::resource_counter(const std::vector<HW_Instruction>& hw_master) const {

        // IO Settings
        std::cout.precision(3);

        // Instantiate HardwareModel
        HardwareModel TI_model;

        // Obtain linear dimensions of the grid (four trapping zone widths per row or column)
        double cell_width = TI_model.get_cell_width() / 1000000; // convert to meters
        double z_dim = nrows_ * cell_width; 
        double x_dim = ncols_ * cell_width;

        // Total number of sites
        unsigned int n_sites = nrows_*ncols_*7;

        // Computation time
        double time = hw_master.back().get_time();
        double time_tmp = 0;

        // Map to compute zone-seconds per op
        std::unordered_map<std::string, double> op_vol;
        for (const auto & [ key, value ] : TI_model.get_ops()) {
            op_vol[key] = 0;
        }

        // Loop over instructions
        for (const HW_Instruction& instruction : hw_master) {

            /* If the instruction is a move, figure out if it is a junction move
            (we have assumed all input instructions are valid, so any non-adjacent move is a junction move) */
            std::string name = instruction.get_name();
            if (name == "Move") {
                std::set<unsigned int> adjacent = get_adjacent(instruction.get_site1());
                if (adjacent.find(instruction.get_site2()) == adjacent.end()) {name = "Junction";}
            }

            // Add contribution to volume per op
            op_vol[name]++;

            // Get longest operation in last slice
            if (instruction.get_time() == time) {
                if (TI_model.get_ops().at(instruction.get_name()) > time_tmp) {
                    time_tmp = TI_model.get_ops().at(instruction.get_name());
                }
            }

        }

        // Print out gate counts and compute zone-seconds for each operation
        std::cout << std::fixed;
        double total = 0;
        for (auto & [ key, value ] : op_vol) {
            if (value != 0) {
                std::cout << key << ": " << int(value) << std::fixed << std::endl; 
                value = value * TI_model.get_ops().at(key) / 1000000;
                if (key == "Move") {
                    value = value * 2;
                }
                else if (key == "Junction") {
                    value = value * 3;
                }
                total += value;
            }
        }

        // Add longest operation in last slice to computation time
        time += time_tmp;

        // Rescale 
        time = time/1000000; // s

        // Print out resource counts
        std::cout << std::scientific;
        std::cout << "Grid area: " << x_dim*z_dim << " m^2." << std::endl;
        std::cout << "Computation time: " << time << " s." << std::endl;
        std::cout << "Space-time volume: " << time * x_dim * z_dim << " s*m^2." << std::endl;
        std::cout << "Trapping zones: " << n_sites << std::endl;
        std::cout << "Trapping zone-seconds: " << n_sites*time << " zone*s" << std::endl;
        std::cout << "Trapping zone-seconds (active): " << total << " zone*s" << std::endl;
        std::cout << "Active zone-seconds (%): " << std::fixed << 100 * total / (n_sites*time) << std::endl;
        for (auto & [ key, value ] : op_vol) {
            if (value != 0) {
                std::cout << key << " zone-seconds (%): " << std::fixed << 100 * value / (n_sites*time) << std::endl;
            }
        }

    }

    // Print out grid
    void GridManager::print_grid() const {

        // I/O settings
        int W = 5;
        std::cout << std::setprecision(1);
        std::cout << std::setiosflags(std::ios::fixed);

        std::cout << std::setw(W) << "Row";
        std::cout << std::setw(W) << "Col";
        std::cout << std::setw(33) << "Qsites (repeating unit: MOMJMOM)" << std::endl;

        for (unsigned int i=0; i<nrows_; i++) {
            for (unsigned int j=0; j<ncols_; j++) {

                // Print out row and column
                std::cout << std::setw(W) << i;
                std::cout << std::setw(W) << j;
                std::cout << std::setw(W);

                // Print out qsites of repeating unit
                for (unsigned int k=0; k<7; k++) {
                    std::cout << std::setw(W-1) << (i*ncols_+j)*7+k << " ";
                }
                std::cout << std::endl;
            }
        }
    }

    // Print out occupied sites
    void GridManager::print_occ_sites() const {
        for (unsigned int i : occupied_sites) {
            std::cout << i << std::endl;
        }
    }

    // Visualize the grid using hard-coded repeating unit
    void GridManager::visualize_grid() const {
        std::vector<std::string> repeating_unit = {
            "J M O M ",
            "M       ",
            "O       ",
            "M       "
        };

        std::cout << std::endl << "  ";
        for (unsigned int j = 0; j < ncols_; j++) {
            std::cout << j << "       ";
        }
        std::cout << std::endl;

        for (unsigned int i=0; i<nrows_; i++) {
            for (unsigned int k=0; k<repeating_unit.size(); k++) {
                if (k==0) std::cout << i << " ";
                else std::cout << "  ";
                for (unsigned int j=0; j<ncols_; j++) {
                    std::cout << repeating_unit[k];
                }
                std::cout << std::endl;
            }
        }
    }

    // Visualize an operator on the grid
    void GridManager::visualize_operator(const std::vector<std::pair<unsigned int,char>>& qsites) const {
        std::vector<std::string> repeating_unit = {
            "J M O M ",
            "M       ",
            "O       ",
            "M       "
        };

        std::vector<std::set<unsigned int>> repeating_unit_idxs = {
            {3, 4, 5, 6},
            {2},
            {1},
            {0}
        };

        // Make sure qsites is sorted
        std::vector<std::pair<unsigned int,char>> sorted_qsites = qsites;
        std::sort(sorted_qsites.begin(), sorted_qsites.end());

        // Collect row, col, and index for each qsite in the operator
        std::vector<unsigned int> rows(sorted_qsites.size());
        std::vector<unsigned int> cols(sorted_qsites.size());
        std::vector<unsigned int> indices(sorted_qsites.size());
        for (unsigned int i=0; i<sorted_qsites.size(); i++) {
            rows[i] = get_row(sorted_qsites[i].first);
            cols[i] = get_col(sorted_qsites[i].first);
            indices[i] = get_idx(sorted_qsites[i].first);
        }

        /* sorted_qsites will be hit in order when looping sequentially over rows and cols in the grid */

        // Loop over rows, cols, rows of repeating unit
        std::cout << std::endl << "  ";
        for (unsigned int j = 0; j < ncols_; j++) {
            std::cout << j << "       ";
        }
        std::cout << std::endl;
        unsigned int qsite_tracker = 0;
        for (unsigned int i = 0; i < nrows_; i++) {
            for (unsigned int k = 0; k < repeating_unit.size(); k++) {
                if (k==0) std::cout << i << " ";
                else std::cout << "  ";
                for (unsigned int j = 0; j < ncols_; j++) {

                    // Check if we have hit the next qsite in the list
                    if (qsite_tracker < sorted_qsites.size()) {
                        if ((i == rows[qsite_tracker]) && (j == cols[qsite_tracker]) && (repeating_unit_idxs[k].count(indices[qsite_tracker]) == 1)) {
                            if (grid_[sorted_qsites[qsite_tracker].first] == SiteType::QSite_Memory_and_Ops) {
                                std::cout << repeating_unit[k].replace(repeating_unit[k].find(grid_[sorted_qsites[qsite_tracker].first]), 1, 1, sorted_qsites[qsite_tracker].second);
                                repeating_unit[k].replace(repeating_unit[k].find(sorted_qsites[qsite_tracker].second), 1, 1, grid_[sorted_qsites[qsite_tracker].first]);
                                qsite_tracker++;
                            }
                            else {
                                std::cerr << "GridManager::visualize_operator: Only SiteType::QSite_Memory_and_Ops are currently used as home bases for qubits in our hardware model." << std::endl;
                                abort();
                            }
                        }
                        else {
                            std::cout << repeating_unit[k];
                        }
                    }
                    else {
                        std::cout << repeating_unit[k];
                    }
                }

                std::cout << std::endl;
            }
        }
    }
}