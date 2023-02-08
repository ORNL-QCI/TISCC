#include <TISCC/gridmanager.hpp>
#include <TISCC/plaquette.hpp>

#include <limits>
#include <cassert>

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

    // Provide a path from one 'O' site to the closest site next to a surrounding 'O' site
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

    // Flip occupation state of two sites i.e. ``move a qubit'' (correctness relies on there only ever being one qubit per site)
    void GridManager::move_qubit(unsigned int site1, unsigned int site2) {

        // Check validity of operation based on current state of occupied_sites
        if (is_occupied(site1) && !is_occupied(site2)) {
            //std::cerr << site1 << " " << is_occupied(site1) << " " << site2 << " " << is_occupied(site2) << std::endl;
            occupied_sites.erase(site1); 
            occupied_sites.insert(site2);
        }
        else {
            //std::cerr << site1 << " " << is_occupied(site1) << " " << site2 << " " << is_occupied(site2) << std::endl;
            std::cerr << "GridManager::move_qubit: operation inconsistent with state of occupied_sites set." << std::endl;
            abort();
        }
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
            return Plaquette(a,b,c,d,m,row,col,shape,type,*this);
        }
        else if (shape== 'n') {
            unsigned int c = index_from_coords(row+1, col-1, 5);
            unsigned int d = index_from_coords(row+1, col, 5);
            std::set<unsigned int> sites{c, d, m};
            occupied_sites.insert(sites.begin(), sites.end());
            return Plaquette(uint_max,uint_max,c,d,m,row,col,shape,type,*this);
        }
        else if (shape== 's') {
            unsigned int a = index_from_coords(row, col-1, 5);
            unsigned int b = index_from_coords(row, col, 5);
            std::set<unsigned int> sites{a, b, m};
            occupied_sites.insert(sites.begin(), sites.end());
            return Plaquette(a,b,uint_max,uint_max,m,row,col,shape,type,*this);
        }
        else if (shape== 'e') {
            unsigned int a = index_from_coords(row, col-1, 5);
            unsigned int c = index_from_coords(row+1, col-1, 5);
            std::set<unsigned int> sites{a, c, m};
            occupied_sites.insert(sites.begin(), sites.end());
            return Plaquette(a,uint_max,c,uint_max,m,row,col,shape,type,*this);
        }
        else if (shape== 'w') {
            unsigned int b = index_from_coords(row, col, 5);
            unsigned int d = index_from_coords(row+1, col, 5);
            std::set<unsigned int> sites{b, d, m};
            occupied_sites.insert(sites.begin(), sites.end());
            return Plaquette(uint_max,b,uint_max,d,m,row,col,shape,type,*this);
        }
        else {std::cerr << "GridManager::get_plaquette: Invalid input given." << std::endl; abort();}
    }

    // Print out grid
    void GridManager::print_grid() const {
        for (unsigned int i=0; i<nrows_; i++) {
            for (unsigned int j=0; j<ncols_; j++) {
                std::cout << i << " " << j << " ";
                for (unsigned int k=0; k<7; k++) {
                    std::cout << grid_[(i*ncols_+j)*7+k];
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
}