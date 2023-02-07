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

    // Provide a path from one 'O' site to the nearest site adjacent to a neighboring 'O' site
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
                std::cout << "GridManager::get_path: Only paths to surrounding Op sites starting from the [1] position are implemented." << std::endl;
                abort();
            }
        }

        return seq;
    }

    // Provide a plaquette object ``pinned" at a particular grid point 
    Plaquette GridManager::get_plaquette(unsigned int row, unsigned int col, char shape, char type) const {
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

        // Construct and return Plaquette objects
        if (shape== 'f') {
            unsigned int a = index_from_coords(row, col-1, 5);
            unsigned int b = index_from_coords(row, col, 5);
            unsigned int c = index_from_coords(row+1, col-1, 5);
            unsigned int d = index_from_coords(row+1, col, 5);
            return Plaquette(a,b,c,d,m,row,col,shape,type,*this);
        }
        else if (shape== 'n') {
            unsigned int c = index_from_coords(row+1, col-1, 5);
            unsigned int d = index_from_coords(row+1, col, 5);
            return Plaquette(uint_max,uint_max,c,d,m,row,col,shape,type,*this);
        }
        else if (shape== 's') {
            unsigned int a = index_from_coords(row, col-1, 5);
            unsigned int b = index_from_coords(row, col, 5);
            return Plaquette(a,b,uint_max,uint_max,m,row,col,shape,type,*this);
        }
        else if (shape== 'e') {
            unsigned int a = index_from_coords(row, col-1, 5);
            unsigned int c = index_from_coords(row+1, col-1, 5);
            return Plaquette(a,uint_max,c,uint_max,m,row,col,shape,type,*this);
        }
        else if (shape== 'w') {
            unsigned int b = index_from_coords(row, col, 5);
            unsigned int d = index_from_coords(row+1, col, 5);
            return Plaquette(uint_max,b,uint_max,d,m,row,col,shape,type,*this);
        }
        else {std::cerr << "GridManager::get_plaquette: Invalid input given." << std::endl; abort();}
    }

    // Print out grid
    void GridManager::print_grid() {
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
}