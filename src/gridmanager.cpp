#include <TISCC/gridmanager.hpp>
#include <TISCC/plaquette.hpp>

#include <limits>

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