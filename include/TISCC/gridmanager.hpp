#ifndef TISCC_GRIDMANAGER_HPP
#define TISCC_GRIDMANAGER_HPP

#include <TISCC/plaquette.hpp>

#include <iostream>
#include <vector>
#include <optional>

namespace TISCC {

class Plaquette;

// Provides access to an array representing a particular trapped-ion hardware architecture 
class GridManager {
public:
    // Define possible types of sites in the hardware layout
    enum SiteType : char {
        QSite_Memory = 'M',
        QSite_Memory_and_Ops = 'O',
        Junction = 'J',
        Nothing = 'X'
    };

    // Constructor for GridManager
    explicit GridManager(unsigned int nrows, unsigned int ncols);

    // Destructor for GridManager
    ~GridManager() {delete[] grid_;}

    // Providing read-only access to array elements 
    const SiteType& operator[] (unsigned int i) const {return grid_[i];}

    // Providing read-only access to array elements from coordinates
    const SiteType& from_coords(unsigned int row, unsigned int col, unsigned int idx) const; 

    // Providing array index from coords
    unsigned int index_from_coords(unsigned int row, unsigned int col, unsigned int idx) const;

    // Providing row, col, and idx from array index
    unsigned int get_idx(unsigned int i) const {return i%7;}
    unsigned int get_col(unsigned int i) const {return ((i-i%7)/7)%ncols_;}
    unsigned int get_row(unsigned int i) const {return (((i-i%7)/7)-get_col(i))/ncols_;}

    // Accessor functions for private member variables
    unsigned int get_nrows() const {return nrows_;}
    unsigned int get_ncols() const {return ncols_;}

    // Provide a path from one site to the site next to another
    std::vector<unsigned int> get_path(unsigned int site1, unsigned int site2) const;

    // Provide a plaquette object ``pinned" at a particular grid point 
    Plaquette get_plaquette(unsigned int row, unsigned int col, char shape, char type) const;

    // Print out grid
    void print_grid();

private:
    SiteType* grid_;
    unsigned int nrows_;
    unsigned int ncols_;
};
}

#endif //TISCC_GRIDMANAGER_HPP