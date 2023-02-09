#ifndef TISCC_GRIDMANAGER_HPP
#define TISCC_GRIDMANAGER_HPP

#include <TISCC/plaquette.hpp>
#include <TISCC/instruction.hpp>

#include <iostream>
#include <vector>
#include <set>

namespace TISCC {

// Need to declare this in advance since the below Class depends on it
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
    const std::set<unsigned int>& get_occ_sites() const {return occupied_sites;}

    // Provide a path from one 'O' site to the closest site next to a surrounding 'O' site
    std::vector<unsigned int> get_path(unsigned int site1, unsigned int site2) const;

    // Return all sites adjacent to a given site
    std::set<unsigned int> get_adjacent(unsigned int site) const;

    // Check if a particular site is occupied
    bool is_occupied(unsigned int site) const {return (occupied_sites.find(site) != occupied_sites.end());}

    // Flip occupation state of two sites i.e. ``move a qubit'' (relies on there only ever being one qubit per site)
    void move_qubit(unsigned int site1, unsigned int site2);

    // Provide a plaquette object ``pinned" at a particular grid point 
    Plaquette get_plaquette(unsigned int row, unsigned int col, char shape, char type);

    // Routine that enforces the validity of hardware instructions
    void check_hw_master_validity(const std::vector<HW_Instruction>& hw_master);

    // Methods for printing
    void print_grid() const;
    void print_occ_sites() const;

private:
    SiteType* grid_;
    unsigned int nrows_;
    unsigned int ncols_;
    std::set<unsigned int> occupied_sites;
};
}

#endif //TISCC_GRIDMANAGER_HPP