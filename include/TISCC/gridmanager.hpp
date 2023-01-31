#ifndef TISCC_GRIDMANAGER_HPP
#define TISCC_GRIDMANAGER_HPP

#include <iostream>

namespace TISCC {

class GridManager {
public:
    // Define possible types of sites in the trapped ion micro layout
    enum SiteType : char {
        QSite_Memory = 'M',
        QSite_Memory_and_Ops = 'O',
        Junction = 'J',
        Nothing = 'X'
    };

    // Constructor for GridManager
    explicit GridManager(unsigned int nrows, unsigned int ncols);

    // Destructor for GridManager
    ~GridManager() {
        for (unsigned int i=0; i<nrows_; i++) {
            for (unsigned int j=0; j<ncols_; j++) {
                delete[] grid_[i][j];
            }
            delete[] grid_[i];
        }
        delete[] grid_;
    }

    // TODO: Provide public (read-only) access to the grid
    // SiteType const*** grid();
    SiteType*** grid();

    // Print out grid
    void print_grid();

private:
    // Should I use std::vector<std::vector<std::vector<SiteType>>> instead?
    // Is three dimensions too unwieldy? Does it matter?
    SiteType*** grid_;
    unsigned int nrows_;
    unsigned int ncols_;

    // TODO: Include a hash table of memory address to qsite indices for later output 
    //std::unordered_map<SiteType*, int> qsite_hash

    // Constructs a particular repeating unit at a grid_[i][j]
    void repeating_unit(unsigned int i, unsigned int j);
};
}

#endif //TISCC_GRIDMANAGER_HPP