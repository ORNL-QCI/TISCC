#include <TISCC/gridmanager.hpp>

namespace TISCC 
{
    // Constructs a particular repeating unit at a grid_[i][j]
    void GridManager::repeating_unit(unsigned int i, unsigned int j) {
        grid_[i][j] = new SiteType[7];
        grid_[i][j][0] = SiteType::QSite_Memory;
        grid_[i][j][1] = SiteType::QSite_Memory_and_Ops;
        grid_[i][j][2] = SiteType::QSite_Memory;
        grid_[i][j][3] = SiteType::Junction;
        grid_[i][j][4] = SiteType::QSite_Memory;
        grid_[i][j][5] = SiteType::QSite_Memory_and_Ops;
        grid_[i][j][6] = SiteType::QSite_Memory;
    }
    
    // Constructor for GridManager
    GridManager::GridManager(unsigned int nrows, unsigned int ncols) {
        nrows_ = nrows;
        ncols_ = ncols;
        grid_ = new SiteType**[nrows]; 
        for (unsigned int i=0; i<nrows; i++) {
            grid_[i] = new SiteType*[ncols];
            // The third level of the array contains the repeating unit of the grid.
            // The private member function repeating_unit(i, j) constructs a particular (the same) repeating unit at grid point.
            // One can easily define multiple/different such functions to use in the below loop.
            for (unsigned int j=0; j<ncols; j++) {
                repeating_unit(i, j);
            }
        }
    }

    // TODO: Provide public (read-only) access to the grid
    // GridManager::SiteType const*** GridManager::grid() {
    //     return grid_;
    // }
    GridManager::SiteType*** GridManager::grid() {
        return grid_;
    }

    // Print out grid
    void GridManager::print_grid() {
        for (unsigned int i = 0; i<nrows_; i++) {
            for (unsigned int j = 0; j<ncols_; j++) {
                std::cout << i << " " << j << std::endl;
                for (unsigned int k=0; k<7; k++) {
                    std::cout << grid_[i][j][k];
                }
                std::cout << std::endl;
            }
        }
    }
}