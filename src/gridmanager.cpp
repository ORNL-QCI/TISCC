#include <TISCC/gridmanager.hpp>

#include <optional>

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