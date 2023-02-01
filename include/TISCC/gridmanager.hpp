#ifndef TISCC_GRIDMANAGER_HPP
#define TISCC_GRIDMANAGER_HPP

#include <iostream>
#include <vector>

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
        delete[] grid_;
    }

    // Providing read-only access to array elements 
    const SiteType& operator[] (unsigned int i) const {return grid_[i];}

    // Providing read-only access to array elements from coordinates
    const SiteType& val_from_coords(unsigned int row, unsigned int col, unsigned int idx) const 
        {return grid_[(row*ncols_+col)*7+idx];}

    // Provide read-only access to array pointers
    SiteType const* grid(unsigned int i) const {return grid_ + i;}
    // TODO: Expand above into an iterator class (see below commented block).
    // Provide indices and/or pointers to the location north, south, east, or west of a given index and/or pointer

    /* All public member functions below depend on our specific choice of repeating unit (see constructor) */
    // Providing row, col, and idx from array index
    unsigned int get_idx(unsigned int i) {return i%7;}
    unsigned int get_col(unsigned int i) {return ((i-i%7)/7)%ncols_;}
    unsigned int get_row(unsigned int i) {return (((i-i%7)/7)-get_col(i))/ncols_;}

    // Accessor functions for private member variables
    unsigned int get_nrows() const {return nrows_;}
    unsigned int get_ncols() const {return ncols_;}

    // Print out grid
    void print_grid();

private:
    SiteType* grid_;
    unsigned int nrows_;
    unsigned int ncols_;
};

// class SiteIterator {
// public:
//     // Default constructor, copy constructor, assignment operator, & destructor
//     SiteIterator() : ptr_(NULL) {}
//     SiteIterator(SiteType* p) : ptr_(p) {}
//     SiteIterator(const SiteIterator& old) : ptr_(old.ptr_) {}
//     SiteIterator& operator= (const SiteIterator& old) {
//         ptr_ = old.ptr_; return *this;}
//     ~SiteIterator() {}

//     // Dereferencing operator
//     SiteType& operator*() {return *ptr_;}

//     // Returns an iterator pointing to the site north of the current one
//     SiteIterator& n() {

//     }

// private:
//     GridManager::SiteType* ptr_;
// }
}

#endif //TISCC_GRIDMANAGER_HPP