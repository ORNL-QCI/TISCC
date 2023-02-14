#ifndef TISCC_PLAQUETTE_HPP
#define TISCC_PLAQUETTE_HPP

#include <TISCC/gridmanager.hpp>
#include <TISCC/instruction.hpp>

#include <iostream>
#include <vector>
#include <set>

namespace TISCC {

// Need to declare this in advance since the below Class depends on it
class GridManager;
class Instruction;

// Stores the grid indices (qsites) occupied by named qubits of a surface code plaquette
class Plaquette {
public:
    // Constructor
    explicit Plaquette(unsigned int a, unsigned int b, unsigned int c, unsigned int d, unsigned int m, unsigned int row,
        unsigned int col, char shape, char type, GridManager& grid) :
        a_(a), b_(b), c_(c), d_(d), m_(m), row_(row), col_(col), shape_(shape), type_(type), grid_(grid) {}

    // Accessor functions
    unsigned int get_qsite(char qubit) const;
    unsigned int get_row() const {return row_;}
    unsigned int get_col() const {return col_;}
    char get_shape() const {return shape_;}
    char get_type() const {return type_;}
    
    // Move qubit to a specified site after checking validity and recording change on grid
    void move_to_site(char q, unsigned int site);

    // Check whether a qubit is currently at its home location on the grid
    bool is_home(char qubit) const;

    bool is_instr_valid(const Instruction& instr) const;

    // Provide const access to the GridManager object for the grid that this plaquette lives on
    const GridManager& grid() const {return grid_;} 

private:
    unsigned int a_;
    unsigned int b_;
    unsigned int c_;
    unsigned int d_;
    unsigned int m_;
    unsigned int row_;
    unsigned int col_;
    char shape_; 
    char type_;
    GridManager& grid_;

    // Private member function that returns an lvalue allowing you to modify the qsite of a qubit
    unsigned int& mod_qsite(char qubit);
};

}
#endif //TISCC_PLAQUETTE_HPP