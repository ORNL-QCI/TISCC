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

    // Accessor functions
    unsigned int get_qsite(char qubit) const;
    unsigned int get_row() const {return row_;}
    unsigned int get_col() const {return col_;}
    char get_shape() const {return shape_;}
    char get_operator_type() const {return operator_type_;}
    char get_circuit_pattern() const {return circuit_pattern_;}

    // Modifier functions
    void change_operator_type(char type) {operator_type_ = type;}
    void set_circuit_pattern(char circuit_pattern) {circuit_pattern_ = circuit_pattern;}
    
    // Move qubit to a specified site after checking validity and recording change on grid
    // **Note: Be careful. Moving one of these plaquette's qubits here does not update its reference in other plaquettes.
    void move_to_site(char q, unsigned int site);

    // Check whether a qubit is currently at its home location on the grid
    bool is_home(char qubit) const;

    // Ask whether the plaquette has support on the site in question and remove if so
    bool remove_supported_qsite(unsigned int site);

   // Check to see if a given instruction is valid on this plaquette 
    bool is_instr_valid(const Instruction& instr) const;

    // Provide const access to the GridManager object for the grid that this plaquette lives on
    GridManager* grid() const {return grid_;} 

private:
    // Plaquette constructor made only accessible to GridManager
    friend class GridManager;
    explicit Plaquette(unsigned int a, unsigned int b, unsigned int c, unsigned int d, unsigned int m, unsigned int row,
        unsigned int col, char shape, char type, GridManager* grid) :
        a_(a), b_(b), c_(c), d_(d), m_(m), row_(row), col_(col), shape_(shape), operator_type_(type), grid_(grid) {}

    unsigned int a_;
    unsigned int b_;
    unsigned int c_;
    unsigned int d_;
    unsigned int m_;
    unsigned int row_;
    unsigned int col_;
    char shape_; 
    char operator_type_;
    char circuit_pattern_;
    GridManager* grid_;

    // Private member function that returns an lvalue allowing you to modify the qsite of a qubit
    unsigned int& mod_qsite(char qubit);
};

}
#endif //TISCC_PLAQUETTE_HPP