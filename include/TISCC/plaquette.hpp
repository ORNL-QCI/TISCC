#ifndef TISCC_PLAQUETTE_HPP
#define TISCC_PLAQUETTE_HPP

#include <TISCC/gridmanager.hpp>
#include <TISCC/instruction.hpp>

#include <iostream>
#include <vector>
#include <set>

namespace TISCC {

// Forward declaration of dependencies
class GridManager;
class Instruction;

/**
 * @brief Stores information about a surface code plaquette.
 * 
 * The Plaquette constructor is private. New Plaquettes are typically requested through GridManager::get_plaquette
 * from within a LogicalQubit member function (such as LogicalQubit::LogicalQubit or LogicalQubit::add_stabilizer).
 */
class Plaquette {
public:

    // Accessor functions

    /**
     * @brief Get the qsite (grid array index) corresponding to a qubit label.
     * 
     * Valid qubit labels are 'a', 'b', 'c', 'd', and 'm'.
     * @return The qsite (grid array index) for the given qubit label.
    */
    unsigned int get_qsite(char qubit) const;

    /**
     * @brief Get the hardware grid row on which this plaquette is "pinned".
     * @return The hardware grid row on which this plaquette is "pinned".
    */
    unsigned int get_row() const {return row_;}

    /**
     * @brief Get the hardware grid column on which this plaquette is "pinned".
     * @return The hardware grid column on which this plaquette is "pinned".
    */
    unsigned int get_col() const {return col_;}

    /**
     * @brief Get the plaquette "shape".
     * 
     * The shape is a character specifying either a four-qubit plaquette ('f') or a two-qubit plaquette with edge type {'n','s','e','w'}.
     * @return The plaquette "shape".
    */
    char get_shape() const {return shape_;}

    /**
     * @brief Get the plaquette type (X or Z).
     * @return The plaquette type (X or Z).
    */
    char get_operator_type() const {return operator_type_;}

    /**
     * @brief Get the circuit pattern for this plaquette (Z or N).
     * @return The circuit pattern (Z or N).
    */
    char get_circuit_pattern() const {return circuit_pattern_;}

private:
    // Declare friend classes
    friend class LogicalQubit;
    friend class HardwareModel;
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

    /**
     * @brief Provide const access to the GridManager object from which this Plaquette was initially requested.
     * @return A pointer to the GridManager object.
    */
    const GridManager* grid() const {return grid_;} 

    /**
     * @brief Check whether a qubit is currently at its home location on the grid
     * 
     * This is primarily used for book-keeping within the CNOT gate compilation method within HardwareModel.
     * Depends on the details of the repeating unit being used by GridManager and the hardware mapping of a surface code plaquette.
     * @return True if the qubit is at its home location on the grid.
    */
    bool is_home(char qubit) const;

    /**
     * @brief Move qubit to a specified site after checking validity and recording change on grid.
     * 
     * This is primarily used for book-keeping within the CNOT gate compilation method within HardwareModel.
    */
    void move_to_site(char q, unsigned int site);

    /**
     * @brief Returns an lvalue allowing you to modify the qsite of a qubit.
     * 
     * Primarily used within move_to_site.
     * @return An lvalue allowing you to modify the qsite of a qubit.
     */
    unsigned int& mod_qsite(char qubit);

    /**
     * @brief Modify the operator type (X or Z).
     * 
     * Primarily used within LogicalQubit::xz_swap.
    */
    void change_operator_type(char type) {operator_type_ = type;}

    /**
     * @brief Modify the circuit pattern (Z or N).
     * 
     * Primarily used in LogicalQubit member functions.
    */
    void set_circuit_pattern(char circuit_pattern) {circuit_pattern_ = circuit_pattern;}

    /**
     * @brief Ask whether the plaquette has support on the site in question and remove qubit if so.
     * 
     * Primarily used within LogicalQubit::add_stabilizer when a corner qubit needs to be removed, leaving a three-qubit plaquette.
     * @return True if a qubit was removed.
    */
    bool remove_supported_qsite(unsigned int site);

    /**
     * @brief Check to see if a given instruction is valid on this plaquette given its shape.
     * @return True if the instruction is valid.
    */
    bool is_instr_valid(const Instruction& instr) const;

};

}
#endif //TISCC_PLAQUETTE_HPP