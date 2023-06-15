#ifndef TISCC_LOGICALQUBIT_HPP
#define TISCC_LOGICALQUBIT_HPP

#include <TISCC/hardwaremodel.hpp>
#include <TISCC/plaquette.hpp>
#include <TISCC/gridmanager.hpp>
#include <TISCC/instruction.hpp>

#include <vector>
#include <set>
#include <map>
#include <optional>

namespace TISCC {

// Consists of the stabilizer plaquettes of a surface code patch and defined operations on them
class LogicalQubit {
public:
    // Constructor
    explicit LogicalQubit(unsigned int dx, unsigned int dz, unsigned int row, unsigned int col, GridManager& grid);

    // Accessors
    unsigned int get_dx() const {return dx_;}
    unsigned int get_dz() const {return dz_;}
    unsigned int get_row() const {return row_;}
    unsigned int get_col() const {return col_;}

    // Operations
    double idle(unsigned int cycles, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);
    double transversal_op(const std::string& op, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);

    // Placeholder function to help implement little test circuits
    double test_circuits(const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);

    // Print functions
    void print_stabilizers();
    void print_parity_check_matrix(const GridManager& grid);

    // Function to output all qsites occupied by the surface code
    std::set<unsigned int> occupied_sites();

    // Function to output all qsites occupied by data qubits on the surface code
    std::set<unsigned int> data_qsites();

    // Helper function to return the data qubits from this patch that are NOT occupied by two others
    std::set<unsigned int> get_strip(LogicalQubit& lq1, LogicalQubit& lq2);

    // Swap roles of x and z for this patch (used during Hadamard and patch rotation)
    void xz_swap();

    // Add new stabilizer plaquette (used in corner movement)
    void add_stabilizer(unsigned int row, unsigned int col, char shape, char type, GridManager& grid, bool debug);

private:
    // Code distances
    unsigned int dx_;
    unsigned int dz_;
    unsigned int row_;
    unsigned int col_;

    // Vectors of X or Z plaquettes 
    std::vector<Plaquette> x_plaquettes;
    std::vector<Plaquette> z_plaquettes;

    // Vectors of instructions to hold syndrome measurement circuits
    std::vector<Instruction> Z_Circuit_Z_Type;
    std::vector<Instruction> X_Circuit_N_Type;
    std::vector<Instruction> Z_Circuit_N_Type;
    std::vector<Instruction> X_Circuit_Z_Type;

    // Boolean matrix to store parity check matrix
    std::optional<std::vector<std::vector<bool>>> parity_check_matrix;

    // Map from qsite to column index in parity check matrix
    std::optional<std::map<unsigned int, unsigned int>> qsite_to_index;

    // Inverse map to above
    std::optional<std::vector<unsigned int>> index_to_qsite;

    // Transform operators from binary representation to pair<qsite unsigned int, Pauli char>
    std::vector<std::pair<unsigned int, char>> binary_operator_to_qsites(const std::vector<bool>& binary_rep);

    // Contains details of hardware native gates and stabilizer circuits
    HardwareModel TI_model;

    // Construct stabilizers and update set of occupied sites in grid
    void init_stabilizers(unsigned int dx, unsigned int dz, unsigned int row, unsigned int col, GridManager& grid); 

    // Test stabilizers
    void test_stabilizers();

    // Set up the circuits that we intend to use
    void init_circuits();

    // Construct parity check matrix from stabilizers
    void construct_parity_check_matrix(const GridManager& grid);

    // Check validity of parity check matrix
    bool validity_parity_check_matrix();

    // Apply a given instruction to all plaquettes in a given vector
    double apply_instruction(const Instruction& instr, std::vector<Plaquette>& plaquettes, double time, unsigned int step,
        const GridManager& grid, std::vector<HW_Instruction>& idle_operation);

};

// Construct and return a logical qubit that represents the merged product of two input logical qubits
LogicalQubit merge(LogicalQubit& lq1, LogicalQubit& lq2, GridManager& grid);

}

#endif //TISCC_LOGICALQUBIT_HPP