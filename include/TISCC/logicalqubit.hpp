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
#include <complex> 

namespace TISCC {

/* Consists of the stabilizer plaquettes of a surface code patch, functions to manipulate them, and functions generate circuits for primitive patch operations.
   Designed to be compatible with grid-like hardware architectures. */
class LogicalQubit {
public:

    // Constructor
    explicit LogicalQubit(unsigned int dx, unsigned int dz, unsigned int row, unsigned int col, GridManager& grid);

    // Accessors
    unsigned int get_dx() const {return dx_;}
    unsigned int get_dz() const {return dz_;}
    unsigned int get_dx_init() const {return dx_init_;}
    unsigned int get_dz_init() const {return dz_init_;}
    unsigned int get_row() const {return row_;}
    unsigned int get_col() const {return col_;}
    const std::vector<std::vector<bool>>& get_parity_check_matrix() const {return parity_check_matrix;}
    const std::map<unsigned int, unsigned int>& get_qsite_to_index() const {return qsite_to_index;}
    const std::vector<unsigned int>& get_index_to_qsite() const {return index_to_qsite;}

    // Derived quantities and other accessors
    std::set<unsigned int> occupied_sites(bool just_data_qubits) const;
    const std::vector<bool>& get_logical_operator_default_edge(char type) const;
    std::vector<bool> get_logical_operator_opposite_edge(char type) const;
    std::vector<bool> get_logical_operator(char type, std::string_view edge_type) const;
    std::vector<unsigned int> get_logical_deformation_qsites(char type) const;
    std::vector<unsigned int> get_logical_deformation_between_edges(char type) const;
    std::vector<unsigned int> get_logical_deformation_operator_movement(char type, int n, const GridManager& grid);
    bool canonical_arrangement() const {return canonical_arrangement_;}
    bool xz_swap_tracker() const {return xz_swap_tracker_;}
    bool flipped_tracker() const {return flipped_tracker_;}

    // Print functions
    void print_stabilizers() const;
    void print_parity_check_matrix() const;

    // Operations
    double idle(unsigned int cycles, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);
    double transversal_op(const std::string& op, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);
    double apply_pauli(char pauli, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);
    double inject_state(char label, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);
    double swap_left(GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) {time = translate_patch(0, -1, grid, hw_master, time); return time;}
    
    // merge: construct and return a logical qubit that represents the merged product of this qubit with an input one
    LogicalQubit* merge(LogicalQubit& lq2, GridManager& grid);

    // A series of corner movements with the resulting strabilizer arrangement the same as if we flipped the patch upside down and then applied xz_swap
    float flip_patch(GridManager& grid, std::vector<HW_Instruction>& hw_master, float time, bool compile_ops, bool debug);

    // Function to return the data qubits from this patch that are NOT occupied by two others (typically used to get the qsites on the intervening strip between lq from a merged product)
    std::set<unsigned int> get_strip(LogicalQubit& lq1, LogicalQubit& lq2);

    // Swap roles of x and z for this patch (typically used during Hadamard and patch rotation)
    void xz_swap(const GridManager& grid);

    /* Corner-movement related operations */

    // Add new stabilizer plaquette (typically used in corner movement)
    double add_stabilizer(unsigned int row, unsigned int col, char shape, char type, GridManager& grid, std::vector<HW_Instruction>& hw_master, double time, bool debug);

    // Corner movement 
    double extend_logical_operator_clockwise(char type, std::string_view edge_type, unsigned int weight_to_add, bool stop_at_patch_corner, GridManager& grid, std::vector<HW_Instruction>& hw_master, double time, bool debug);

    // Reset stabilizer circuits to default values
    void reset_stabilizer_circuit_patterns();

    // Swap stabilizer circuits 
    void swap_stabilizer_circuit_patterns();

    /* Other helpful functions */

    // Transform operators from binary representation to pair<qsite unsigned int, Pauli char>
    std::vector<std::pair<unsigned int, char>> binary_operator_to_qsites(const std::vector<bool>& binary_rep);

    // Obtain pair<qsite unsigned int, Pauli char> for each stabilizer measure qubit for the purpose of labeling on the grid
    std::vector<std::pair<unsigned int, char>> syndrome_measurement_qsites();

private:
    // Code distances and location on grid
    unsigned int dx_;
    unsigned int dx_init_;
    unsigned int dz_;
    unsigned int dz_init_;
    unsigned int row_;
    unsigned int col_;

    // Track whether deformations have taken place
    bool canonical_arrangement_;

    // Boolean to toggle each time xz_swap is applied
    bool xz_swap_tracker_;

    // Boolean to toggle each time flip_patch is applied
    bool flipped_tracker_;

    // Contains details of hardware native gates and stabilizer circuits
    HardwareModel TI_model;

    // Vectors of X or Z plaquettes 
    std::vector<Plaquette> x_plaquettes;
    std::vector<Plaquette> z_plaquettes;

    // Boolean matrix to store parity check matrix
    std::vector<std::vector<bool>> parity_check_matrix;

    // Map from qsite to column index in parity check matrix
    std::map<unsigned int, unsigned int> qsite_to_index;
    std::vector<unsigned int> index_to_qsite;

    // Vectors of instructions to hold syndrome measurement circuits
    std::vector<Instruction> Z_Circuit_Z_Type;
    std::vector<Instruction> X_Circuit_N_Type;
    std::vector<Instruction> Z_Circuit_N_Type;
    std::vector<Instruction> X_Circuit_Z_Type;

    // We track stabilizer measurement qsites corresponding with logical operator deformation
    std::vector<unsigned int> x_deformation_qsites;
    std::vector<unsigned int> z_deformation_qsites;

    // Construct stabilizers and update set of occupied sites in grid
    void init_stabilizers(unsigned int dx, unsigned int dz, unsigned int row, unsigned int col, GridManager& grid); 

    // Test stabilizers
    void test_stabilizers();

    // Set up the circuits that we intend to use
    void init_circuits();

    // Construct parity check matrix from stabilizers
    void construct_parity_check_matrix(const GridManager& grid);

    // Check validity of parity check matrix
    bool validity_parity_check_matrix() const;

    // Recalculate code distance using logical operators
    void recalculate_code_distance(); 

    // Update the logical deformation vectors with qsite
    void add_logical_deformation_qsites(char type, unsigned int qsite);

    // Apply a given instruction to all plaquettes in a given vector
    double apply_instruction(const Instruction& instr, Plaquette& p, double time, unsigned int step,
        const GridManager& grid, std::vector<HW_Instruction>& hw_master);

    // Translate patch s rows "South" or e rows "East" on the underlying grid (not fully implemented)
    double translate_patch(int s, int e, GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);

};

// Helper functions related to binary vector math
bool bin_dot_prod_mod_2(const std::vector<bool>& v1, const std::vector<bool> v2);
std::pair<std::vector<bool>, int> operator_product_binary_format(const std::vector<bool>& v1, const std::vector<bool> v2);
std::vector<bool> symplectic_transform(const std::vector<bool>& v);
unsigned int pauli_weight(const std::vector<bool>& v);
std::pair<std::string, std::complex<double>> binary_operator_to_pauli_string(const std::vector<bool>& binary_rep);
std::pair<std::vector<bool>, std::complex<double>> pauli_string_to_binary_operator(const std::string& pauli_string);

}

#endif //TISCC_LOGICALQUBIT_HPP