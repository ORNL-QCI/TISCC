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

/**
 * @brief Manages the Plaquettes of a surface code patch and offers functionality both to manipulate them and to generate circuits for a set of primitive patch operations. 
 * 
 * Designed to be compatible with grid-like hardware architectures where data and syndrome measurement qubits, even if mobile, have dedicated "home" qsites.
 * Its construction on the hardware grid assumes the definition for logical tile specified in the TISCC paper.
*/
class LogicalQubit {
public:

    /**
     * @brief Constructor for the LogicalQubit class.
     * 
     * Requests Plaquettes from GridManager and constructs other members such as the parity check matrix and syndrome extraction circuits.
     * @param dx The initial x code distance (width) of the surface code patch.
     * @param dz The initial z code distance (height) of the surface code patch.
     * @param row The hardware grid row specifying the top-left corner of the logical tile on which this patch lives.
     * @param col The hardware grid column specifying the top-left corner of the logical tile on which this patch lives.
     * @param grid The GridManager from which Plaquettes should be requested for this LogicalQubit.
     */
    explicit LogicalQubit(unsigned int dx, unsigned int dz, unsigned int row, unsigned int col, GridManager& grid);

    // Accessors

    /**
     * @brief Get current x code distance.
     * @return Current x code distance.
    */
    unsigned int get_dx() const {return dx_;}

    /**
     * @brief Get current z code distance.
     * @return Current z code distance.
    */
    unsigned int get_dz() const {return dz_;}

    /**
     * @brief Get initial x code distance.
     * @return Initial x code distance.
    */
    unsigned int get_dx_init() const {return dx_init_;}

    /**
     * @brief Get initial z code distance.
     * @return Initial z code distance.
    */
    unsigned int get_dz_init() const {return dz_init_;}

    /**
     * @brief Get hardware grid row specifying the top-left corner of the logical tile on which this patch lives.
     * @return Hardware grid row specifying the top-left corner of the logical tile on which this patch lives.
    */
    unsigned int get_row() const {return row_;}

    /**
     * @brief Get hardware grid column specifying the top-left corner of the logical tile on which this patch lives.
     * @return Hardware grid column specifying the top-left corner of the logical tile on which this patch lives.
    */
    unsigned int get_col() const {return col_;}

    /**
     * @brief Get parity check matrix.
     * @return Const reference to parity check matrix.
    */
    const std::vector<std::vector<bool>>& get_parity_check_matrix() const {return parity_check_matrix;}

    /**
     * @brief Get map from qsite (grid array index) to parity check matrix column.
     * 
     * The parity check matrix actually has twice as many columns. 
     * The returned column corresponds with the Z operator.
     * To obtain the column corresponding with the X operator, use lq.get_qsite_to_index()[qsite] + lq.get_qsite_to_index().size().
     * @return Const reference to map from qsite (grid array index) to parity check matrix column.
    */
    const std::map<unsigned int, unsigned int>& get_qsite_to_index() const {return qsite_to_index;}

    /**
     * @brief Get a vector that maps parity check matrix column to qsite.
     * @return Const reference to a vector that maps parity check matrix column to qsite.
    */
    const std::vector<unsigned int>& get_index_to_qsite() const {return index_to_qsite;}

    // Derived quantities and other accessors

    /** 
     * @brief Constructs and returns a set of qsites occupied by the plaquettes of this patch.
     * @param just_data_qubits Specifies whether only data qubits (not syndrome measurement qubits) should be included.
     * @return A set of qsites occupied by the plaquettes of patch.
    */
    std::set<unsigned int> occupied_sites(bool just_data_qubits) const;

    /**
     * @brief Returns a boolean vector representing the default-edge logical operator of given type.
     * @param type Specifies either the X logical operator ('X') or the Z logical operator ('Z').
     * @return A const reference to the boolean vector representing the default-edge logical operator of given type.
    */
    const std::vector<bool>& get_logical_operator_default_edge(char type) const;

    /**
     * @brief Constructs and returns a boolean vector representing the opposite-edge logical operator of given type.
     * @param type Specifies either the X logical operator ('X') or the Z logical operator ('Z').
     * @return A boolean vector representing the opposite-edge logical operator of given type.
    */
    std::vector<bool> get_logical_operator_opposite_edge(char type) const;

    /**
     * @brief Calls get_logical_operator_default_edge or get_logical_operator_opposite_edge depending on given edge_type.
     * @param type Specifies either the X logical operator ('X') or the Z logical operator ('Z').
     * @param edge_type Specifies either the default-edge logical operator ("default") or the opposite-edge logical operator ("opposite").
     * @return A boolean vector representing the logical operator of given type and edge_type.
    */
    std::vector<bool> get_logical_operator(char type, std::string_view edge_type) const;

    /**
     * @brief Returns a vector containing qsites representing deformations that have been made to the default-edge logical operator of given type.
     * 
     * Primarily used to track logical operator deformations resulting from corner movements.
     * In corner movements, newly-measured boundary stabilizers may anti-commute with logical operators and/or stabilizers, requiring the logical operators to be re-defined.
     * TISCC also updates logical operator definitions when newly measured stabilizers overlap with them in order that they retain a minimal Pauli weight.
     * Lastly, logical operators may be re-defined in the case that corner data qubits of a patch are measured out.
     * All of this information must be incorporated into a sign update of the logical operator following the corner movement operation. 
     * @param type Specifies either the X logical operator deformation ('X') or the Z logical operator deformation ('Z').
     * @return A vector containing qsites representing deformations that have been made to the default-edge logical operator of given type.
    */
    std::vector<unsigned int> get_logical_deformation_qsites(char type) const;

    /**
     * @brief Clears vector containing qsites representing deformations that have been made to the default-edge logical operator of given type.
     * @param type Specifies either the X logical operator deformation ('X') or the Z logical operator deformation ('Z').
    */
    void clear_logical_deformation_qsites(char type);

     /**
     * @brief Returns a vector containing the syndrome measurement qsites required to deform the default-edge logical operator of given type into the corresponding opposite-edge logical operator.
     * 
     * @param type Specifies either the X logical operator ('X') or the Z logical operator ('Z').
     * @return A vector containing the syndrome measurement qsites required to deform the default-edge logical operator of given type into the corresponding opposite-edge logical operator.
    */   
    std::vector<unsigned int> get_logical_deformation_between_edges(char type) const;

     /**
     * @brief Returns a vector containing the syndrome measurement qsites required to move the default-edge logical operator of given type by n columns or rows.
     * 
     * Only valid for operators with support all on one column or all on one row. 
     * @param type Specifies either the X logical operator ('X') or the Z logical operator ('Z').
     * @param n Specifies the number of rows or columsn to shift the operator by.
     * @param grid The GridManager object for reference.
     * @return A vector containing the syndrome measurement qsites required to move the default-edge logical operator of given type by n columns or rows.
    */   
    std::vector<unsigned int> get_logical_deformation_operator_movement(char type, int n, const GridManager& grid);

    /**
     * @brief Returns True if the patch is in a canonical stabilizer arrangement (see Fig. 2 of TISCC paper), and False otherwise.
     * @return True if the patch is in a canonical stabilizer arrangement (see Fig. 2 of TISCC paper), and False otherwise.
    */
    bool canonical_arrangement() const {return canonical_arrangement_;}

    /**
     * @brief Returns True if xz_swap has been applied to this patch and False otherwise.
     * @return True if xz_swap has been applied to this patch and False otherwise.
    */
    bool xz_swap_tracker() const {return xz_swap_tracker_;}

    /**
     * @brief Returns True if flip_patch has been applied to this patch and False otherwise.
     * @return True if flip_patch has been applied to this patch and False otherwise.
    */
    bool flipped_tracker() const {return flipped_tracker_;}

    // Print functions

    /**
     * @brief Prints out some details of the stabilizers stored by this patch. 
    */
    void print_stabilizers() const;

    /**
     * @brief Prints out the parity check matrix stored by this patch and some ancillary data.
     * 
     * Also checks and returns validity of parity check matrix.
    */
    void print_parity_check_matrix() const;

    // Operations (see TISCC paper Table 2 and surrounding text for details).

    /**
     * @brief Generates the hardware circuit for a logical Idle operation (which is a specified number of cycles of syndrome extraction).
     * 
     * @param cycles The number of cycles of syndrome extraction.
     * @param grid The GridManager object for reference.
     * @param hw_master The vector of hardware instructions to append this operation to.
     * @param time The nominal time at which this operation begins.
     * @return The nominal time at which this operation ends.
    */
    double idle(unsigned int cycles, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);

    /**
     * @brief Generates the hardware circuit for a given transversal operation (acts on all data qubits of the patch).
     * 
     * @param op Operation label. Options: "prepz", "prepx", "measz", "measx", and "hadamard".
     * @param grid The GridManager object for reference.
     * @param hw_master The vector of hardware instructions to append this operation to.
     * @param time The nominal time at which this operation begins.
     * @return The nominal time at which this operation ends.
    */
    double transversal_op(const std::string& op, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);

    /**
     * @brief Generates the hardware circuit for the application of a given logical Pauli operator.
     * 
     * @param pauli Operator label. Options: 'X', 'Y', 'Z'.
     * @param grid The GridManager object for reference.
     * @param hw_master The vector of hardware instructions to append this operation to.
     * @param time The nominal time at which this operation begins.
     * @return The nominal time at which this operation ends.
    */
    double apply_pauli(char pauli, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);

    /**
     * @brief Generates a hardware circuit preparing the data qubits for state injection (the state is encoded upon a subsequent Idle operation).
     * 
     * Prepares the top-left corner qubit in the appropriate single-qubit state depending on the given label.
     * Otherwise, prepares data qubits in the |+> state along the default-edge logical X operator and |0> everywhere else.
     * Does not employ post-selection and patch expansion techniques of protocols such as the one proposed in https://arxiv.org/abs/1410.7808.
     * @param label State label. Options: 'y', 't'.
     * @param grid The GridManager object for reference.
     * @param hw_master The vector of hardware instructions to append this operation to.
     * @param time The nominal time at which this operation begins.
     * @return The nominal time at which this operation ends.
    */
    double inject_state(char label, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);

    /**
     * @brief Generates the hardware circuit for the Swap Left operation (see Fig. 4 of the TISCC paper).
     * 
     * @param grid The GridManager object for reference.
     * @param hw_master The vector of hardware instructions to append this operation to.
     * @param time The nominal time at which this operation begins.
     * @return The nominal time at which this operation ends.
    */
    double swap_left(GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) {time = translate_patch(0, -1, grid, hw_master, time); return time;}

    /**
     * @brief Generates the hardware circuit for the Move Right operation (see Fig. 4 of the TISCC paper).
     * 
     * @param cycles The number of cycles of syndrome extraction.
     * @param lq_extended Exchanges a nullptr for a pointer to the LogicalQubit corresponding with the intermediate extended patch.
     * @param lq_contracted Exchanges a nullptr for a pointer to the LogicalQubit corresponding with the final contracted patch.
     * @param grid The GridManager object for reference.
     * @param hw_master The vector of hardware instructions to append this operation to.
     * @param time The nominal time at which this operation begins.
     * @return The nominal time at which this operation ends.
    */
    double move_right(unsigned int cycles, LogicalQubit*& lq_extended, LogicalQubit*& lq_contracted, GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);
    
    /**
     * @brief Construct and return a logical qubit that represents the merged product of this qubit with an input one.
     * 
     * Currently limited to patches in the standard arrangement.
     * lq2 must be oriented either one tile east-ward or one tile south-ward from the present LogicalQubit.
     * @param lq2 LogicalQubit object to be merged with this one.
     * @param grid The GridManager object to use in constructing the merged patch.
     * @return A pointer to the resultant merged patch. 
    */
    LogicalQubit* get_merged_lq(LogicalQubit& lq2, GridManager& grid);

    /**
     * @brief Generates the hardware circuit for a merge operation.
     * 
     * This is only a valid operation on LogicalQubit objects that had been returned by get_merged_lq.
     * @param cycles The number of cycles of syndrome extraction.
     * @param grid The GridManager object for reference.
     * @param hw_master The vector of hardware instructions to append this operation to.
     * @param time The nominal time at which this operation begins.
     * @return The nominal time at which this operation ends.
    */
    double merge(unsigned int cycles, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);

    /**
     * @brief Generates the hardware circuit for a split operation.
     * 
     * This is only a valid operation on LogicalQubit objects that had been returned by get_merged_lq.
     * @param grid The GridManager object for reference.
     * @param hw_master The vector of hardware instructions to append this operation to.
     * @param time The nominal time at which this operation begins.
     * @return The nominal time at which this operation ends.
    */
    double split(GridManager& grid, std::vector<HW_Instruction>& hw_master, double time);

    /**
     * @brief Constructs and returns the set of data qsites from this patch that are not occupied by two other input patches.
     * 
     * Typically used to obtain the data qsites between two patches involved in a merge or split operation.
     * @param lq1 
     * @param lq2
     * @return The set of data qsites from this patch that are not occupied by two other input patches. 
    */
    std::set<unsigned int> get_strip(LogicalQubit& lq1, LogicalQubit& lq2);

    /**
     * @brief Swap roles of X and Z operators for this patch.
     * 
     * Typically used to transform operators following a transversal Hadamard or to obtain a LogicalQubit object with a rotated stabilizer arrangement.
    */
    void xz_swap(const GridManager& grid);

    /**
     * @brief Adds a new boundary stabilizer plaquette and handles the consequences of doing so (typically used in corner movement).
     * 
     * Removes anti-commuting stabilizers, updates logical operators, and measures/prepares corner qubits as needed.
     * Checks validity of the final parity check matrix.
     * A subsequent Idle operation is needed to actually measure the newly added stabilizer.
     * @param row The hardware grid row at which the new stabilizer plaquette is "pinned".
     * @param col The hardware grid col at which the new stabilizer plaquette is "pinned".
     * @param shape The shape of the added stabilizer.
     * @param type The type of the added stabilizer.
     * @param grid The GridManager object to obtain the new stabilizer from.
     * @param hw_master The vector of hardware instructions to append any needed hardware operations to.
     * @param time The nominal time at which this operation begins.
     * @param debug A boolean to trigger debugging output.
     * @return The nominal time at which this operation ends.
    */
    double add_stabilizer(unsigned int row, unsigned int col, char shape, char type, GridManager& grid, std::vector<HW_Instruction>& hw_master, double time, bool debug);

    /**
     * @brief Performs a corner movement by adding a specified weight to a specified logical operator.
     * 
     * Figures out which stabilizers to measure when extending a logical operator `clockwise' and calls add_stabilizer until the specified weight has been added.
     * Warning: add_stabilizer will update the default-edge (rather than opposite-edge) logical operator that has support on the added stabilizer in ambiguous cases.  
     * This may cause unexpected results where the operator to be extended has support on a single qubit.
     * @param type The type of logical operator to extend ('X' or 'Z').
     * @param edge_type The edge_type of logical operator to extend ("default" or "opposite").
     * @param weight_to_add The (minimum) weight to add in extending the operator (will add stabilizers until the final added weight is at least this amount)
     * @param stop_at_patch_corner Causes extension to stop if a patch corner is hit.
     * @param grid The GridManager object to obtain new stabilizers from.
     * @param hw_master The vector of hardware instructions to append any needed hardware operations to.
     * @param time The nominal time at which this operation begins.
     * @param debug A boolean to trigger debugging output.
     * @return The nominal time at which this operation ends.
    */
    double extend_logical_operator_clockwise(char type, std::string_view edge_type, unsigned int weight_to_add, bool stop_at_patch_corner, GridManager& grid, std::vector<HW_Instruction>& hw_master, double time, bool debug);

    /**
     * @brief A series of corner movements with the resulting stabilizer arrangement the same as if we flipped the patch upside down and then applied xz_swap.
     * 
     * See Fig. 3 of the TISCC paper for details.
     * Only implemented for patches in the standard and rotated arrangements.
     * Works for all odd code distances and even code distances >= 6.
     * Needn't be used as a surface code operation in itself; can be used to transform the stabilizers of this LogicalQubit to a flipped or rotated-flipped arrangement.
     * @param grid The GridManager object to obtain new stabilizers from in corner movements.
     * @param hw_master The vector of hardware instructions to append this operation to (depending on compile_ops).
     * @param time The nominal time at which this operation begins.
     * @param compile_ops A boolean specifying whether the hardware operations involved in Flip Patch should be added to the circuit.
     * @param debug A boolean to trigger debugging output.
     * @return The nominal time at which this operation ends.
    */
    double flip_patch(GridManager& grid, std::vector<HW_Instruction>& hw_master, float time, bool compile_ops, bool debug);

    /**
     * @brief Resets stabilizer circuit patterns to default values.
    */ 
    void reset_stabilizer_circuit_patterns();

    /**
     * @brief Toggles all stabilizer circuit patterns.
    */ 
    void swap_stabilizer_circuit_patterns();

    /* Other helpful functions */

    /**
     * @brief Transform operators from binary representation to pair<qsite unsigned int, Pauli char>
     * 
     * Note: this function will convert ZX directly to Y without tracking any phase.
    */
    std::vector<std::pair<unsigned int, char>> binary_operator_to_qsites(const std::vector<bool>& binary_rep);

    /**
     * @brief Obtain pair<qsite unsigned int, Pauli char> for each stabilizer measure qubit for the purpose of labeling on the grid.
     */
    std::vector<std::pair<unsigned int, char>> syndrome_measurement_qsites();

// Helper functions related to binary vector math

    /**
     * @brief Calculates the dot product mod 2 of two binary vectors.
     * @return A Boolean variable representing the result.
    */
    static bool bin_dot_prod_mod_2(const std::vector<bool>& v1, const std::vector<bool> v2);

    /**
     * @brief Calculates the product of two operators expressed in binary symplectic format.
     * @return A pair with the binary symplectic vector on one hand and a sign on the other.
    */
    static std::pair<std::vector<bool>, int> operator_product_binary_format(const std::vector<bool>& v1, const std::vector<bool> v2);

    /**
     * @brief Symplectic transform of a binary symplectic vector.
     * @return The transformed vector.
    */
    static std::vector<bool> symplectic_transform(const std::vector<bool>& v);

    /**
     * @brief Returns the Hamming weight of a given binary vector.
     * @return The Hamming weight of a given binary vector.
    */
    static unsigned int hamming_weight(const std::vector<bool>& v);

    /**
     * @brief Transform operators from binary symplectic to string representation while tracking phase e.g. 11000101 = Z(ZX)IX = i*ZYIX.
     * @return A pair with the operator string on one hand and a phase on the other.
    */
    static std::pair<std::string, std::complex<double>> binary_operator_to_pauli_string(const std::vector<bool>& binary_rep);

    /**
     * @brief Transform operators from string to binary symplectic representation, tracking phase e.g. ZYIX -> -i*(11000101)
     * @return A pair with the binary symplectic vector on one hand and a phase on the other.
    */
    static std::pair<std::vector<bool>, std::complex<double>> pauli_string_to_binary_operator(const std::string& pauli_string);

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

    // If this is a merged product, we employ pointers to track which lq were merged to produce it
    LogicalQubit* lq1;
    LogicalQubit* lq2;

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

}

#endif //TISCC_LOGICALQUBIT_HPP