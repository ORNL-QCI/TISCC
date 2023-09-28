#ifndef TISCC_HARDWAREMODEL_HPP
#define TISCC_HARDWAREMODEL_HPP

#include <TISCC/instruction.hpp>
#include <TISCC/gridmanager.hpp>
#include <TISCC/plaquette.hpp>

#include<map>
#include<string>
#include<vector>
#include<set>

namespace TISCC {

// Forward declaration of dependencies
class Plaquette;
class GridManager;
class HW_Instruction;

/**
 * @brief Defines a set of native hardware operations and provides an interface for constructing circuits out of them.
 */
class HardwareModel {
public:
    /**
     * @brief Constructs a HardwareModel object with default trap width and initializes a map of trapped-ion operations to their default durations.
     * 
     * Trap width, available operations, and their durations are hard-coded but can be straightforwardly adapted within the code.
     */
    explicit HardwareModel();

    /**
     * @brief Returns a map containing the durations for implemented native trapped-ion gates.
     * @return A const reference to the map containing the durations for implemented native trapped-ion gates.
     */
    const std::map<std::string, double>& get_ops() const {return TI_ops;}

    /* Functions for adding native trapped-ion gates to a master hardware circuit */

    // The below are meant to be used within the context of a syndrome extraction circuit (i.e. give a plaquette and qubit label)

    /**
     * @brief Adds a Prepare_Z operation to the circuit for a specific qubit label on a plaquette.
     *
     * @param p The plaquette on which the operation is applied.
     * @param qubit The qubit label.
     * @param time The current time in microseconds.
     * @param step The step of the syndrome extraction circuit.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */
    double add_init(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const;
    
    /**
     * @brief Adds a Measure_Z operation to the circuit for a specific qubit label on a plaquette.
     *
     * @param p The plaquette on which the operation is applied.
     * @param qubit The qubit label.
     * @param time The current time in microseconds.
     * @param step The step of the syndrome extraction circuit.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */
    double add_meas(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const;
    
    /**
     * @brief Compiles a Hadamard gate into hardware instructions for a specific qubit label on a plaquette.
     *
     * @param p The plaquette on which the operation is applied.
     * @param qubit The qubit label.
     * @param time The current time in microseconds.
     * @param step The step of the syndrome extraction circuit.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */
    double add_H(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Compiles a CNOT gate into hardware instructions for a specific qubit label on a plaquette.
     *
     * @param p The plaquette on which the operation is applied.
     * @param control The control qubit label.
     * @param target The target qubit label.
     * @param time The current time in microseconds.
     * @param step The step of the syndrome extraction circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */
    double add_CNOT(Plaquette& p, char control, char target, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    // The below are meant to be used in general (i.e. give a qsite)

    /**
     * @brief Adds a Prepare_Z operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */
    double add_init(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Adds a Measure_Z operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */
    double add_meas(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Compiles a Hadamard gate into hardware instructions for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */
    double add_H(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Adds an X_pi/2 operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */
    double add_X(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Adds an X_pi/4 operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */
    double add_X_pi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Adds an X_-pi/4 operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */
    double add_X_mpi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const; 

    /**
     * @brief Adds a Y_pi/2 operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */  
    double add_Y(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Adds a Y_pi/4 operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */ 
    double add_Y_pi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Adds a Y_-pi/4 operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */ 
    double add_Y_mpi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Adds a Z_pi/2 operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */ 
    double add_Z(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Adds a Z_pi/4 operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */ 
    double add_Z_pi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Adds a Z_-pi/4 operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */ 
    double add_Z_mpi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Adds a Z_pi/8 operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */ 
    double add_Z_pi8(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Adds a Z_-pi/8 operation to the circuit for a specific site.
     *
     * @param site The site on which the operation is applied.
     * @param time The current time in microseconds.
     * @param step The step of the sub-circuit.
     * @param grid The GridManager object for reference.
     * @param circuit The hardware circuit to which the operation is added.
     * @return The updated time after adding the operation.
     */ 
    double add_Z_mpi8(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    /**
     * @brief Shifts a group of data qubits one column left-ward on the hardware grid.
     *
     * This function is responsible for moving a set of qubits one column left-ward on the hardware grid.
     * It depends on the details of our chosen hardware architecture specified in GridManager.
     * It assumes that all input qubits are on 'O' qsites and that there are no occupied 'M' sites along the relevant paths.
     * It first moves the input qubits out of the way (adjacent to the south-west 'O' site), then shuttles the strip 
     * of qsites (if occupied) from the left-ward column to the other end of the tile.
     * Finally, it moves the input qubits to their final positions.
     *
     * @param qsites A set of qubits to be shifted left-ward.
     * @param patch_cols The number of columns in the surface code patch.
     * @param grid The hardware grid on which to perform qubit movements.
     * @param hw_master The master circuit to which the Move operations are added.
     * @param time The current time in microseconds.
     * @return The updated time after performing the Move operations.
     */
    double shift_left(const std::set<unsigned int>& qubits, unsigned int patch_cols, GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) const;

    /**
     * @brief Prints map from trapped-ion operations to their default durations.
     * @return Map from trapped-ion operations to their default durations.
    */
    void print_TI_ops() const;

    /**
     * @brief Returns the default trap width.
     * @return Default trap width.
    */
    double get_trap_width() const {return trap_width;}

    /**
     * @brief Returns the default cell width.
     * 
     * It depends on the details of our chosen hardware architecture specified in GridManager.
     * @return Default cell width.
    */
    double get_cell_width() const {return cell_width;}

private:
    // Parameters of the hardware model
    double trap_width; // us
    double cell_width; // us

    // Hash table to map trapped-ion instructions to time (in microseconds)
    std::map<std::string, double> TI_ops;

    // Initialize hash table to define trapped-ion instruction set and map instructions to time (in microseconds)
    void init_TI_ops();

    // HW circuit helper functions to compile Moves
    void move_along_path(Plaquette& p, unsigned int step, std::vector<HW_Instruction>& circuit, double& time,
        const std::vector<unsigned int>& path, const GridManager& grid) const;
    double move_along_path(unsigned int qsite, unsigned int step, std::vector<HW_Instruction>& circuit, double time,
        const std::vector<unsigned int>& path, GridManager& grid) const;
};

}
#endif //TISCC_HARDWAREMODEL_HPP