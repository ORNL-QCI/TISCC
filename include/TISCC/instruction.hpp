#ifndef TISCC_INSTRUCTION_HPP
#define TISCC_INSTRUCTION_HPP

#include<TISCC/gridmanager.hpp>
#include<string>

namespace TISCC {

// Forward declaration of dependencies
class GridManager;

/**
 * @brief Represents a named instruction acting on the qubit labels (a, b, c, d, or m) of a general surface code plaquette.
 * 
 * Vectors of Instructions are currently used by LogicalQubits to apply parallel layers of gates to all their plaquettes.
 * These are intended to be QASM-like instructions that act on qubit labels rather than qsites.
 * There are no guard-rails as to acceptable instruction names or qubit labels. 
 */
class Instruction {
public:
    /**
     * @brief Constructor for the Instruction class.
     * @param name The name of the instruction.
     * @param q1 The first qubit label.
     * @param q2 The second qubit label.
     */
    explicit Instruction(std::string name, char q1, char q2) : name_(name), q1_(q1), q2_(q2) {};
    /**
     * @brief Gets the name of the instruction.
     * @return A const reference to the instruction name.
     */
    const std::string& get_name() const {return name_;}
    /**
     * @brief Gets the first qubit label.
     * @return The first qubit label as a character.
     */
    char get_q1() const {return q1_;}
    /**
     * @brief Gets the second qubit label.
     * @return The second qubit label as a character.
     */
    char get_q2() const {return q2_;}
private:
    std::string name_;
    char q1_;
    char q2_;
};

/**
 * @brief Represents a named instruction acting at a particular qsite(s) on the hardware grid at a specific time.
 * 
 * Instructions are compiled into HW_Instructions by LogicalQubit using the member functions of HardwareModel.
 * The names are intended to be members of a hardware gate set defined within HardwareModel.
 * HW_Instructions act on qsites rather than qubits.
 * There are no guard-rails as to acceptable instruction names or grid array indices (qsites). 
 */
class HW_Instruction {
public:
    /**
     * @brief Primary constructor for the HW_Instruction class.
     * @param name The name of the hardware instruction.
     * @param site1 The first target qsite on the grid.
     * @param site2 The second target qsite on the grid.
     * @param time The time at which the instruction is applied.
     * @param step The step number of the instruction (useful to specify order within a time-step).
     * @param q1 The first qubit label within the relevant plaquette (for debugging, if applicable).
     * @param q2 The second qubit label within the relevant plaquette (for debugging, if applicable).
     * @param shape The shape of the relevant plaquette (for debugging, if applicable).
     * @param type The type of the relevant plaquette (for debugging, if applicable).
     */
    explicit HW_Instruction(std::string name, unsigned int site1, unsigned int site2, double time, unsigned int step, char q1, char q2, char shape, char type) : 
        name_(name), site1_(site1), site2_(site2), time_(time), step_(step), q1_(q1), q2_(q2), shape_(shape), type_(type) {};

    /**
     * @brief Constructor for shifting the time of an instruction by a given time offset.
     * @param a The original HW_Instruction.
     * @param time_offset The time offset to apply.
     */
    explicit HW_Instruction(const HW_Instruction& a, double time_offset);

    /**
     * @brief Constructor for shifting the sites targeted by an instruction by a number of rows and columns on the grid.
     * @param a The original HW_Instruction.
     * @param row_offset The row offset to apply.
     * @param col_offset The column offset to apply.
     * @param grid The GridManager object for reference.
     */
    explicit HW_Instruction(const HW_Instruction& a, int row_offset, int col_offset, const GridManager& grid);

    // Accessors
    /**
     * @brief Gets the name of the instruction.
     * @return A const reference to the instruction name.
     */
    const std::string& get_name() const {return name_;}
    /**
     * @brief Gets the first qsite acted upon by the instruction.
     * @return The grid array index of the first qsite.
    */
    unsigned int get_site1() const {return site1_;}
    /**
     * @brief Gets the second qsite acted upon by the instruction.
     * 
     * This is usually set to std::numeric_limits<unsigned int>::max() for instructions acting on one qsite.
     * @return The grid array index of the second qsite.
    */
    unsigned int get_site2() const {return site2_;}
    /**
     * @brief Gets the time at which the instruction is applied.
     * @return The time at which the instruction is applied.
    */
    double get_time() const {return time_;}
    /**
     * @brief Gets the step number (specifying order within a time-step) of the instruction.
     * 
     * Usually set to 0 if order does not need to be specified.
     * @return The step number of the instruction.
    */
    unsigned int get_step() const {return step_;}

    // Accessors for debugging output
    /**
     * @brief Gets the first qubit label within the relevant plaquette (for debugging, if applicable).
     * 
     * Usually set to 'X' if HW_Instruction not added as part of a syndrome extraction circuit.
     * @return The first qubit label as a character.
    */
    char get_q1() const {return q1_;}
    /**
     * @brief Gets the second qubit label within the relevant plaquette (for debugging, if applicable).
     * 
     * Usually set to 'X' if HW_Instruction not added as part of a syndrome extraction circuit and ' ' for instructions acting on one qsite.
     * @return The second qubit label as a character.
    */
    char get_q2() const {return q2_;}
    /**
     * @brief Gets the shape of the relevant plaquette (for debugging, if applicable).
     * 
     * Usually set to 'X' if HW_Instruction not added as part of a syndrome extraction circuit.
     * @return The plaquette shape as a character.
    */
    char get_shape() const {return shape_;}
    /**
     * @brief Gets the type of the relevant plaquette (for debugging, if applicable).
     * 
     * Usually set to 'X' if HW_Instruction not added as part of a syndrome extraction circuit.
     * @return The plaquette type as a character.
    */
    char get_type() const {return type_;}

private:
    std::string name_;
    unsigned int site1_;
    unsigned int site2_;
    double time_;
    unsigned int step_;
    char q1_;
    char q2_;
    char shape_;
    char type_;
};

/**
 * @brief Comparison operator for use in sorting hardware instructions before printing them.
 * @param i1 The first HW_Instruction.
 * @param i2 The second HW_Instruction.
 * @return True if i1 acts at an earlier time than i2 or if they act at the same time and i1 has a smaller step number than i2.
 */
bool operator<(const HW_Instruction& i1, const HW_Instruction& i2);

}

#endif //TISCC_INSTRUCTION_HPP