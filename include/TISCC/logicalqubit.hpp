#ifndef TISCC_LOGICALQUBIT_HPP
#define TISCC_LOGICALQUBIT_HPP

#include <TISCC/hardwaremodel.hpp>
#include <TISCC/plaquette.hpp>
#include <TISCC/gridmanager.hpp>
#include <TISCC/instruction.hpp>

#include <vector>
#include <set>
#include <unordered_map>

namespace TISCC {

// Consists of the stabilizer plaquettes of a surface code patch and defined operations on them
class LogicalQubit {
public:
    // Constructor
    explicit LogicalQubit(unsigned int dx, unsigned int dz, GridManager& grid);

    // Operations
    float idle(unsigned int cycles, const GridManager& grid, std::vector<HW_Instruction>& hw_master, float time);
    float prepz(const GridManager& grid, std::vector<HW_Instruction>& hw_master, float time);
    float prepx(const GridManager& grid, std::vector<HW_Instruction>& hw_master, float time);
    float measz(const GridManager& grid, std::vector<HW_Instruction>& hw_master, float time);
    float measx(const GridManager& grid, std::vector<HW_Instruction>& hw_master, float time);

    // Print functions
    void print_stabilizers();

    // Function to output all qsites occupied by the surface code
    std::set<unsigned int> occupied_sites();

    // Function to output all qsites occupied by data qubits on the surface code
    std::set<unsigned int> data_qsites();

private:
    // Vectors of X or Z plaquettes 
    std::vector<Plaquette> x_plaquettes;
    std::vector<Plaquette> z_plaquettes;

    // Vectors of instructions to hold syndrome measurement circuits
    std::vector<Instruction> Z_Circuit_Z_Type;
    std::vector<Instruction> X_Circuit_N_Type;

    // Contains details of hardware native gates and stabilizer circuits
    HardwareModel TI_model;

    // Construct stabilizers and update set of occupied sites in grid
    void init_stabilizers(unsigned int dx, unsigned int dz, GridManager& grid); 

    // Test stabilizers (not fully implemented)
    void test_stabilizers(unsigned int dx, unsigned int dz);

    // Set up the circuits that we intend to use
    void init_circuits();

    // Check to see if a given instruction is valid on a given plaquette 
    bool is_instr_valid(const Instruction& instr, const Plaquette& p);

    // Apply a given instruction to all plaquettes in a given vector
    float apply_instruction(const Instruction& instr, std::vector<Plaquette>& plaquettes, float time, unsigned int step,
        const GridManager& grid, std::vector<HW_Instruction>& idle_operation);

};
}

#endif //TISCC_LOGICALQUBIT_HPP