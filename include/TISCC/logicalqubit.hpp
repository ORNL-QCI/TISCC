#ifndef TISCC_LOGICALQUBIT_HPP
#define TISCC_LOGICALQUBIT_HPP

#include <TISCC/hardwaremodel.hpp>
#include <TISCC/plaquette.hpp>
#include <TISCC/gridmanager.hpp>
#include <TISCC/instruction.hpp>

#include <vector>
#include <unordered_map>

namespace TISCC {

// Consists of the stabilizer plaquettes of a surface code patch and defined operations on them
class LogicalQubit {
public:
    // Constructor
    explicit LogicalQubit(unsigned int dx, unsigned int dz, const GridManager& grid);

    // Operations
    void idle(unsigned int cycles);

    // 
    void print_stabilizers();

private:
    // Vectors of X or Z plaquettes 
    std::vector<Plaquette> x_plaquettes;
    std::vector<Plaquette> z_plaquettes;

    // Contains details of hardware native gates and stabilizer circuits
    HardwareModel TI_model;

    // Define stabilizers
    void init_stabilizers(unsigned int dx, unsigned int dz, const GridManager& grid); 

    // Test stabilizers
    void test_stabilizers(unsigned int dx, unsigned int dz);

    // Check to see if a given instruction is valid on a given plaquette 
    bool is_instr_valid(const Instruction& instr, const Plaquette& p);

    // Apply a given instruction to all plaquettes in a given vector
    void apply_instruction(const Instruction& instr, std::vector<Plaquette>& plaquettes, float time, unsigned int step,
        std::vector<HW_Instruction>& idle_operation);
};
}

#endif //TISCC_LOGICALQUBIT_HPP