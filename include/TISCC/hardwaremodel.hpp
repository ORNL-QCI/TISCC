#ifndef TISCC_HARDWAREMODEL_HPP
#define TISCC_HARDWAREMODEL_HPP

#include <TISCC/instruction.hpp>
#include <TISCC/gridmanager.hpp>
#include <TISCC/plaquette.hpp>

#include<unordered_map>
#include<string>
#include<vector>

namespace TISCC {

// This class defines a set of native hardware operations and constructs circuits out of them
class HardwareModel {
public:
    // Constructor
    explicit HardwareModel();

    // Accessor functions
    const std::unordered_map<std::string, float>& get_ops() const {return TI_ops;}
    const std::vector<Instruction>& get_Z_circuit_Z_type() const {return Z_Circuit_Z_Type;}
    const std::vector<Instruction>& get_X_circuit_N_type() const {return X_Circuit_N_Type;}

    // Compile gates to hardware operations
    float add_init(const Plaquette& p, char qubit, float time, unsigned int step, std::vector<HW_Instruction>& circuit);
    float add_H(const Plaquette& p, char qubit, float time, unsigned int step, std::vector<HW_Instruction>& circuit);
    float add_meas(const Plaquette& p, char qubit, float time, unsigned int step, std::vector<HW_Instruction>& circuit);
    float add_CNOT(Plaquette& p, char control, char target, float time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit);


private:
    // Hash table to map trapped-ion instructions to time (in microseconds)
    std::unordered_map<std::string, float> TI_ops;

    // Vectors of instructions to hold syndrome measurement circuits
    std::vector<Instruction> Z_Circuit_Z_Type;
    std::vector<Instruction> X_Circuit_N_Type;

    // Initialize hash table to define trapped-ion instruction set and map instructions to time (in microseconds)
    void init_TI_ops();

    // Set up the circuits that we intend to use
    void init_circuits();

    // HW circuit helper function
    std::vector<unsigned int> move_next_to(Plaquette& p, char control, char target, unsigned int step, const GridManager& grid, 
        std::vector<HW_Instruction>& circuit, float& time);
};

}
#endif //TISCC_HARDWAREMODEL_HPP