#ifndef TISCC_HARDWAREMODEL_HPP
#define TISCC_HARDWAREMODEL_HPP

#include <TISCC/instruction.hpp>

#include<unordered_map>
#include<string>
#include<vector>

namespace TISCC {

// This class stores indices to grid_ that define the site locations of qubits in a surface code plaquette
class HardwareModel {
public:
    // Constructor
    explicit HardwareModel();

    // Accessor functions
    const std::unordered_map<std::string, float>& get_ops() const {return TI_ops;}
    const std::vector<Instruction>& get_Z_circuit_Z_type() const {return Z_Circuit_Z_Type;}
    const std::vector<Instruction>& get_X_circuit_N_type() const {return X_Circuit_N_Type;}

private:
    // Hash table to map trapped-ion instructions to time (in microseconds)
    std::unordered_map<std::string, float> TI_ops;

    // Vectors of instructions to hold stabilizer circuits
    std::vector<Instruction> Z_Circuit_Z_Type;
    std::vector<Instruction> X_Circuit_N_Type;

    // Set up the circuits that we intend to use
    void add_H(char qubit, std::vector<Instruction>& circuit);
    void idle_while_H(std::vector<Instruction>& circuit);
    void add_CNOT(char control, char target, std::vector<Instruction>& circuit);
    void init_circuits();

    // Initialize hash table to define trapped-ion instruction set and map instructions to time (in microseconds)
    void init_TI_ops();
};

}
#endif //TISCC_HARDWAREMODEL_HPP