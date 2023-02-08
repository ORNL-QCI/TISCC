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

    // Compile gates to hardware operations
    float add_init(const Plaquette& p, char qubit, float time, unsigned int step, std::vector<HW_Instruction>& circuit) const;
    float add_H(const Plaquette& p, char qubit, float time, unsigned int step, std::vector<HW_Instruction>& circuit) const;
    float add_meas(const Plaquette& p, char qubit, float time, unsigned int step, std::vector<HW_Instruction>& circuit) const;
    float add_CNOT(Plaquette& p, char control, char target, float time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    // Print TI_ops
    void print_TI_ops() const;

private:
    // Hash table to map trapped-ion instructions to time (in microseconds)
    std::unordered_map<std::string, float> TI_ops;

    // Initialize hash table to define trapped-ion instruction set and map instructions to time (in microseconds)
    void init_TI_ops();

    // HW circuit helper function
    void move_along_path(Plaquette& p, unsigned int step, std::vector<HW_Instruction>& circuit, float& time,
        const std::vector<unsigned int>& path) const;
};

}
#endif //TISCC_HARDWAREMODEL_HPP