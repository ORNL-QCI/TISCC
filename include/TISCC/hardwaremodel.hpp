#ifndef TISCC_HARDWAREMODEL_HPP
#define TISCC_HARDWAREMODEL_HPP

#include <TISCC/instruction.hpp>
#include <TISCC/gridmanager.hpp>
#include <TISCC/plaquette.hpp>

#include<unordered_map>
#include<string>
#include<vector>

namespace TISCC {

// Need to declare this in advance since the below Class depends on it
class Plaquette;
class GridManager;
class HW_Instruction;

// This class defines a set of native hardware operations and constructs circuits out of them
class HardwareModel {
public:
    // Constructor
    explicit HardwareModel();

    // Accessor functions
    const std::unordered_map<std::string, double>& get_ops() const {return TI_ops;}

    // Compile gates to hardware operations (one each for plaquette&qubit vs. site)
    double add_init(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const;
    double add_init(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_X(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_Z(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_H(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const;
    double add_H(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_meas(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const;
    double add_meas(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_CNOT(Plaquette& p, char control, char target, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    // This is just a placeholder to apply a test gate with
    double add_test(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const;

    // Print TI_ops
    void print_TI_ops() const;

    // Accessors
    double get_trap_width() const {return trap_width;}
    double get_cell_width() const {return cell_width;}

private:
    // Parameters of the hardware model
    double trap_width; // us
    double cell_width; // us

    // Hash table to map trapped-ion instructions to time (in microseconds)
    std::unordered_map<std::string, double> TI_ops;

    // Initialize hash table to define trapped-ion instruction set and map instructions to time (in microseconds)
    void init_TI_ops();

    // HW circuit helper function
    void move_along_path(Plaquette& p, unsigned int step, std::vector<HW_Instruction>& circuit, double& time,
        const std::vector<unsigned int>& path, const GridManager& grid) const;
};

}
#endif //TISCC_HARDWAREMODEL_HPP