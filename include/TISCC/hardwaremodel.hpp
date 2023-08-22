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
    const std::map<std::string, double>& get_ops() const {return TI_ops;}

    // Compile gates to hardware operations (one each for plaquette&qubit vs. site)
    double add_init(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const;
    double add_init(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_X(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_X_pi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_X_mpi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;   
    double add_Y(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_Y_pi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_Y_mpi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_Z(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_Z_pi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_Z_mpi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_Z_pi8(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_Z_mpi8(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_H(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const;
    double add_H(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_meas(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const;
    double add_meas(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;
    double add_CNOT(Plaquette& p, char control, char target, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const;

    // Translate all qubits left one column on the grid
    double shift_left(const std::set<unsigned int>& qubits, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) const;

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
    std::map<std::string, double> TI_ops;

    // Initialize hash table to define trapped-ion instruction set and map instructions to time (in microseconds)
    void init_TI_ops();

    // HW circuit helper functions to compile Moves
    void move_along_path_for_CNOT(Plaquette& p, unsigned int step, std::vector<HW_Instruction>& circuit, double& time,
        const std::vector<unsigned int>& path, const GridManager& grid) const;
    double move_along_path_for_shift(unsigned int qsite, std::vector<HW_Instruction>& circuit, double time,
        const std::vector<unsigned int>& path, const GridManager& grid) const;
};

}
#endif //TISCC_HARDWAREMODEL_HPP