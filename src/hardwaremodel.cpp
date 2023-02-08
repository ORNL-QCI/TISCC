#include <TISCC/hardwaremodel.hpp>

#include <limits>
#include <cassert>
#include <algorithm>


namespace TISCC 
{
    // Initialize hash table to map native TI operations to times (in microseconds)
    void HardwareModel::init_TI_ops() {

        // Single-site ops
        TI_ops["Prepare_Z"] = 10;
        TI_ops["Measure_Z"] = 120;
        TI_ops["X_pi/2"] = 10;
        TI_ops["Y_pi/2"] = 10;
        TI_ops["Z_pi/2"] = 0;
        TI_ops["X_pi/4"] = 10;
        TI_ops["Y_pi/4"] = 10;
        TI_ops["Z_pi/4"] = 0;
        TI_ops["X_-pi/4"] = 10;
        TI_ops["Y_-pi/4"] = 10;
        TI_ops["Z_-pi/4"] = 0;

        // Move is currently per-site. It probably needs to be revised and split into cases (i.e. across junctions)
        TI_ops["Move"] = 1.4;

        // Includes Merge, Cool, Interact, and Split 
        TI_ops["ZZ"] = 2000;

    }

    // Helper function to add the Prepare_Z HW_Instruction to a circuit
    float HardwareModel::add_init(const Plaquette& p, char qubit, float time, unsigned int step, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (p.grid()[p.get_qsite(qubit)] != 'O') {
            std::cerr << "HardwareModel::add_init: Can only Prepare Z at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instruction onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Prepare_Z", p.get_qsite(qubit), uint_max, time, step, qubit, ' ', p.get_shape(), p.get_type()));

        // Return updated time
        return time + TI_ops.at("Prepare_Z");

    }

   // Helper function to add H gate in terms of native TI gates to a circuit
    float HardwareModel::add_H(const Plaquette& p, char qubit, float time, unsigned int step, std::vector<HW_Instruction>& circuit) const {
        
        // Perform validity check
        if (p.grid()[p.get_qsite(qubit)] != 'O') {
            std::cerr << "HardwareModel::add_H: Can only apply Hadamard gates at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Y_-pi/4", p.get_qsite(qubit), uint_max, time, step, qubit, ' ', p.get_shape(), p.get_type()));
        circuit.push_back(HW_Instruction("Z_pi/4", p.get_qsite(qubit), uint_max, time + TI_ops.at("Y_-pi/4"), step, qubit, ' ', p.get_shape(), p.get_type()));

        // Return updated time
        return time + TI_ops.at("Y_-pi/4") + TI_ops.at("Z_pi/4");

    }

    // Helper function to add the Measure_Z HW_Instruction to a circuit
    float HardwareModel::add_meas(const Plaquette& p, char qubit, float time, unsigned int step, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (p.grid()[p.get_qsite(qubit)] != 'O') {
            std::cerr << "HardwareModel::add_meas: Can only apply Measure Z at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Measure_Z", p.get_qsite(qubit), uint_max, time, step, qubit, ' ', p.get_shape(), p.get_type()));

        // Return updated time
        return time + TI_ops.at("Measure_Z");

    }

    // Move the measure qubit to the closest site adjacent to a data qubit
    void HardwareModel::move_along_path(Plaquette& p, unsigned int step, std::vector<HW_Instruction>& circuit, float& time,
        const std::vector<unsigned int>& path) const {

        // Validity check 
        if (p.get_qsite('m') != path[0]) {
            std::cerr << "HardwareModel::move_along_path: current site inconsistent with path given." << std::endl;
            abort();
        }

        // Construct HW_Instructions to Move 'm' along the path while also applying the move_to_site plaquette member function (this ensures validity)
        for (unsigned int i=1; i<path.size(); i++) {
            circuit.push_back(HW_Instruction("Move", path[i-1], path[i], time + (i-1)*TI_ops.at("Move"), step, 'm', ' ', p.get_shape(), p.get_type()));
            p.move_to_site('m', path[i]);
        }

        // Update time 
        time += (path.size()-1)*TI_ops.at("Move");

    }

    // Helper function to add CNOT gate in terms of native TI gates to a circuit
    float HardwareModel::add_CNOT(Plaquette& p, char control, char target, float time, unsigned int step, const GridManager& grid, 
        std::vector<HW_Instruction>& circuit) const {
        
        // This function currently requires one of the qubits to be 'm'
        if (control != 'm' && target != 'm') {
            std::cerr << "HardwareModel::add_CNOT: Currently requires either the control or the target to be 'm'." << std::endl;
            abort();
        }

        // This function currently depends on both control and target beginning at their 'home' sites on the plaquette
        if ((!p.is_home(control)) || (!p.is_home(target))) {
            std::cerr << "HardwareModel::add_CNOT: Currently requires control and target to be at their home site on the plaquette." << std::endl;
            abort();
        }

        // A qubit cannot perform a CNOT gate with itself
        if (control == target) {
            std::cerr << "HardwareModel::add_CNOT: A qubit cannot perform a CNOT gate with itself." << std::endl;
            abort();
        }

        // Obtain data qubit label
        char data_qubit = control;
        if (control == 'm') {
            data_qubit = target;
        }

        // Perform initial rotation
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Y_-pi/4", p.get_qsite(target), uint_max, time, step, target, ' ', p.get_shape(), p.get_type()));
        time += TI_ops.at("Y_-pi/4");

        // Retrieve path to data qubit and apply move_along_path
        std::vector<unsigned int> path = grid.get_path(p.get_qsite('m'), p.get_qsite(data_qubit));
        move_along_path(p, step, circuit, time, path);

        // Apply ZZ operation
        circuit.push_back(HW_Instruction("ZZ", p.get_qsite(data_qubit), p.get_qsite('m'), time, step, control, target, p.get_shape(), p.get_type()));
        time += TI_ops.at("ZZ");

        // Move the measure qubit back to its home base for further operations
        reverse(path.begin(), path.end());
        move_along_path(p, step, circuit, time, path);

        // Before final rotations, ensure control and target are home
        assert(p.is_home(control) && p.is_home(target)); 

        // Perform final rotations
        circuit.push_back(HW_Instruction("Z_-pi/4", p.get_qsite(control), uint_max, time, step, control, ' ', p.get_shape(), p.get_type()));
        circuit.push_back(HW_Instruction("X_-pi/4", p.get_qsite(target), uint_max, time, step, target, ' ', p.get_shape(), p.get_type()));
        time += TI_ops.at("Z_-pi/4") + TI_ops.at("X_-pi/4");
        circuit.push_back(HW_Instruction("Z_-pi/4", p.get_qsite(target), uint_max, time, step, target, ' ', p.get_shape(), p.get_type()));

        // Return updated time
        return time + TI_ops.at("Z_-pi/4");
    }

    HardwareModel::HardwareModel() {
        init_TI_ops();
    }

    void HardwareModel::print_TI_ops() const {
        for (auto const &op_pair : TI_ops) {
            std::cout << op_pair.first << ": " << op_pair.second << std::endl;
        }
    }
}