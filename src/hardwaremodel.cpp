#include <TISCC/hardwaremodel.hpp>

#include <limits>


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
        TI_ops["Move"] = 1.4;
        
        // Bundle of Merge, Cool, Interact (exp{-i*pi*ZZ/4}), Split operations that can occur between an 'O' site and any site adjacent to it
        TI_ops["ZZ"] = 2000;

        // Possible ops
        // {"Qubit_At", "Prepare_Z", "Measure_Z", "X_pi/2", "Y_pi/2", "Z_pi/2", "X_pi/4", "Y_pi/4", "Z_pi/4", "X_-pi/4", "Y_-pi/4", "Z_-pi/4", "Move", "ZZ"}
    }

   // Helper function to add the Prepare_Z HW_Instruction to a circuit
    float HardwareModel::add_init(const Plaquette& p, char qubit, float time, unsigned int step, std::vector<HW_Instruction>& circuit) {
        // TODO: Insert a check here to make sure there is an 'O' on the grid at the relevant qsite
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        // circuit.push_back(HW_Instruction("Prepare_Z", p.get_qsite(qubit), uint_max, time, step))
        circuit.push_back(HW_Instruction("Prepare_Z", p.get_qsite(qubit), uint_max, time, step, qubit, ' ', p.get_shape(), p.get_type()));
        return time + TI_ops["Prepare_Z"];
    }

   // Helper function to add H gate in terms of native TI gates to a circuit
    float HardwareModel::add_H(const Plaquette& p, char qubit, float time, unsigned int step, std::vector<HW_Instruction>& circuit) {
        // TODO: Insert a check here to make sure there is an 'O' on the grid at the relevant qsite
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        // circuit.push_back(HW_Instruction("Y_-pi/4", p.get_qsite(qubit), uint_max, time, step));
        // circuit.push_back(HW_Instruction("Z_pi/4", p.get_qsite(qubit), uint_max, time + TI_ops["Y_-pi/4"], step));
        circuit.push_back(HW_Instruction("Y_-pi/4", p.get_qsite(qubit), uint_max, time, step, qubit, ' ', p.get_shape(), p.get_type()));
        circuit.push_back(HW_Instruction("Z_pi/4", p.get_qsite(qubit), uint_max, time + TI_ops["Y_-pi/4"], step, qubit, ' ', p.get_shape(), p.get_type()));
        return time + TI_ops["Y_-pi/4"] + TI_ops["Z_pi/4"];
    }

    // Helper function to add the Measure_Z HW_Instruction to a circuit
    float HardwareModel::add_meas(const Plaquette& p, char qubit, float time, unsigned int step, std::vector<HW_Instruction>& circuit) {
        // TODO: Insert a check here to make sure there is an 'O' on the grid at the relevant qsite
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        // circuit.push_back(HW_Instruction("Measure_Z", p.get_qsite(qubit), uint_max, time, step));
        circuit.push_back(HW_Instruction("Measure_Z", p.get_qsite(qubit), uint_max, time, step, qubit, ' ', p.get_shape(), p.get_type()));
        return time + TI_ops["Measure_Z"];
    }

    // Move the measure qubit to the closest site adjacent to the data qubit
    std::vector<unsigned int> HardwareModel::move_next_to(Plaquette& p, char control, char target, unsigned int step, 
        const GridManager& grid, std::vector<HW_Instruction>& circuit, float& time) {

        if (control == 'm') {
            std::vector<unsigned int> path = grid.get_path(p.get_qsite(control), p.get_qsite(target));                 
            for (unsigned int i=1; i<path.size(); i++) {
                // circuit.push_back(HW_Instruction("Move", path[i-1], path[i], time + (i-1)*TI_ops["Move"], step));
                circuit.push_back(HW_Instruction("Move", path[i-1], path[i], time + (i-1)*TI_ops["Move"], step, control, ' ', p.get_shape(), p.get_type()));
            }
            time += (path.size()-1)*TI_ops["Move"];
            p.move_to_site(control, path.back());
            return path;
        }

        else if (target == 'm') {
            // Note we flip the entries to get_path in this case because we always move 'm'
            std::vector<unsigned int> path = grid.get_path(p.get_qsite(target), p.get_qsite(control));
            for (unsigned int i=1; i<path.size(); i++) {
                // circuit.push_back(HW_Instruction("Move", path[i-1], path[i], time + (i-1)*TI_ops["Move"], step));
                circuit.push_back(HW_Instruction("Move", path[i-1], path[i], time + (i-1)*TI_ops["Move"], step, target, ' ', p.get_shape(), p.get_type()));
            }
            time += (path.size()-1)*TI_ops["Move"];
            p.move_to_site(target, path.back());
            return path;
        }
        else  {
            std::cout << "HardwareModel::move_next_to: CNOT is only implemented where either control or target has label 'm'." << std::endl;
            abort();
        }
    }

    // Helper function to add CNOT gate in terms of native TI gates to a circuit
    float HardwareModel::add_CNOT(Plaquette& p, char control, char target, float time, unsigned int step, const GridManager& grid, 
        std::vector<HW_Instruction>& circuit) {
        // TODO: Insert a check here to make sure there is an 'O' on the grid at the target qsite

        // Dummy qsite for q2 in single-qubit operations
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();

        // Perform initial rotation
        // circuit.push_back(HW_Instruction("Y_-pi/4", p.get_qsite(target), uint_max, time, step));
        circuit.push_back(HW_Instruction("Y_-pi/4", p.get_qsite(target), uint_max, time, step, target, ' ', p.get_shape(), p.get_type()));
        time += TI_ops["Y_-pi/4"];

        // Move the measure qubit to the closest site adjacent to the data qubit (note that time is updated within this function)
        std::vector<unsigned int> path = move_next_to(p, control, target, step, grid, circuit, time);

        // Apply ZZ operation
        // circuit.push_back(HW_Instruction("ZZ", p.get_qsite(control), p.get_qsite(target), time, step));
        circuit.push_back(HW_Instruction("ZZ", p.get_qsite(control), p.get_qsite(target), time, step, control, target, p.get_shape(), p.get_type()));
        time += TI_ops["ZZ"];

        // Move the measure qubit back to its home base for further operations
        p.move_home('m');
        for (unsigned int i=1; i<path.size(); i++) {
            // circuit.push_back(HW_Instruction("Move", path[path.size()-i], path[path.size()-i-1], time + (i-1)*TI_ops["Move"], step));
            circuit.push_back(HW_Instruction("Move", path[path.size()-i], path[path.size()-i-1], time + (i-1)*TI_ops["Move"], step, 'm', ' ', p.get_shape(), p.get_type()));
        }
        time += (path.size()-1)*TI_ops["Move"];

        // Perform final rotations
        // circuit.push_back(HW_Instruction("Z_-pi/4", p.get_qsite(control), uint_max, time, step));
        // circuit.push_back(HW_Instruction("X_-pi/4", p.get_qsite(target), uint_max, time, step));
        circuit.push_back(HW_Instruction("Z_-pi/4", p.get_qsite(control), uint_max, time, step, control, ' ', p.get_shape(), p.get_type()));
        circuit.push_back(HW_Instruction("X_-pi/4", p.get_qsite(target), uint_max, time, step, target, ' ', p.get_shape(), p.get_type()));
        time += TI_ops["Z_-pi/4"] + TI_ops["X_-pi/4"];
        // circuit.push_back(HW_Instruction("Z_-pi/4", p.get_qsite(target), uint_max, time, step));
        circuit.push_back(HW_Instruction("Z_-pi/4", p.get_qsite(target), uint_max, time, step, target, ' ', p.get_shape(), p.get_type()));
        return time + TI_ops["Z_-pi/4"];
    }

    // Set up the circuits that we intend to use
    // TODO: Place this in the LogicalQubit class
    void HardwareModel::init_circuits() {

        // Add instructions to the Z plaquette measurement circuit
        Z_Circuit_Z_Type.push_back(Instruction("Prepare_Z", 'm', ' ')); 
        Z_Circuit_Z_Type.push_back(Instruction("Idle", ' ', ' '));   
        Z_Circuit_Z_Type.push_back(Instruction("CNOT", 'a', 'm'));
        Z_Circuit_Z_Type.push_back(Instruction("CNOT", 'b', 'm'));
        Z_Circuit_Z_Type.push_back(Instruction("CNOT", 'c', 'm'));
        Z_Circuit_Z_Type.push_back(Instruction("CNOT", 'd', 'm'));
        Z_Circuit_Z_Type.push_back(Instruction("Idle", ' ', ' ')); 
        Z_Circuit_Z_Type.push_back(Instruction("Measure_Z", 'm', ' '));

        // Add instructions to the X plaquette measurement circuit
        X_Circuit_N_Type.push_back(Instruction("Prepare_Z", 'm', ' '));  
        X_Circuit_N_Type.push_back(Instruction("Hadamard", 'm', ' '));   
        X_Circuit_N_Type.push_back(Instruction("CNOT", 'm', 'a'));
        X_Circuit_N_Type.push_back(Instruction("CNOT", 'm', 'c'));
        X_Circuit_N_Type.push_back(Instruction("CNOT", 'm', 'b'));
        X_Circuit_N_Type.push_back(Instruction("CNOT", 'm', 'd'));
        X_Circuit_N_Type.push_back(Instruction("Hadamard", 'm', ' '));   
        X_Circuit_N_Type.push_back(Instruction("Measure_Z", 'm', ' '));
    }

    HardwareModel::HardwareModel() {
        init_TI_ops();
        init_circuits();
    }
}