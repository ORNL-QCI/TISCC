#include <TISCC/hardwaremodel.hpp>


namespace TISCC 
{
    // Initialize hash table to map native TI operations to times (in microseconds)
    void HardwareModel::init_TI_ops() {
        TI_ops["Initialize"] = 10;
        TI_ops["Measure"] = 120;
        // Move is eight times the move time of 1.4us, since 8 is the longest distance traveled by our circuits
        TI_ops["Move"] = 11.2;
        TI_ops["X_pi/2"] = 10;
        TI_ops["Y_pi/2"] = 10;
        TI_ops["Z_pi/2"] = 0;
        TI_ops["X_pi/4"] = 10;
        TI_ops["Y_pi/4"] = 10;
        TI_ops["Z_pi/4"] = 0;
        TI_ops["X_-pi/4"] = 10;
        TI_ops["Y_-pi/4"] = 10;
        TI_ops["Z_-pi/4"] = 0;
        // Bundle of Merge, Cool, Interact, Split
        TI_ops["ZZ"] = 2000;
    }

   // Helper function to add H gate in terms of native TI gates to a circuit
    void HardwareModel::add_H(char qubit, std::vector<Instruction>& circuit) {
        circuit.push_back(Instruction("Y_-pi/4", qubit, ' ', TI_ops["Y_-pi/4"]));
        circuit.push_back(Instruction("Z_pi/4", qubit, ' ', TI_ops["Z_pi/4"]));
    }

    // Helper function to add the idle ops which complement an H gate (for use when waiting on another circuit)
    void HardwareModel::idle_while_H(std::vector<Instruction>& circuit) {
        circuit.push_back(Instruction("Idle", ' ', ' ', TI_ops["Y_-pi/4"]));
        circuit.push_back(Instruction("Idle", ' ', ' ', TI_ops["Z_pi/4"]));
    }

    // Helper function to add CNOT gate in terms of native TI gates to a circuit
    void HardwareModel::add_CNOT(char control, char target, std::vector<Instruction>& circuit) {
        circuit.push_back(Instruction("Y_-pi/4", target, ' ', TI_ops["Y_-pi/4"]));
        circuit.push_back(Instruction("ZZ", control, target, TI_ops["ZZ"]));
        circuit.push_back(Instruction("Z_-pi/4", control, ' ', TI_ops["Z_-pi/4"]));
        circuit.push_back(Instruction("X_-pi/4", target, ' ', TI_ops["X_-pi/4"]));
        circuit.push_back(Instruction("Z_-pi/4", target, ' ', TI_ops["Z_-pi/4"]));
    }

    // Set up the circuits that we intend to use
        /* TODO:
            - Fix error wherein single-qubit rotations on 'm' will be applied to a plaquette 
                even if the overarching CNOT does not apply to its shape
            - Break Move operations down into the corresponding 1-site moves
                - NOTE that right now Move destinations are always qubits not sites
                - Should include them WITHIN the CNOT function above
                    - Apply single-qubit gates with qubits at their homes
                    - Move target to site adjacent to control
                    - Apply ZZ gate to these two sites (there is an implied merge, cool, interact, split)
                    - Move target 'home' and then apply remaining single-qubit gates
                    - Note that two of them can be applied in parallel
        */ 
    void HardwareModel::init_circuits() {

        // Add instructions to the Z plaquette measurement circuit
        Z_Circuit_Z_Type.push_back(Instruction("Initialize", 'm', ' ', TI_ops["Initialize"]));  
        idle_while_H(Z_Circuit_Z_Type);
        Z_Circuit_Z_Type.push_back(Instruction("Move", 'm', 'a', TI_ops["Move"]));
        add_CNOT('a', 'm', Z_Circuit_Z_Type);
        Z_Circuit_Z_Type.push_back(Instruction("Move", 'm', 'b', TI_ops["Move"]));
        add_CNOT('b', 'm', Z_Circuit_Z_Type); 
        Z_Circuit_Z_Type.push_back(Instruction("Move", 'm', 'c', TI_ops["Move"]));
        add_CNOT('c', 'm', Z_Circuit_Z_Type);
        Z_Circuit_Z_Type.push_back(Instruction("Move", 'm', 'd', TI_ops["Move"]));
        add_CNOT('d', 'm', Z_Circuit_Z_Type);
        Z_Circuit_Z_Type.push_back(Instruction("Move", 'm', 'h', TI_ops["Move"]));
        idle_while_H(Z_Circuit_Z_Type);
        Z_Circuit_Z_Type.push_back(Instruction("Measure", 'm', ' ', TI_ops["Measure"]));

        // Add instructions to the X plaquette measurement circuit
        X_Circuit_N_Type.push_back(Instruction("Initialize", 'm', ' ', TI_ops["Initialize"]));
        add_H('m', X_Circuit_N_Type);
        X_Circuit_N_Type.push_back(Instruction("Move", 'm', 'a', TI_ops["Move"]));
        add_CNOT('m', 'a', X_Circuit_N_Type);
        X_Circuit_N_Type.push_back(Instruction("Move", 'm', 'c', TI_ops["Move"]));
        add_CNOT('m', 'c', X_Circuit_N_Type);
        X_Circuit_N_Type.push_back(Instruction("Move", 'm', 'b', TI_ops["Move"]));
        add_CNOT('m', 'b', X_Circuit_N_Type);
        X_Circuit_N_Type.push_back(Instruction("Move", 'm', 'd', TI_ops["Move"]));
        add_CNOT('m', 'd', X_Circuit_N_Type);
        X_Circuit_N_Type.push_back(Instruction("Move", 'm', 'h', TI_ops["Move"]));
        add_H('m', X_Circuit_N_Type);
        X_Circuit_N_Type.push_back(Instruction("Measure", 'm', ' ', TI_ops["Measure"])); 
    }

    HardwareModel::HardwareModel() {
        init_TI_ops();
        init_circuits();
    }
}