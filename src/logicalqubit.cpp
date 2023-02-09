#include <TISCC/logicalqubit.hpp>

#include <unordered_map>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <limits>
#include <set>

namespace TISCC 
{
    // Construct stabilizers and update set of occupied sites in grid
    void LogicalQubit::init_stabilizers(unsigned int dx, unsigned int dz, GridManager& grid) {
        
        // Construct the 4-qubit plaquettes of the surface code
        for (unsigned int i=1; i<dz; i++) {
            for (unsigned int j=1; j<dx; j++) {
                if ((i+j)%2 == 0) {
                    z_plaquettes.push_back(grid.get_plaquette(i, j, 'f', 'Z'));
                }
                else if ((i+j)%2 == 1) {
                    x_plaquettes.push_back(grid.get_plaquette(i, j, 'f', 'X'));
                }
            }
        }

        /* Construct the 2-qubit plaquettes of the surface code */
        int gap_tmp = 0;
        // Left edge of surface code
        for (unsigned int i=1; i<dz; i+=2) {
            x_plaquettes.push_back(grid.get_plaquette(i, 0, 'w', 'X'));
            if (i==dz-1) {gap_tmp = 1;}
        }
        int gap = gap_tmp; gap_tmp = 0;
        for (unsigned int i=1; i<dx; i+=2) {
            // Top edge of surface code 
            if (i!=dx-1) {
                z_plaquettes.push_back(grid.get_plaquette(0, i+1, 'n', 'Z'));
                if (i+1==dx-1) {gap_tmp = 1;}
            }
            // Bottom edge of surface code
            if ((gap==1) && (i==dx-1)) {continue;}
            z_plaquettes.push_back(grid.get_plaquette(dz, i+gap, 's', 'Z'));
        }
        gap = gap_tmp;
        // Right edge of surface code
        for (unsigned int i=1+gap; i<dz; i+=2) {
            x_plaquettes.push_back(grid.get_plaquette(i, dx, 'e', 'X'));
        }
    }

    // Set up the circuits that we intend to use
    void LogicalQubit::init_circuits() {

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

    void LogicalQubit::print_stabilizers() {

        // Print all x stabilizers
        for (const Plaquette& p : x_plaquettes) {
            std::cout << p.get_type() << " " << p.get_row() << " " << p.get_col() << " " << p.get_shape() << std::endl;
        }

        // Print all z stabilizers
        for (const Plaquette& p : z_plaquettes) {
            std::cout << p.get_type() << " " << p.get_row() << " " << p.get_col() << " " << p.get_shape() << std::endl;
        }
    }

    // Test stabilizers (not fully implemented)
    void LogicalQubit::test_stabilizers(unsigned int dx, unsigned int dz) {
        unsigned int num_data_qubits = dx*dz;
        unsigned int num_measure_qubits = num_data_qubits-1;

        // Check that number of stabilizers is equal to the number expected
        unsigned int num_plaquettes = z_plaquettes.size() + x_plaquettes.size();
        if (num_plaquettes != num_measure_qubits) {
            std::cerr << "Incorrect number of plaquettes: " << num_plaquettes << " " << num_measure_qubits << std::endl; 
            abort();
        }

        // Check that measure qubits are all unique
        std::vector<unsigned int> measure_qubits;
        for (Plaquette p : z_plaquettes) {
            for (unsigned int q : measure_qubits) {
                if (p.get_qsite('m') == q) {std::cerr << "Found duplicate measure qubit among plaquettes."; abort();}
            }
            measure_qubits.push_back(p.get_qsite('m'));
        }
        for (Plaquette p : x_plaquettes) {
            for (unsigned int q : measure_qubits) {
                if (p.get_qsite('m') == q) {std::cerr << "Found duplicate measure qubit among plaquettes."; abort();}
            }
            measure_qubits.push_back(p.get_qsite('m'));
        }

        // TODO: Create a way to find all plaquettes containing a particular qubit
    }

    LogicalQubit::LogicalQubit(unsigned int dx, unsigned int dz, GridManager& grid) : TI_model() { 
        init_stabilizers(dx, dz, grid);
        init_circuits();
        test_stabilizers(dx, dz);
    }

    // Apply a given ``qubit-level'' instruction to all plaquettes in a given vector and add corresponding HW_Instructions to hw_master
    float LogicalQubit::apply_instruction(const Instruction& instr, std::vector<Plaquette>& plaquettes, float time, unsigned int step, 
        const GridManager& grid, std::vector<HW_Instruction>& hw_master) {

        // Create tmp variable for time
        float time_tmp = 0;

        // Don't explicitly apply an "Idle" operation
        if (instr.get_name() != "Idle") {

            // Loop over plaquettes
            for (Plaquette& p : plaquettes) {

                // Certain stabilizer shapes do not accomodate certain qubit labels
                if (p.is_instr_valid(instr)) 
                {
                    
                    // Add single-qubit instructions
                    if ((instr.get_name() == "Prepare_Z") && (instr.get_q2() == ' ')) {
                        time_tmp = TI_model.add_init(p, instr.get_q1(), time, step, hw_master);
                    }

                    else if ((instr.get_name() == "Hadamard") && (instr.get_q2() == ' ')) {
                        time_tmp = TI_model.add_H(p, instr.get_q1(), time, step, hw_master);
                    }

                    else if ((instr.get_name() == "Measure_Z") && (instr.get_q2() == ' ')) {
                        time_tmp = TI_model.add_meas(p, instr.get_q1(), time, step, hw_master);
                    }

                    // Add CNOT gate
                    else if (instr.get_name() == "CNOT") {
                        time_tmp = TI_model.add_CNOT(p, instr.get_q1(), instr.get_q2(), time, step, grid, hw_master);
                    }

                    else {std::cerr << "LogicalQubit::apply_instruction: Invalid instruction given." << std::endl; abort();}
                    
                }
            }
        }

        return time_tmp;
    }

    // Function to return all qsites occupied by the surface code
    std::set<unsigned int> LogicalQubit::occupied_sites() {
        std::set<unsigned int> sites;
        for (char qubit : {'a', 'b', 'c', 'd', 'm'}) {
            for (const Plaquette& p : z_plaquettes) {
                sites.insert(p.get_qsite(qubit));
            }
            for (const Plaquette& p : x_plaquettes) {
                sites.insert(p.get_qsite(qubit));
            }
        }
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();  
        sites.erase(uint_max);
        return sites;
    }

    void LogicalQubit::idle(unsigned int cycles, const GridManager& grid, std::vector<HW_Instruction>& hw_master) {

        // Initialize time counter
        float time = 0;
        
        // Loop over surface code cycles
        for (unsigned int cycle=0; cycle < cycles; cycle++) {

            // Enforce that the two circuits contain the same number of instructions
            assert(Z_Circuit_Z_Type.size() == X_Circuit_N_Type.size());
            unsigned int num_instructions = Z_Circuit_Z_Type.size();

            // Loop over ``qubit-level'' instructions and apply them to plaquettes while adding HW_Instructions to hw_master
            for (unsigned int i=0; i<num_instructions; i++) {
                float t1 = apply_instruction(Z_Circuit_Z_Type[i], z_plaquettes, time, i, grid, hw_master);
                float t2 = apply_instruction(X_Circuit_N_Type[i], x_plaquettes, time, i, grid, hw_master);

                // Increment time counter
                if (t1 == 0) {time = t2;}
                else if (t2 == 0) {time = t1;}
                else {assert(t1==t2); time = t1;}  

            }
            
        }

        // Sort master list of hardware instructions according to overloaded operator<
        std::stable_sort(hw_master.begin(), hw_master.end());

    }

}