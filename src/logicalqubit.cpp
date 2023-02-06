#include <TISCC/logicalqubit.hpp>

#include <unordered_map>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <limits>

#define DEBUG_OUTPUT_INSTRUCTIONS

namespace TISCC 
{
    void LogicalQubit::init_stabilizers(unsigned int dx, unsigned int dz, const GridManager& grid) {
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

    LogicalQubit::LogicalQubit(unsigned int dx, unsigned int dz, const GridManager& grid) : TI_model() { 
        init_stabilizers(dx, dz, grid);
        test_stabilizers(dx, dz);
        //print_stabilizers();
    }

   // Check to see if a given instruction is valid on a given plaquette 
    bool LogicalQubit::is_instr_valid(const Instruction& instr, const Plaquette& p) {
        bool n_valid = !(p.get_shape() == 'n' && ((instr.get_q1() == 'a') || (instr.get_q1() == 'b') || (instr.get_q2() == 'a') || (instr.get_q2() == 'b')));
        bool s_valid = !(p.get_shape() == 's' && ((instr.get_q1() == 'c') || (instr.get_q1() == 'd') || (instr.get_q2() == 'c') || (instr.get_q2() == 'd')));
        bool e_valid = !(p.get_shape() == 'e' && ((instr.get_q1() == 'b') || (instr.get_q1() == 'd') || (instr.get_q2() == 'b') || (instr.get_q2() == 'd')));
        bool w_valid = !(p.get_shape() == 'w' && ((instr.get_q1() == 'a') || (instr.get_q1() == 'c') || (instr.get_q2() == 'a') || (instr.get_q2() == 'c')));
        bool return_val = n_valid && s_valid && e_valid && w_valid;
        return return_val;
    }

    // Apply a given instruction to all plaquettes in a given vector
    void LogicalQubit::apply_instruction(const Instruction& instr, std::vector<Plaquette>& plaquettes, float time, unsigned int step, 
        std::vector<HW_Instruction>& idle_operation) {

        // Use the maximum possible unsigned int to designate dummy qsites
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();

        // Don't explicitly apply an "Idle" operation
        if (instr.get_name() != "Idle") {

            // Loop over plaquettes
            for (Plaquette& p : plaquettes) {

                // Implementing multi-step move operation
                // if (instr.get_name() == "CNOT") {
                //     TI_model.CNOT(instr.get_q1(), instr.get_q2(), idle_operation);
                // }

                // Certain stabilizer shapes do not accomodate certain qubit labels
                if (is_instr_valid(instr, p)) 
                {
                    // Apply ops to plaquettes and create hardware instruction
                    // TODO: Add logic to enforce grid-based conditions on operations i.e. that 'O' is necessary for most ops.
                    if (instr.get_q2() == ' ') { 
                        HW_Instruction instruction = HW_Instruction(instr.get_name(), p.get_qsite(instr.get_q1()), uint_max, time, step);
                        idle_operation.push_back(instruction);
                    }
                    // The 'h' designation tells a qubit to go to its home on the grid
                    else if ((instr.get_q2() == 'h') && (instr.get_name() == "Move")) {
                        unsigned int old_site = p.get_qsite(instr.get_q1()); 
                        p.move_home(instr.get_q1());
                        HW_Instruction instruction = HW_Instruction(instr.get_name(), old_site, p.get_qsite(instr.get_q1()), time, step);
                        idle_operation.push_back(instruction);
                    }
                    // Otherwise a move operation will target another qubit's qsite
                    else if (instr.get_name() == "Move") {
                        unsigned int old_site = p.get_qsite(instr.get_q1());
                        p.apply_move(instr.get_q1(), instr.get_q2());
                        HW_Instruction instruction = HW_Instruction(instr.get_name(), old_site, p.get_qsite(instr.get_q1()), time, step);
                        idle_operation.push_back(instruction);
                    }
                    // The only other case currently is ZZ, which targets an adjacent site
                    else if (instr.get_name() == "ZZ") {
                        HW_Instruction instruction = HW_Instruction(instr.get_name(), p.get_qsite(instr.get_q1()), uint_max, time, step);
                        idle_operation.push_back(instruction);                        
                    }
                    else {std::cerr << "LogicalQubit::apply_instruction: Invalid instruction given." << std::endl; abort();}
                    
                    // #ifdef DEBUG_OUTPUT_INSTRUCTIONS
                    // std::cout << std::setw(W) << instr.get_q1();
                    // std::cout << std::setw(W) << instr.get_q2();
                    // std::cout << std::setw(W) << p.get_shape();
                    // std::cout << std::setw(W) << p.get_type();
                    // #endif
                }
            }
        }
    }

    void LogicalQubit::idle(unsigned int cycles) {
        
        // I/O settings
        int W = 15;
        std::cout << std::setprecision(1);
        std::cout << std::setiosflags(std::ios::fixed);

        // Initialize master list of hardware instructions
        std::vector<HW_Instruction> idle_operation;

        // Initialize time counter
        float time = 0;
        
        // Loop over surface code cycles
        for (unsigned int cycle=0; cycle < cycles; cycle++) {
            // TODO: Make sure that all plaquettes are re-set

            // To ensure synchronicity of the two different types of circuits we employ, we enforce that they contain the same number of instructions
            // (and later that parallel instructions take the same amount of time)
            assert(TI_model.get_Z_circuit_Z_type().size() == TI_model.get_X_circuit_N_type().size());
            unsigned int num_instructions = TI_model.get_Z_circuit_Z_type().size();

            // Loop over instructions and apply them to plaquettes
            for (unsigned int i=0; i<num_instructions; i++) {
                assert(TI_model.get_Z_circuit_Z_type()[i].get_time() == TI_model.get_X_circuit_N_type()[i].get_time());
                apply_instruction(TI_model.get_Z_circuit_Z_type()[i], z_plaquettes, time, i, idle_operation);
                apply_instruction(TI_model.get_X_circuit_N_type()[i], x_plaquettes, time, i, idle_operation);

                // Increment time counter
                time += TI_model.get_Z_circuit_Z_type()[i].get_time();
            }
        }

        // Sort master list of hardware instructions according to overloaded operator<
        std::sort(idle_operation.begin(), idle_operation.end());

        // Output HW instructions to file
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();  
        for (const HW_Instruction& instruction : idle_operation) {
            std::cout << std::setw(W) << instruction.get_time();
            std::cout << std::setw(W) << instruction.get_name();
            std::cout << std::setw(W) << instruction.get_site1();
            if (instruction.get_site2() != uint_max) {std::cout << std::setw(W) << instruction.get_site2();}
            else {std::cout << std::setw(W) << " ";}
            std::cout << std::endl;
        }

    }

}