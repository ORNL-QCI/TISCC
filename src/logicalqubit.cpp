#include <TISCC/logicalqubit.hpp>
#include <unordered_map>
#include <iomanip>
#include <cassert>

#define DEBUG_OUTPUT_INSTRUCTIONS

namespace TISCC 
{
    void LogicalQubit::init_stabilizers(unsigned int dx, unsigned int dz, const GridManager& a) {
        // Construct the 4-qubit plaquettes of the surface code
        for (unsigned int i=1; i<dz; i++) {
            for (unsigned int j=1; j<dx; j++) {
                if ((i+j)%2 == 0) {
                    z_plaquettes.push_back(a.get_plaquette(i, j, 'f', 'Z'));
                }
                else if ((i+j)%2 == 1) {
                    x_plaquettes.push_back(a.get_plaquette(i, j, 'f', 'X'));
                }
            }
        }

        /* Construct the 2-qubit plaquettes of the surface code */
        int gap_tmp = 0;
        // Left edge of surface code
        for (unsigned int i=1; i<dz; i+=2) {
            x_plaquettes.push_back(a.get_plaquette(i, 0, 'w', 'X'));
            if (i==dz-1) {gap_tmp = 1;}
        }
        int gap = gap_tmp; gap_tmp = 0;
        for (unsigned int i=1; i<dx; i+=2) {
            // Top edge of surface code 
            if (i!=dx-1) {
                z_plaquettes.push_back(a.get_plaquette(0, i+1, 'n', 'Z'));
                if (i+1==dx-1) {gap_tmp = 1;}
            }
            // Bottom edge of surface code
            if ((gap==1) && (i==dx-1)) {continue;}
            z_plaquettes.push_back(a.get_plaquette(dz, i+gap, 's', 'Z'));
        }
        gap = gap_tmp;
        // Right edge of surface code
        for (unsigned int i=1+gap; i<dz; i+=2) {
            x_plaquettes.push_back(a.get_plaquette(i, dx, 'e', 'X'));
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
        // for (Plaquette p : z_plaquettes) {
        //     for (unsigned int q : measure_qubits) {
        //         if (p.m == q) {std::cerr << "Found duplicate measure qubit among plaquettes."; abort();}
        //     }
        //     measure_qubits.push_back(p.m);
        // }
        // for (Plaquette p : x_plaquettes) {
        //     for (unsigned int q : measure_qubits) {
        //         if (p.m == q) {std::cerr << "Found duplicate measure qubit among plaquettes."; abort();}
        //     }
        //     measure_qubits.push_back(p.m);
        // }

        // TODO: Create a way to find all plaquettes containing a particular qubit
    }

    LogicalQubit::LogicalQubit(unsigned int dx, unsigned int dz, const GridManager& a) : TI_model() { 
        init_stabilizers(dx, dz, a);
        test_stabilizers(dx, dz);
        //print_stabilizers();
    }

    void LogicalQubit::idle(unsigned int cycles, const GridManager& a) {
        // I/O settings
        int W = 15;
        float time = 0;
        std::cout << std::setprecision(1);
        std::cout << std::setiosflags(std::ios::fixed);
        
        // Loop over surface code cycles
        for (unsigned int cycle=0; cycle < cycles; cycle++) {
            // TODO: Make sure that all plaquettes are re-set

            // Loop over instructions and apply them to the plaquettes
            // TODO: Since most of the code is the same between x and z stabilizers, we should throw it into a separate function
            assert(TI_model.get_Z_circuit_Z_type().size() == TI_model.get_X_circuit_N_type().size());
            for (unsigned int i=0; i<TI_model.get_Z_circuit_Z_type().size(); i++) {
                Instruction instruction1 = TI_model.get_Z_circuit_Z_type()[i];
                Instruction instruction2 = TI_model.get_X_circuit_N_type()[i];
                if (instruction1.get_name() != "Idle") {
                    for (Plaquette& p : z_plaquettes) {
                        // Not all instructions are applied to all types of stabilizers. The conditionality is given by the following logic.
                        if (!(p.get_shape() == 'n' && ((instruction1.get_q1() == 'a') || (instruction1.get_q1() == 'b') || (instruction1.get_q2() == 'a') || (instruction1.get_q2() == 'b')))
                            && !(p.get_shape() == 's' && ((instruction1.get_q1() == 'c') || (instruction1.get_q1() == 'd') || (instruction1.get_q2() == 'c') || (instruction1.get_q2() == 'd')))
                            && !(p.get_shape() == 'e' && ((instruction1.get_q1() == 'b') || (instruction1.get_q1() == 'd') || (instruction1.get_q2() == 'b') || (instruction1.get_q2() == 'd')))
                            && !(p.get_shape() == 'w' && ((instruction1.get_q1() == 'a') || (instruction1.get_q1() == 'c') || (instruction1.get_q2() == 'a') || (instruction1.get_q2() == 'c'))))
                        {
                            std::cout << std::setw(W) << time;
                            std::cout << std::setw(W) << instruction1.get_name();
                            std::cout << std::setw(W) << p.get_qsite(instruction1.get_q1());
                            if (instruction1.get_q2() == ' ') {
                                std::cout << std::setw(W) << ' ';  
                            }
                            else if ((instruction1.get_q2() == 'h') && (instruction1.get_name() == "Move")) {
                                p.move_home(instruction1.get_q1());
                                std::cout << std::setw(W) << p.get_qsite(instruction1.get_q1());
                            }
                            else if (instruction1.get_name() == "Move") {
                                std::cout << std::setw(W) << p.get_qsite(instruction1.get_q2());
                                p.apply_move(instruction1.get_q1(), instruction1.get_q2());
                            }
                            else {
                                std::cout << std::setw(W) << p.get_qsite(instruction1.get_q2());
                            }
                            #ifdef DEBUG_OUTPUT_INSTRUCTIONS
                            std::cout << std::setw(W) << instruction1.get_q1();
                            std::cout << std::setw(W) << instruction1.get_q2();
                            std::cout << std::setw(W) << p.get_shape();
                            std::cout << std::setw(W) << p.get_type();
                            #endif
                            std::cout << std::endl;
                        }
                    }
                }
                if (instruction2.get_name() != "Idle") {
                    for (Plaquette& p : x_plaquettes) {
                        if (!(p.get_shape() == 'n' && ((instruction2.get_q1() == 'a') || (instruction2.get_q1() == 'b') || (instruction2.get_q2() == 'a') || (instruction2.get_q2() == 'b')))
                            && !(p.get_shape() == 's' && ((instruction2.get_q1() == 'c') || (instruction2.get_q1() == 'd') || (instruction2.get_q2() == 'c') || (instruction2.get_q2() == 'd')))
                            && !(p.get_shape() == 'e' && ((instruction2.get_q1() == 'b') || (instruction2.get_q1() == 'd') || (instruction2.get_q2() == 'b') || (instruction2.get_q2() == 'd')))
                            && !(p.get_shape() == 'w' && ((instruction2.get_q1() == 'a') || (instruction2.get_q1() == 'c') || (instruction2.get_q2() == 'a') || (instruction2.get_q2() == 'c')))) 
                        {
                            std::cout << std::setw(W) << time;
                            std::cout << std::setw(W) << instruction2.get_name();
                            std::cout << std::setw(W) << p.get_qsite(instruction2.get_q1());
                            if (instruction2.get_q2() == ' ') {
                                std::cout << std::setw(W) << ' ';  
                            }
                            else if ((instruction2.get_q2() == 'h') && (instruction2.get_name() == "Move")) {
                                p.move_home(instruction2.get_q1());
                                std::cout << std::setw(W) << p.get_qsite(instruction2.get_q1());
                            }
                            else if (instruction2.get_name() == "Move") {
                                std::cout << std::setw(W) << p.get_qsite(instruction2.get_q2());
                                p.apply_move(instruction2.get_q1(), instruction2.get_q2());
                            }
                            else {
                                std::cout << std::setw(W) << p.get_qsite(instruction2.get_q2());
                            }
                            #ifdef DEBUG_OUTPUT_INSTRUCTIONS
                            std::cout << std::setw(W) << instruction2.get_q1();
                            std::cout << std::setw(W) << instruction2.get_q2();
                            std::cout << std::setw(W) << p.get_shape();
                            std::cout << std::setw(W) << p.get_type();
                            #endif
                            std::cout << std::endl;
                        }
                    } 
                }
                assert(instruction1.get_time() == instruction2.get_time());
                time += instruction1.get_time();
            }
        }
    }

}