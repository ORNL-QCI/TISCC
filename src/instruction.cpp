#include <TISCC/instruction.hpp>

namespace TISCC 
{
    // This constructor is for shifting the time of an instruction by a given time offset
    HW_Instruction::HW_Instruction(const HW_Instruction& a, double time_offset) : name_(a.get_name()), site1_(a.get_site1()), site2_(a.get_site2()), time_(a.get_time() + time_offset),
        step_(a.get_step()), q1_(a.get_q1()), q2_(a.get_q2()), shape_(a.get_shape()), type_(a.get_type()) {};

    // This constructor is for shifting the sites targeted by all instructions by a number of rows and columns on the grid
    HW_Instruction::HW_Instruction(const HW_Instruction& a, int row_offset, int col_offset, const GridManager& grid) : name_(a.get_name()), 
        site1_(grid.shift_qsite(a.get_site1(), row_offset, col_offset)), site2_(grid.shift_qsite(a.get_site2(), row_offset, col_offset)), 
        time_(a.get_time()), step_(a.get_step()), q1_(a.get_q1()), q2_(a.get_q2()), shape_(a.get_shape()), type_(a.get_type()) {};

    // Once operations are compiled into hardware instructions, this function is used to output them
    std::optional<std::vector<unsigned int>> HW_Instruction::print_hw_master(std::ostream& output, const std::vector<HW_Instruction>& hw_master, const std::set<unsigned int>& occupied_sites, bool debug, bool stim) 
    {
        // I/O settings
        int W = 15;
        output << std::setprecision(2);
        output << std::setiosflags(std::ios::fixed);

        // Initialize HardwareModel (TODO: have a HardwareModel be input)
        HardwareModel TI_model;

        // Process occupied qsites
        std::map<unsigned int, unsigned int> qsites_to_qubits;
        std::vector<unsigned int> meas_idx_to_qsite;
        if (!stim)
        {
            // For typical TISCC circuits, dump qsites to front of circuit
            for (unsigned int site : occupied_sites) {
                output << std::setw(W) << -1.0;
                output << std::setw(W) << "Qubit_at";
                output << std::setw(W) << site;
                output << std::endl;
            }
        }

        else
        {
            // For Stim circuits, we associate qubit indices to occupied sites 
            unsigned int counter = 0;
            for (unsigned int site : occupied_sites) {
                qsites_to_qubits[site] = counter;
                counter++;
            }
        }

        // Output HW instructions to file
        double current_time = 0;
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();  
        for (const HW_Instruction& instruction : hw_master) {

            if (!stim)
            {
                output << std::setw(W) << instruction.get_time();
                output << std::setw(W) << instruction.get_name();
                if (instruction.get_site2() != uint_max) {
                    output << std::setw(W) << ((std::to_string(instruction.get_site1()) += ",") += std::to_string(instruction.get_site2()));
                }
                else {
                    output << std::setw(W) << instruction.get_site1();
                }
                if (debug) {
                    output << std::setw(W) << instruction.get_step();
                    output << std::setw(W) << instruction.get_q1();
                    output << std::setw(W) << instruction.get_q2();
                    output << std::setw(W) << instruction.get_shape();
                    output << std::setw(W) << instruction.get_type();
                }
                output << std::endl;
            }

            else
            {
                // Check if time has changed and put a tick if so
                if (instruction.get_time() != current_time)
                {
                    double time_diff = instruction.get_time() - current_time;
                    current_time = instruction.get_time();
                    output << "TICK[" << time_diff << "]" << std::endl;
                }

                // Move instructions update the qsites_to_qubits map
                if (instruction.get_name() == "Move")
                {
                    qsites_to_qubits[instruction.get_site2()] = qsites_to_qubits[instruction.get_site1()];
                    qsites_to_qubits.erase(instruction.get_site1());
                }

                /* We handle this by gate names instead of as above because we need to discriminate between Clifford and non-Clifford gates */
                // Handle single-qubit operations
                else if ((instruction.get_name() == "Prepare_Z") ||
                         (instruction.get_name() == "Measure_Z") ||
                         (instruction.get_name() == "X_pi/2")    || 
                         (instruction.get_name() == "Y_pi/2")    || 
                         (instruction.get_name() == "Z_pi/2")    || 
                         (instruction.get_name() == "X_pi/4")    ||
                         (instruction.get_name() == "Y_pi/4")    ||  
                         (instruction.get_name() == "Z_pi/4")    || 
                         (instruction.get_name() == "X_-pi/4")   ||
                         (instruction.get_name() == "Y_-pi/4")   || 
                         (instruction.get_name() == "Z_-pi/4"))
                {
                    if (instruction.get_name() == "Measure_Z")
                    {
                        meas_idx_to_qsite.push_back(instruction.get_site1());
                    }
                    output << TI_model.get_TI_ops_to_stim().at(instruction.get_name()) << "[" << instruction.get_time() << "] " << qsites_to_qubits[instruction.get_site1()] << std::endl;
                }

                // Handle two-qubit operations
                else if (instruction.get_name() == "ZZ")
                {
                    output << TI_model.get_TI_ops_to_stim().at(instruction.get_name())  << "[" << instruction.get_time() << "] " << qsites_to_qubits[instruction.get_site1()] << " " << qsites_to_qubits[instruction.get_site2()] << std::endl;
                }

                else if ((instruction.get_name() == "Z_pi/8") || (instruction.get_name() == "Z_-pi/8"))
                {
                    throw std::invalid_argument("Non-Clifford gate cannot be printed in Stim circuit: " + instruction.get_name());
                }

                else
                {
                    throw std::invalid_argument("Instruction not recognized: " + instruction.get_name());
                }

            }
        }

        return stim ? std::optional<std::vector<unsigned int>>{meas_idx_to_qsite} : std::nullopt;
    }

    // Comparison operator for use in sorting hardware instructions before printing them 
    bool operator<(const HW_Instruction& i1, const HW_Instruction& i2) {
        return (((i1.get_time() < i2.get_time())) ||
                ((i1.get_time() == i2.get_time()) && i1.get_step() < i2.get_step()));
    }
}

/* 
std::stable_sort should keep elements in order UNLESS:
1. i1 precedes i2 in time
2. i1 equals i2 in time and i1 precedes i2 in circuit steps

Assumptions:
1. The time of a later step is always greater than or equal to the time of a previous step

*/