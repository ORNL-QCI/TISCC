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
    void HW_Instruction::print_hw_master(std::ostream& output, const std::vector<HW_Instruction>& hw_master, const std::set<unsigned int>& occupied_sites, bool debug) 
    {
        // I/O settings
        int W = 15;
        output << std::setprecision(1);
        output << std::setiosflags(std::ios::fixed);

        // Dump qsites 
        for (unsigned int site : occupied_sites) {
            output << std::setw(W) << -1.0;
            output << std::setw(W) << "Qubit_at";
            output << std::setw(W) << site;
            output << std::endl;
        }

        // Output HW instructions to file
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();  
        for (const HW_Instruction& instruction : hw_master) {
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