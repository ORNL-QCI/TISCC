#include <TISCC/instruction.hpp>

namespace TISCC 
{
    // Comparison operator for use in sorting hardware instructions before printing them 
    bool operator<(const HW_Instruction& i1, const HW_Instruction& i2) {
        return ((i1.get_time() < i2.get_time()));
    }
}

/* 
std::stable_sort should keep elements in order UNLESS:
1. i1 precedes i2 in time

Assumptions:
1. Hardware ops were added "step by step" through the original circuit
2. The time of a later step is always greater than or equal to the time of a previous step

*/