#include <TISCC/instruction.hpp>


namespace TISCC 
{
    // Comparison operator for use in sorting hardware instructions before printing them 
    bool operator<(const HW_Instruction& i1, const HW_Instruction& i2) {
    return ((i1.get_time() < i2.get_time()));
    }
}