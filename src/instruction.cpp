#include <TISCC/instruction.hpp>


namespace TISCC 
{
    // Comparison operator for use in sorting hardware instructions before printing them 
    // TODO: modify sorting to ensure that it won't commute non-commuting operations
    bool operator<(const HW_Instruction& i1, const HW_Instruction& i2) {
    return ((i1.get_time() < i2.get_time()) || 
    (i1.get_time() == i2.get_time() &&
    i1.get_step() < i2.get_step()) ||
    (i1.get_time() == i2.get_time() &&
    i1.get_step() < i2.get_step() &&
    i1.get_site1() < i2.get_site1()) ||
    (i1.get_time() == i2.get_time() &&
    i1.get_step() < i2.get_step() &&
    i1.get_site1() == i2.get_site1() &&
    i1.get_site2() < i2.get_site2()));
    }
}