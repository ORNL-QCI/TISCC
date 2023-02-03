#ifndef TISCC_LOGICALQUBIT_HPP
#define TISCC_LOGICALQUBIT_HPP

#include <TISCC/hardwaremodel.hpp>
#include <TISCC/plaquette.hpp>
#include <TISCC/gridmanager.hpp>
#include <TISCC/instruction.hpp>

#include <vector>
#include <unordered_map>

namespace TISCC {

class LogicalQubit {
public:
    explicit LogicalQubit(unsigned int dx, unsigned int dz, const GridManager& a);
    void idle(unsigned int cycles, const GridManager& a);
    void print_stabilizers();

private:
    // Vectors of X or Z plaquettes 
    std::vector<Plaquette> x_plaquettes;
    std::vector<Plaquette> z_plaquettes;

    // Contains details of hardware native gates and stabilizer circuits
    HardwareModel TI_model;

    // Define stabilizers
    void init_stabilizers(unsigned int dx, unsigned int dz, const GridManager& a); 

    // Test stabilizers
    void test_stabilizers(unsigned int dx, unsigned int dz);

    // std::vector<unsigned int> data_qubits;
    // std::vector<unsigned int> measure_qubits;
};
}

#endif //TISCC_LOGICALQUBIT_HPP